// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Yijia Chen $
// $Authors: Yijia Chen $
// --------------------------------------------------------------------------

// This file contains the FeatureFinderIonMobility TOPP tool. It is meant to be used with the IonMobilityBinner TOPP
// tool as part of the ion mobility feature finding workflow. No Doxygen documentation has been added because I was
// unable to build documentation on my local Windows machine.

// This TOPP tool takes <num_bins> mzML files as input (preferably the output files from IonMobilityBinner), runs the
// ion mobility feature finding algorithm on them, and outputs a single featureXML file, holding the resulting
// features.

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include "IonMobilityBinnerCore.hpp"

using namespace OpenMS;

class TOPPFeatureFinderIonMobility : public TOPPBase
{
public:
  TOPPFeatureFinderIonMobility() :
    TOPPBase("FeatureFinderIonMobility", "Detects three-dimensional features in LC-IMS-MS data."),
    debugMode_(false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // The input prefix for the <num_bins> mzML input files. The required format is
    // "<in_prefix>-<pass number>-<bin number>.mzML".
    registerStringOption_("in_prefix", "<string>", "", "input prefix");
    registerIntOption_("num_bins", "<int>", 50, "number of ion mobility bins", false);
    // Specifies whether or not a second binning pass was performed, and if it should be considered when running the
    // feature finder.
    registerFlag_("use_offset", "offset ion mobility bins");

    registerOutputFile_("out", "<file>", "", "output file", false);
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    // Specifies whether or not to dump the intermediate debug files generated when running the feature finder
    registerFlag_("dump", "dump intermediate files");

    // Options for testing only. Remove when no longer needed
    registerFlag_("pick_only", "[DEBUG] stop the tool after peak picking");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    return {};
  }

  // The main entry point.
  TOPPBase::ExitCodes main_(int argc, const char** argv) override
  {
    // printf("[FeatureFinderIonMobility::main_] Start\n");

    String inputPrefix = getStringOption_("in_prefix");
    Int numBins = getIntOption_("num_bins");
    bool useOffset = getFlag_("use_offset");
    String outputFile = getStringOption_("out");
    debugMode_ = getFlag_("dump");

    if (outputFile.empty())
    {
      outputFile = inputPrefix + ".featureXML";
    }

    // TODO: adjust the parameters for peak picking
    PeakPickerHiRes pphr;
    pphr.setLogType(log_type_);
    Param pphrParam;
    // Do not touch this parameter; it is required to force peak picking
    pphrParam.setValue("ms_levels", ListUtils::create<Int>("1"));
    // The following parameters can be played with to see what gets better results
    // pphrParam.setValue("spacing_difference_gap", 0.0);
    // pphrParam.setValue("spacing_difference", 4.0);
    // pphrParam.setValue("missing", 4);
    pphr.setParameters(pphrParam);

    FeatureFinder ffim;
    ffim.setLogType(log_type_);

    // Holds the resulting features for each binning pass. They will be collapsed into a single FeatureMap later
    std::vector<FeatureMap> features;
    std::vector<FeatureMap> offsetFeatures;

    auto constructFilename = [&](Size pass, Size bin) -> String
    {
      String filename = inputPrefix + "-";
      filename += std::to_string(pass);
      filename += "-";
      filename += std::to_string(bin);
      // filename += ".mzML";
      return filename;
    };

    // Go through each bin and run the existing feature finding workflow
    // TODO: these loops absolutely need to be parallelized
    for (Size i = 1; i <= numBins; i++)
    {
      String filename = constructFilename(1, i) + ".mzML";
      FeatureMap binFeatures;
      processBin_(filename, pphr, ffim, binFeatures);
      /*try
      {
        FeatureXMLFile().load(constructFilename(1, i) + ".featureXML", binFeatures);
      }
      catch (Exception::FileNotFound e)
      {
        continue;
      }*/
      features.push_back(std::move(binFeatures));
    }

    if (useOffset)
    {
      for (Size i = 1; i <= numBins + 1; i++)
      {
        String filename = constructFilename(2, i) + ".mzML";
        FeatureMap binFeatures;
        processBin_(filename, pphr, ffim, binFeatures);
        /*try
        {
          FeatureXMLFile().load(constructFilename(2, i) + ".featureXML", binFeatures);
        }
        catch (Exception::FileNotFound e)
        {
          continue;
        }*/
        offsetFeatures.push_back(std::move(binFeatures));
      }
    }

    // If only one bin was used, then there is only one featureXML file; nothing more to do
    if (numBins == 1)
    {
      FeatureXMLFile().store(outputFile, features[0]);
      return EXECUTION_OK;
    }

    // If a second binning pass was required, try to map the features from one pass to the other, followed by mapping
    // the features within the resulting pass. See the comments for the respective helper methods for a more in-depth
    // description of the mapping algorithms.
    if (useOffset)
    {
      /*FeatureMap result = removeDuplicateFeaturesInterPass_(features, offsetFeatures);
      FeatureXMLFile().store(outputFile, result);*/

      std::vector<FeatureMap> mergedPass = std::move(findFeaturesInAdjacentPass_(features, offsetFeatures));
      FeatureMap result = findFeaturesInAdjacentBins_(mergedPass);
      FeatureXMLFile().store(outputFile, result);
    }
    // Just map the features within the single pass
    else
    {
      FeatureMap result = findFeaturesInAdjacentBins_(features);
      FeatureXMLFile().store(outputFile, result);
    }

    return EXECUTION_OK;
  }

private:
  // Runs an existing feature finding workflow on a given bin.
  // <binName>: the path to the bin (input mzML file) to process.
  // <pphr>: the peak picker to run on the bin.
  // <ffim>: the feature finder to run on the bin.
  // <features>: the found features for the given bin.
  void processBin_(const String& binName, PeakPickerHiRes& pphr, FeatureFinder& ffim, FeatureMap& features)
  {
    // printf("[FeatureFinderIonMobility::processBin_] Start for \"%s\"\n", binName.c_str());

    if (!File::exists(binName))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __FUNCTION__, binName);
    }

    // Only load MS1 scans to save time and memory (since MS2 scans are currently not supported)
    PeakFileOptions options;
    options.setMSLevels(std::vector<Int>(1, 1));
    
    MzMLFile file;
    file.getOptions() = options;
    file.setLogType(log_type_);

    MSExperiment exp;
    file.load(binName, exp);
    exp.updateRanges();

    // Peak picking step
    MSExperiment pickedExp;
    pphr.pickExperiment(exp, pickedExp, false);
    /*for (MSSpectrum spec : exp.getSpectra())
    {
      MSSpectrum newSpec;
      pphr.pick(spec, newSpec);
      pickedExp.addSpectrum(newSpec);
    }*/

    if (debugMode_)
    {
      addDataProcessing_(pickedExp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
      Size dot = binName.rfind(".");
      String pickedBinName = binName.substr(0, dot) + "-picked.mzML";
      file.store(pickedBinName, pickedExp);
    }

    if (getFlag_("pick_only"))
    {
      return;
    }

    // Feature finding step
    FeatureMap seeds;
    FeatureMap tempFeatures;
    tempFeatures.setPrimaryMSRunPath({binName}, pickedExp);

    // Parameters for FeatureFinderCentroided. Can be played with
    Param ffimParam = getParam_().copy("algorithm:", true);
    ffimParam.setValue("mass_trace:min_spectra", 7);
    ffimParam.setValue("mass_trace:max_missing", 1);
    ffimParam.setValue("seed:min_score", 0.65);
    ffimParam.setValue("feature:min_score", 0.6);

    pickedExp.updateRanges();
    // If no peaks exist after the peak picking step (maybe this is near the first or last bin, which may not have much
    // data), there is nothing more to do. Trying to keep going will result in a crash.
    if (!experimentHasPeaks_(pickedExp))
    {
      return;
    }

    // Actually run FeatureFinderCentroided here and make sure the resulting features have unique IDs
    ffim.run(FeatureFinderAlgorithmPicked::getProductName(), pickedExp, tempFeatures, ffimParam, seeds);
    tempFeatures.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // Do a bit of cleaning here, removing any potentially "duplicate" features. See the comments of this helper method
    // for a more in-depth explanation.
    removeDuplicateFeaturesBin_(tempFeatures, features);
    features.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    if (debugMode_)
    {
      addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));
      Size dot = binName.rfind(".");
      String binFeaturesName = binName.substr(0, dot) + ".featureXML";

      FeatureXMLFile fFile;
      // To reduce file size, expand each feature's convex hull to its bounding box
      if (debug_level_ < 5)
      {
        for (Feature& f : features)
        {
          f.getConvexHull().expandToBoundingBox();
          for (auto& hull : f.getConvexHulls())
          {
            hull.expandToBoundingBox();
          }
          f.getSubordinates().clear();
        }
      }

      fFile.store(binFeaturesName, features);
    }
  }

  // Checks to see if an MSExperiment object has any data.
  // <exp>: the MSExperiment object to check.
  // Returns true if any data points are present in the experiment; false otherwise.
  bool experimentHasPeaks_(MSExperiment& exp) const
  {
    for (Size i = 0; i < exp.getNrSpectra(); i++)
    {
      const MSSpectrum& spec = exp.getSpectrum(i);
      if (spec.size() > 0)
      {
        return true;
      }
    }
    return false;
  }

  // "Cleans" the given bin by removing "duplicate" (or "similar") features. For each set of similar features, the
  // single feature with the greatest convex hull by size is chosen to "represent" the set, remaining in the bin while
  // all the other features are removed.
  // <tempFeatures>: the bin to clean.
  // <features>: the cleaned bin of features (instead of removing redundant features, the representatives will be added
  //   to this feature map).
  void removeDuplicateFeaturesBin_(FeatureMap& tempFeatures, FeatureMap& features)
  {
    // printf("[FeatureFinderIonMobility::removeDuplicateFeaturesBin_] Start\n");

    // Ensure the features are sorted by ascending RT so that they can be binary searched.
    tempFeatures.sortByRT();

    for (Size i = 0; i < tempFeatures.size(); i++)
    {
      // For each feature in the set, look for similar features
      const auto& feature = tempFeatures[i];
      Size firstIdx = findLeftRT_(tempFeatures, feature.getRT() - RT_THRESHOLD_);

      // Track the feature with the greatest convex hull
      double maxArea = polygonArea_(feature.getConvexHull().getHullPoints());
      Feature maxFeature = feature;

      // Looking for similar features
      for (Size j = firstIdx; j < tempFeatures.size(); j++)
      {
        // Skip the same feature
        if (i == j)
        {
          continue;
        }

        const auto& secondFeature = tempFeatures[j];
        // Check to see if we are no longer looking at similar features
        if (secondFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
        {
          break;
        }

        if (similarFeatures_(feature, secondFeature))
        {
          // Update the representative feature if necessary
          double area = polygonArea_(secondFeature.getConvexHull().getHullPoints());
          if (area > maxArea)
          {
            maxArea = area;
            maxFeature = secondFeature;
          }
        }
      }
      
      // Add the representative feature to the new set of features if it hasn't already been added
      // TODO: use an unordered set to check for duplicates instead?
      if (std::find(features.begin(), features.end(), maxFeature) == features.end())
      {
        features.push_back(maxFeature);
      }
    }
  }

  // Planned for deprecation
  FeatureMap removeDuplicateFeaturesInterPass_(std::vector<FeatureMap>& features, std::vector<FeatureMap>& offsetFeatures)
  {
    printf("[FeatureFinderIonMobility::removeDuplicateFeaturesInterPass_] Start\n");

    std::vector<Size> firstIdxs;
    printf("[FeatureFinderIonMobility::removeDuplicateFeaturesInterPass_] Start intra 1\n");
    FeatureMap firstPass = removeDuplicateFeaturesIntraPass_(features, firstIdxs);

    if (debugMode_)
    {
      String filename = getStringOption_("in") + "-pass1.featureXML";
      FeatureXMLFile().store(filename, firstPass);
    }

    if (offsetFeatures.empty())
    {
      firstPass.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      return firstPass;
    }

    std::vector<Size> secondIdxs;
    printf("[FeatureFinderIonMobility::removeDuplicateFeaturesInterPass_] Start intra 2\n");
    FeatureMap secondPass = removeDuplicateFeaturesIntraPass_(offsetFeatures, secondIdxs);

    if (debugMode_)
    {
      String filename = getStringOption_("in") + "-pass2.featureXML";
      FeatureXMLFile().store(filename, secondPass);
    }

    // Match between the two passes
    FeatureMap result;

    auto matchAgainst = [&](FeatureMap& refFeatures, std::vector<Size>& refIdxs,
                            FeatureMap& againstFeatures, std::vector<Size>& againstIdxs,
                            bool isFirst)
    {
      for (Size i = 0; i < refFeatures.size(); i++)
      {
        Feature firstFeature = refFeatures[i];
        Size firstIdx = refIdxs[i];
        double maxArea = polygonArea_(firstFeature.getConvexHull().getHullPoints());
        Feature maxFeature = firstFeature;

        Size j = findLeftRT_(againstFeatures, firstFeature.getRT() - RT_THRESHOLD_);
        for (; j < againstFeatures.size(); j++)
        {
          Feature secondFeature = againstFeatures[j];
          Size secondIdx = againstIdxs[j];
          if (secondFeature.getRT() > firstFeature.getRT() + RT_THRESHOLD_)
          {
            break;
          }

          if (similarFeatures_(firstFeature, secondFeature))
          {
            if ((isFirst && (firstIdx == secondIdx || firstIdx + 1 == secondIdx)) ||
                (!isFirst && (firstIdx == secondIdx || firstIdx - 1 == secondIdx)))
            {
              double area = polygonArea_(secondFeature.getConvexHull().getHullPoints());
              if (area > maxArea)
              {
                maxArea = area;
                maxFeature = secondFeature;
              }
            }
          }
        }
        
        if (std::find(result.begin(), result.end(), maxFeature) == result.end())
        {
          result.push_back(maxFeature);
        }
      }
    };

    printf("[FeatureFinderIonMobility::removeDuplicateFeaturesInterPass_] Start inter 1\n");
    matchAgainst(firstPass, firstIdxs, secondPass, secondIdxs, true);
    printf("[FeatureFinderIonMobility::removeDuplicateFeaturesInterPass_] Start inter 2\n");
    matchAgainst(secondPass, secondIdxs, firstPass, firstIdxs, false);

    result.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    return result;
  }

  // Maps the features from one pass to the features in the second, if the second binning pass (offset pass) is being
  // used. In order for a feature in one pass to map to a feature in the second, the two features must be similar and
  // be found in "overlapping" bins. For example, for a feature in bin i in the first pass to be mapped to a feature in
  // the second pass, the feature must be in either bin i or bin (i + 1).
  // <features>: the features for the first pass.
  // <offsetFeatures>: the features for the second pass.
  // Returns a vector of <num_bins> feature maps, where each feature map holds the features mapped across passes for
  //   that particular bin.
  std::vector<FeatureMap> findFeaturesInAdjacentPass_(std::vector<FeatureMap>& features, std::vector<FeatureMap>& offsetFeatures)
  {
    // printf("[FeatureFinderIonMobility::findFeaturesInAdjacentPass_] Start\n");

    std::vector<FeatureMap> results(features.size());

    // Make sure the features in each bin are sorted by ascending RT so the bins can be binary searched
    for (Size i = 0; i < features.size(); i++)
    {
      features[i].sortByRT();
    }
    for (Size i = 0; i < offsetFeatures.size(); i++)
    {
      offsetFeatures[i].sortByRT();
    }

    // For each bin in the first pass
    for (Size bIdx = 0; bIdx < features.size(); bIdx++)
    {
      // For each feature in the current bin
      for (Size fIdx = 0; fIdx < features[bIdx].size(); fIdx++)
      {
        Feature feature = features[bIdx][fIdx];
        double maxArea = polygonArea_(feature.getConvexHull().getHullPoints());
        Feature maxFeature = feature;

        // Finds similar features to the current feature in the given feature map. It will continuously update the
        // feature with the largest convex hull by size, which will be representative of the set.
        auto checkOtherBin = [&](FeatureMap otherBin) -> void
        {
          // Find the first similar feature in the other bin and iterate through ascending RT
          Size i = findLeftRT_(otherBin, feature.getRT() - RT_THRESHOLD_);
          for (; i < otherBin.size(); i++)
          {
            Feature otherFeature = otherBin[i];
            // Stop when the other features are no longer similar
            if (otherFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
            {
              break;
            }

            // Update the feature with the largest convex hull if necessary
            if (similarFeatures_(feature, otherFeature))
            {
              double area = polygonArea_(otherFeature.getConvexHull().getHullPoints());
              if (area > maxArea)
              {
                maxArea = area;
                maxFeature = otherFeature;
              }
            }
          }
        };

        // Look for similar features in the "overlapping" bins
        checkOtherBin(offsetFeatures[bIdx]);
        checkOtherBin(offsetFeatures[bIdx + 1]);

        // Add the feature to the mapped bin if it has not already been added
        if (std::find(results[bIdx].begin(), results[bIdx].end(), maxFeature) == results[bIdx].end())
        {
          results[bIdx].push_back(maxFeature);
        }
      }
    }

    // Do the same for the bins and features in the second pass. The idea is that any features found here were not
    // found in the first pass, possibly due to those features being cut apart by the bin boundaries, which are not
    // the same here.
    for (Size bIdx = 0; bIdx < offsetFeatures.size(); bIdx++)
    {
      for (Size fIdx = 0; fIdx < offsetFeatures[bIdx].size(); fIdx++)
      {
        Feature feature = offsetFeatures[bIdx][fIdx];
        bool found = false;

        auto checkOtherBin = [&](FeatureMap otherBin) -> void
        {
          Size i = findLeftRT_(otherBin, feature.getRT() - RT_THRESHOLD_);
          for (; i < otherBin.size(); i++)
          {
            Feature otherFeature = otherBin[i];
            if (otherFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
            {
              break;
            }
            if (similarFeatures_(feature, otherFeature))
            {
              found = true;
              break;
            }
          }
        };

        if (bIdx - 1 >= 0)
        {
          checkOtherBin(features[bIdx - 1]);
        }
        if (bIdx < features.size())
        {
          checkOtherBin(features[bIdx]);
        }

        if (!found)
        {
          if (bIdx - 1 >= 0)
          {
            results[bIdx - 1].push_back(feature);
          }
          else
          {
            results[bIdx].push_back(feature);
          }
        }
      }
    }

    return results;
  }

  // Planned for deprecation
  FeatureMap removeDuplicateFeaturesIntraPass_(std::vector<FeatureMap>& features, std::vector<Size>& binIdxs)
  {
    FeatureMap result;

    for (Size i = 0; i < features.size(); i++)
    {
      features[i].sortByRT();
    }

    for (Size bIdx = 0; bIdx < features.size(); bIdx++)
    {
      for (Size fIdx = 0; fIdx < features[bIdx].size(); fIdx++)
      {
        Feature feature = features[bIdx][fIdx];
        std::vector<std::pair<Feature, Size>> similarFeatures;

        for (Size nextBIdx = bIdx + 1; nextBIdx < features.size(); nextBIdx++)
        {
          double maxArea = -1.0;
          Feature maxFeature;

          Size nextFIdx = findLeftRT_(features[nextBIdx], feature.getRT() - RT_THRESHOLD_);
          for (; nextFIdx < features[nextBIdx].size(); nextFIdx++)
          {
            Feature nextFeature = features[nextBIdx][nextFIdx];
            if (nextFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
            {
              break;
            }

            if (similarFeatures_(feature, nextFeature))
            {
              double area = polygonArea_(nextFeature.getConvexHull().getHullPoints());
              if (area > maxArea)
              {
                maxArea = area;
                maxFeature = nextFeature;
              }
            }
          }

          if (maxArea < 0)
          {
            break;
          }

          similarFeatures.emplace_back(maxFeature, nextBIdx);
        }

        for (Size nextBIdx = bIdx - 1; nextBIdx >= 0; nextBIdx--)
        {
          double maxArea = -1.0;
          Feature maxFeature;

          Size nextFIdx = findLeftRT_(features[nextBIdx], feature.getRT() - RT_THRESHOLD_);
          for (; nextFIdx < features[nextBIdx].size(); nextFIdx++)
          {
            Feature nextFeature = features[nextBIdx][nextFIdx];
            if (nextFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
            {
              break;
            }

            if (similarFeatures_(feature, nextFeature))
            {
              double area = polygonArea_(nextFeature.getConvexHull().getHullPoints());
              if (area > maxArea)
              {
                maxArea = area;
                maxFeature = nextFeature;
              }
            }
          }

          if (maxArea < 0)
          {
            break;
          }

          similarFeatures.emplace_back(maxFeature, nextBIdx);
        }

        double maxIntensity = feature.getIntensity();
        Size maxIdx = bIdx;
        Feature maxFeature = feature;

        for (auto& similarFeature : similarFeatures)
        {
          double intensity = similarFeature.first.getIntensity();
          if (intensity > maxIntensity)
          {
            maxIntensity = intensity;
            maxIdx = similarFeature.second;
            maxFeature = similarFeature.first;
          }
        }

        if (std::find(result.begin(), result.end(), maxFeature) == result.end())
        {
          result.push_back(maxFeature);
          binIdxs.push_back(maxIdx);
        }
      }
    }

    return result;
  }

  // Maps the features in a single pass together across bins. In order for a feature in one bin to map to a feature in
  // another, the two features must be similar and be found in adjacent bins. For example, a feature in bin i can only
  // be mapped to similar features found in bins (i - 1) or (i + 1). This is applied iteratively, meaning that a
  // "chain" of mapped features can be created; e.g., the same feature may be able to be mapped to a similar feature in
  // bins (i + 2), (i + 3), etc.
  // <features>: the features for the pass.
  // Returns a single, final feature map containing the mapped features.
  FeatureMap findFeaturesInAdjacentBins_(std::vector<FeatureMap>& features)
  {
    // printf("[FeatureFinderIonMobility::findFeaturesInAdjacentBins] Start\n");

    FeatureMap result;

    // Make sure the features in each bin are sorted by ascending RT so that they can bin binary searched.
    for (Size i = 0; i < features.size(); i++)
    {
      features[i].sortByRT();
    }

    // For each bin in the pass
    for (Size bIdx = 0; bIdx < features.size(); bIdx++)
    {
      // For each feature in the bin
      for (Size fIdx = 0; fIdx < features[bIdx].size(); fIdx++)
      {
        Feature feature = features[bIdx][fIdx];
        std::vector<Feature> similarFeatures;

        // Check the next bin (try to extend the chain forward)
        for (Size nextBIdx = bIdx + 1; nextBIdx < features.size(); nextBIdx++)
        {
          double maxArea = -1.0;
          Feature maxFeature;

          // Find the first similar feature in this bin
          Size nextFIdx = findLeftRT_(features[nextBIdx], feature.getRT() - RT_THRESHOLD_);
          for (; nextFIdx < features[nextBIdx].size(); nextFIdx++)
          {
            Feature nextFeature = features[nextBIdx][nextFIdx];
            // If there are no more similar features, stop looking in this bin
            if (nextFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
            {
              break;
            }

            if (similarFeatures_(feature, nextFeature))
            {
              double area = polygonArea_(nextFeature.getConvexHull().getHullPoints());
              if (area > maxArea)
              {
                maxArea = area;
                maxFeature = nextFeature;
              }
            }
          }

          // If no similar features were found in this bin, do not look in next bins
          if (maxArea < 0)
          {
            break;
          }

          similarFeatures.push_back(maxFeature);
        }

        // Do the same for previous bins (try to extend the chain backward)
        for (long long nextBIdx = bIdx - 1; nextBIdx >= 0; nextBIdx--)
        {
          double maxArea = -1.0;
          Feature maxFeature;

          Size nextFIdx = findLeftRT_(features[nextBIdx], feature.getRT() - RT_THRESHOLD_);
          for (; nextFIdx < features[nextBIdx].size(); nextFIdx++)
          {
            Feature nextFeature = features[nextBIdx][nextFIdx];
            if (nextFeature.getRT() > feature.getRT() + RT_THRESHOLD_)
            {
              break;
            }

            if (similarFeatures_(feature, nextFeature))
            {
              double area = polygonArea_(nextFeature.getConvexHull().getHullPoints());
              if (area > maxArea)
              {
                maxArea = area;
                maxFeature = nextFeature;
              }
            }
          }

          if (maxArea < 0)
          {
            break;
          }

          similarFeatures.push_back(maxFeature);
        }

        // Out of all the similar features found, find the representative feature and add it to the final list of
        // features
        double maxArea = polygonArea_(feature.getConvexHull().getHullPoints());
        Feature maxFeature = feature;

        for (auto& similarFeature : similarFeatures)
        {
          double area = polygonArea_(similarFeature.getConvexHull().getHullPoints());
          if (area > maxArea)
          {
            maxArea = area;
            maxFeature = similarFeature;
          }
        }

        if (std::find(result.begin(), result.end(), maxFeature) == result.end())
        {
          result.push_back(maxFeature);
        }
      }
    }

    return result;
  }

  // Runs a binary search on a given list of features to find the index of the feature with the target RT. If multiple
  // features with the target RT are present, the index of the first one (the lowest index) will be returned. If no
  // feature with the target RT is present, the number of features with RT less than the target RT will be returned.
  // <features>: the list of features to search. It must be sorted by ascending RT.
  // <target>: the target RT to search for.
  // Returns the index of the first feature with the target RT (if present), or the number of features with RT less
  //   than the target (if not present).
  Size findLeftRT_(const FeatureMap& features, double target) const
  {
    Int lo = 0;
    Int hi = features.size();

    while (lo < hi)
    {
      Int mid = Int((lo + hi) / 2);
      double result = features[mid].getRT();
      if (result < target)
      {
        lo = mid + 1;
      }
      else
      {
        hi = mid;
      }

      // if (features[lo].getRT() > target)
      // {
      //   return (lo > 0) ? lo - 1 : 0;
      // }
      // lo++;
    }

    return lo;
  }

  // Calculates the area of a convex polygon, using the shoelace formula. Intended to be used to get the size of a
  // feature's convex hull.
  // <points>: The coordinates of the vertices of the polygon.
  // Returns the area of the given polygon.
  double polygonArea_(const std::vector<DPosition<2>>& points) const
  {
    double area = 0.0;
    for (Size i = 0; i < points.size(); i++)
    {
      area += points[i].getX() * points[(i + 1) % points.size()].getY();
      area -= points[i].getY() * points[(i + 1) % points.size()].getX();
    }
    return std::abs(area) / 2.0;
  }

  // Checks if two features are "similar"; i.e., their RT and m/z are within some fixed, user-defined thresholds of
  // each other. If so, we may consider them to be "identical" or essentially the same feature.
  // <firstFeature>: the first feature to check.
  // <secondFeature>: the second feature to check.
  // Returns true if the two features are similar; false otherwise.
  bool similarFeatures_(const Feature& firstFeature, const Feature& secondFeature) const
  {
    bool inRT = (std::abs(firstFeature.getRT() - secondFeature.getRT()) < RT_THRESHOLD_);
    bool inMZ = (std::abs(firstFeature.getMZ() - secondFeature.getMZ()) < MZ_THRESHOLD_);
    return (inRT && inMZ);
  }

  bool debugMode_; // Determines if intermediate files should be dumped

  // The thresholds to use for determining if features are similar
  const double RT_THRESHOLD_ = 5.0;
  const double MZ_THRESHOLD_ = 0.01;
};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderIonMobility tool;
  return tool.main(argc, argv);
}
