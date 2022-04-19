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

// This header file contains all of the necessary functionality for the IonMobilityBinner TOPP tool. This code should
// either be moved into the OpenMS library if it needs to be used elsewhere, or moved back into the corresponding
// source file (IonMobilityBinner.cpp). No Doxygen documentation has been added because I was unable to build
// documentation on my local Windows machine.

// This TOPP tool takes a single mzML file as input, bins the data points in the ion mobility (IM) dimension, and
// outputs <numBins> mzML files, each corresponding to an individual bin.

#ifndef ION_MOBILITY_BINNER_H
#define ION_MOBILITY_BINNER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{
class IonMobilityBinner
{
public:
  IonMobilityBinner::IonMobilityBinner() :
      binSize_(0.0),
      imStart_(0.0),
      imEnd_(0.0),
      imDelta_(0.0),
      imOffset_(0.0)
  {
    bins_.clear();
    offsetBins_.clear();
  }

  // The main entry point. This method is invoked by the IonMobilityBinner TOPP tool.
  // <inputFile>: the path to the input mzML file (to be binned).
  // <outputPrefix>: the string to use as a prefix for the <numBins> output mzML files. The output format is
  //   "<outputPrefix>-<pass number>-<bin number>.mzML".
  // <numBins>: the number of bins to use.
  // <useOffset>: specifies whether or not to perform a second binning pass with the bin boundaries shifted by 50%. If
  //   enabled, the second pass will have (<numBins> + 1) bins.
  void IonMobilityBinner::run(const String& inputFile, const String& outputPrefix, Size numBins, bool useOffset)
  {
    // printf("[IonMobilityBinner::run] Start\n");

    if (!File::exists(inputFile))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __FUNCTION__, inputFile);
    }
    if (!exp_.openFile(inputFile))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, __FUNCTION__, inputFile);
    }

    bins_.resize(numBins);
    if (useOffset)
    {
      offsetBins_.resize(numBins + 1);
    }

    setupBins_();
    // If files with the same required names already exist (e.g., from old runs), those files will be removed. This is
    // because this tool continuously updates the output mzML files by adding new data to them, so this prevents adding
    // to old data.
    removeOldExperiments_(outputPrefix);

    // Process the mzML file by binning one spectrum at a time
    // TODO: look into parallelizing this. It may be difficult since each bin needs its own lock
    for (Size i = 0; i < exp_.getNrSpectra(); i++)
    {
      MSSpectrum spec = exp_.getSpectrum(i);
      // The tool currently only supports MS1 spectra
      // TODO: add support for MS2 spectra
      if (spec.getMSLevel() != 1)
      {
        continue;
      }
      binSpectrum_(spec);

      // Every "this many" spectra are binned, dump the in-memory MSExperiment objects to the corresponding output mzML
      // files and then clear the in-memory objects. This prevents the tool from using too much memory at once.
      // TODO: benchmark this to find an appropriate value to use
      if ((i + 1) % 10000 == 0)
      {
        writeExperiments_(outputPrefix);
      }
    }

    // Need to dump any remaining bins to the output mzML files
    writeExperiments_(outputPrefix);
  }

protected:
  // Prepares for binning by getting the appropriate boundaries for the bins, ensuring that all are of equal size.
  void IonMobilityBinner::setupBins_()
  {
    // printf("[IonMobilityBinner::setupBins_] Start\n");

    // Get the smallest and largest IM values in the entire experiment
    imStart_ = std::numeric_limits<double>::max();
    imEnd_ = std::numeric_limits<double>::min();

    for (Size i = 0; i < exp_.getNrSpectra(); i++)
    {
      MSSpectrum spec = exp_.getSpectrum(i);
      const auto& fda = spec.getFloatDataArrays();

      // All spectra should have IM data. Maybe error out here if no IM data is present?
      if (!spec.containsIMData())
      {
        continue;
      }
      // This is the appropriate index in the float data arrays for IM data
      Size imDataIdx = spec.getIMData().first;

      for (float im : fda[imDataIdx])
      {
        if (im < imStart_)
        {
          imStart_ = im;
        }
        if (im > imEnd_)
        {
          imEnd_ = im;
        }
      }
    }

    // See the comments for the private attributes of this class for an explanation as to what these are for
    imDelta_ = imEnd_ - imStart_;
    binSize_ = imDelta_ / bins_.size();
    imOffset_ = imStart_ + (binSize_ / 2.0);
  }

  // Removes old output mzML files to ensure that running the tool again does not add to useless data.
  // <outputPrefix>: the output prefix of the files to check for and remove.
  void IonMobilityBinner::removeOldExperiments_(const String& outputPrefix) const
  {
    // printf("[IonMobilityBinner::removeOldExperiments_] Start\n");

    // Generates the correct filename format for a given pass and bin number and removes that file.
    auto removeExps = [&](Size passNo, Size numBins) -> void {
      for (Size i = 1; i <= numBins; i++)
      {
        String filename = outputPrefix + "-";
        filename += std::to_string(passNo);
        filename += "-";
        filename += std::to_string(i);
        filename += ".mzML";

        // This is safe, even if the file does not exist
        File::remove(filename);
      }
    };

    removeExps(1, bins_.size());
    removeExps(2, offsetBins_.size());
  }

  // Bins a given spectrum. New spectra are created (one for each bin), and data points in the original spectrum are
  // assigned to one of the new spectra. Then, each new spectra is checked for duplicate points (combining them if
  // necessary), and at the end, the new spectra are added to the private in-memory MSExperiment objects (which
  // represent the output mzML files).
  // <spec>: the spectrum to bin.
  void IonMobilityBinner::binSpectrum_(MSSpectrum& spec)
  {
    // printf("[IonMobilityBinner::binSpectrum_] Start for RT %.6f\n", spec.getRT());

    // Ensure that the data points in the spectrum are sorted by ascending m/z
    spec.sortByPosition();

    // All spectra should have IM data. Maybe error out here if no IM data is present?
    if (!spec.containsIMData())
    {
      return;
    }
    Size imDataIdx = spec.getIMData().first;
    // imData[i] is the IM value for the data point in <spec> at position i; spec[i]
    const auto& imData = spec.getFloatDataArrays()[imDataIdx];

    // Create new temporary "bins" to split the data points to
    std::vector<MSSpectrum> tempBins(bins_.size());
    std::vector<MSSpectrum> tempOffsetBins(offsetBins_.size());
    // Also store the data points' IM values in parallel. This is implemented as a lookup table, where each interior
    // vector v represents a bin, v[i] for the i-th data point in v represents the index j in <spec> for the
    // corresponding data point, and imData[j] represents its corresponding IM value.
    std::vector<std::vector<Size>> imIndexes(bins_.size());
    std::vector<std::vector<Size>> imOffsetIndexes(offsetBins_.size());

    // Iterate through each data point (by ascending m/z)
    for (Size i = 0; i < imData.size(); i++)
    {
      // Calculate the appropriate bin index and add the data point to the corresponding temporary bin
      Size binIdx = Size((imData[i] - imStart_) / binSize_);
      if (binIdx >= bins_.size())
      {
        binIdx = bins_.size() - 1;
      }
      tempBins[binIdx].push_back(spec[i]);
      imIndexes[binIdx].push_back(i);

      // Do the same for the second binning pass, if required
      if (!offsetBins_.empty())
      {
        binIdx = Size((imData[i] - imOffset_) / binSize_) + 1;
        if (imData[i] < imOffset_)
        {
          binIdx = 0;
        }
        else if (binIdx > bins_.size())
        {
          binIdx = bins_.size();
        }
        tempOffsetBins[binIdx].push_back(spec[i]);
        imOffsetIndexes[binIdx].push_back(i);
      }
    }

    // Each temporary bin (represented as a spectrum) may have multiple data points at the same m/z (which will cause
    // problems), which need to be combined into a single data point. For such data points, the resulting single data
    // point will have its intensity equal to the sum of intensities for all affected data points, and its IM value
    // will be an intensity-weighted average of those points.
    auto combinePeaks = [&](const std::vector<MSSpectrum>& bins,
                            const std::vector<std::vector<Size>>& indexes,
                            std::vector<MSExperiment>& exps) -> void {
      // We will consider points with an m/z difference of this epsilon to actually have the same m/z. This probably
      // doesn't make a difference since the points are already binned in m/z.
      const double MZ_EPSILON = 0.001;

      for (Size i = 0; i < bins.size(); i++)
      {
        const MSSpectrum& tempSpec = bins[i];
        if (tempSpec.empty())
        {
          continue;
        }

        MSSpectrum newSpec;
        newSpec.setRT(spec.getRT());
        DataArrays::FloatDataArray newImData;

        Size mzStart = 0;                    // The index of the first peak with the current m/z
        double mzCurr = tempSpec[0].getMZ(); // The current m/z
        double runningIntensity = 0.0;       // The sum of intensities for all peaks with the current m/z
        float weightedIM = 0.0;              // The intensity-weighted IM value for all peaks with the current m/z

        // Go through each data point in the current temporary bin
        for (Size j = 0; j < tempSpec.size(); j++)
        {
          const Peak1D& peak = tempSpec[j];
          // Found a data point with the same m/z, so we need to combine it with the current one
          if (mzCurr - MZ_EPSILON <= peak.getMZ() && peak.getMZ() <= mzCurr + MZ_EPSILON)
          {
            runningIntensity += peak.getIntensity();
            weightedIM += peak.getIntensity() * imData[indexes[i][j]];
          }
          // Reached a data point with a new m/z, so save the set of previous data points
          else
          {
            Peak1D newPeak(tempSpec[mzStart]);
            newPeak.setIntensity(runningIntensity);
            newSpec.push_back(newPeak);

            weightedIM /= runningIntensity;
            newImData.push_back(weightedIM);

            // Reset the counters for the new set
            mzStart = j;
            mzCurr = peak.getMZ();
            runningIntensity = peak.getIntensity();
            weightedIM = peak.getIntensity() * imData[indexes[i][j]];
          }
        }

        // Take care of the last m/z set
        Peak1D newPeak(tempSpec[mzStart]);
        newPeak.setIntensity(runningIntensity);
        newSpec.push_back(newPeak);

        weightedIM /= runningIntensity;
        newImData.push_back(weightedIM);

        // Add the ion mobility values as a new float data array
        newImData.setName("Ion Mobility");
        std::vector<DataArrays::FloatDataArray> fda = {newImData};
        newSpec.setFloatDataArrays(fda);
        // Add the spectrum to the in-memory MSExperiment (to become part of the appropriate output mzML file)
        exps[i].addSpectrum(newSpec);
      }
    };

    combinePeaks(tempBins, imIndexes, bins_);
    combinePeaks(tempOffsetBins, imOffsetIndexes, offsetBins_);
  }

  // Writes the current state of the in-memory MSExperiment objects to their corresponding output mzML files.
  // <outputPrefix>: the output prefix to use when writing the files.
  void IonMobilityBinner::writeExperiments_(const String& outputPrefix)
  {
    // printf("[IonMobilityBinner::writeExperiments_] Start\n");

    auto writeExps = [&](std::vector<MSExperiment>& bins, Size passNo) -> void {
      for (Size i = 0; i < bins.size(); i++)
      {
        String filename = outputPrefix + "-";
        filename += std::to_string(passNo);
        filename += "-";
        filename += std::to_string(i + 1);
        filename += ".mzML";

        MSExperiment exp;
        MzMLFile file;

        // Load the existing mzML files because we need to add to them
        try
        {
          file.load(filename, exp);
        }
        // File not found, so it will be created
        catch (Exception::FileNotFound e)
        {
        }

        // Add the spectra from the current experiment to the loaded experiment
        for (const MSSpectrum& spec : bins[i])
        {
          exp.addSpectrum(spec);
        }
        // Store the new resulting experiment and clear the current/in-memory one
        file.store(filename, exp);
        bins[i].clear(true);
      }
    };

    writeExps(bins_, 1);
    writeExps(offsetBins_, 2);
  }

private:
  OnDiscMSExperiment exp_;               // The mzML file to bin. It must be indexed
  std::vector<MSExperiment> bins_;       // Holds binning results
  std::vector<MSExperiment> offsetBins_; // Holds binning results for the second pass (optional; 50% offset)

  double binSize_;  // Size of a normal bin (the first and last offset bins are half this size)
  double imStart_;  // Smallest IM value in the input file
  double imEnd_;    // Largest IM value in the input file
  double imDelta_;  // Difference between the largest and smallest IM values
  double imOffset_; // IM value defining where the second bin starts when using the second pass (optional; 50% offset)

}; // class IonMobilityBinner

} // namespace OpenMS

#endif // ION_MOBILITY_BINNER_H
