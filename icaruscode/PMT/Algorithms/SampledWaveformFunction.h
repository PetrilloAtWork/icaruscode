/**
 * @file   icaruscode/PMT/Algorithms/SampledWaveformFunction.h
 * @brief  Pulse from one photoelectron as a train of samples.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 6, 2022
 *
 * This library is header only.
 *
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_SAMPLEDWAVEFORMFUNCTION_H
#define ICARUSCODE_PMT_ALGORITHMS_SAMPLEDWAVEFORMFUNCTION_H

// library header
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds

// C++ standard library
#include <ostream> // std::ostream
#include <vector>
#include <string>
#include <algorithm> // std::transform()
#include <numeric> // std::reduce()
#include <cmath> // std::round(), std::floor(), ...
#include <cstddef> // std::ptrdiff_t, std::size_t
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  using namespace util::quantities::time_literals; // ""_ns
  template <typename T> class SampledWaveformFunction;
}

// -----------------------------------------------------------------------------
/**
 * @brief Describes the waveform from a single photoelectron.
 * @tparam T type of time unit to be used
 *
 * This functor (class behaving like a function) describes the shape of the
 * response to a single photoelectron in a non-analytical form from a sequence
 * of samples.
 * 
 * The peak time is assigned to the sample with the largest value.
 * 
 * See more details in the constructor.
 * 
 * @note Currently the shape is hard-coded; it is possible to extend this object
 *       to receive the sample data from a text file with proper format
 *       (which needs to include all the information for `WaveformSpecs_t`).
 */
template <typename T>
class icarus::opdet::SampledWaveformFunction
  : public icarus::opdet::PhotoelectronPulseFunction<T>
{
  using Base_t = icarus::opdet::PhotoelectronPulseFunction<T>;

    public:
  /// Type for ADC counts (floating point).
  using ADCcount = typename Base_t::ADCcount;

  using Time = typename Base_t::Time; ///< Type of time being used.

  /**
   * @brief Constructor: initializes from an hard-coded shape.
   * @param peakTime time to assign to the peak
   * @param gain the gain of the optical detector
   * 
   * The hard-coded shape is shifted in time so that evaluation at `peakTime`
   * (`evaluateAt(peakTime)`) returns the peak amplitude; more precisely,
   * `peakTime` is set to match the start of the largest sample.
   * 
   * The polarity of the waveform is deduced by the value at peak.
   * 
   * All the rest is beautifully hard-coded (see `Waveform`) so far.
   */
  SampledWaveformFunction(Time peakTime, float gain);

  /// @{
  /// @name Parameter accessors.

  /// Returns the gain the waveform is representing.
  float gain() const { return fGain; }

  /// @}

    private:
  
  /// Internal structure to specify the waveform shape and features.
  struct WaveformSpecs_t {
    std::string name;           ///< Name or description of this waveform.
    std::vector<float> samples; ///< Samples [mV]
    Time sampleDuration;        ///< Time extension of one sample.
    float gain;                 ///< Gain of this waveform.
  };
  
  WaveformSpecs_t const* fSource = nullptr; ///< Which waveform to import.
  
  std::vector<ADCcount> const fSamples; ///< All samples.
  
  float const fGain; ///< The gain this waveform represents.
  
  std::size_t const fPeakSample; ///< The sample with the absolute peak.
  
  Time const fRefTime; ///< The time of the start of sample #0.
  
  /// The duration of each sample.
  Time sampleDuration() const { return fSource->sampleDuration; }
  
  
  // --- BEGIN -- Interface implementation -------------------------------------
  /**
   * @brief Evaluates the pulse at the given time.
   * @param time time to evaluate the shape at
   *
   * The scale of the time is defined by the peak time passed at construction.
   */
  virtual ADCcount doEvaluateAt(Time time) const override;

  /// Returns the time at which the first peak is found.
  virtual Time doPeakTime() const override
    { return fRefTime + fPeakSample * sampleDuration(); }

  /// Returns the amplitude of the first peak in ADC.
  virtual ADCcount doPeakAmplitude() const override
    { return fSamples[fPeakSample]; }

  /**
   * @brief Prints on stream the parameters of this shape.
   * @param out the stream to write into
   * @param indent indentation string, prepended to all lines except first
   * @param indentFirst indentation string prepended to the first line
   */
  virtual void doDump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const override;

  // --- END -- Interface implementation -------------------------------------
  
  
  /// Returns whether a sample with the specified index is within the range.
  bool hasSample(std::ptrdiff_t index) const
    { return (index >= 0) && (std::size_t(index) < fSamples.size()); }
  
  
  /// Returns the integral of the waveform.
  ADCcount integral() const;
  
  /**
   * @brief Transforms the input waveform.
   * @param waveform input waveform (in millivolt and for a known gain)
   * @param targetGain the desired gain
   * @return a sequence of samples in ADC
   * 
   * The returned waveform has the same time domain as the input one, but is
   * expressed in ADC instead of voltage (conversion is perfectly linear,
   * 2 V across 14 bits), and rescaled to meet the target gain.
   */
  std::vector<ADCcount> buildSamples
    (WaveformSpecs_t const& waveform, float targetGain) const;
  
  /// Returns the index of the sample under the peak of the waveform.
  static std::size_t findPeak(std::vector<ADCcount> const& samples);
  
  /**
   * @brief The waveform currently used for the shape.
   * 
   * This waveform was extracted by Andrea Scarpelli in Winter 2022 from channel
   * 258 of run TODO.
   * The waveform is in millivolt, normalized for gain @f$ 7\cdot 10^{6} @f$,
   * and sampled at 1 ns.
   */
  static const WaveformSpecs_t Waveform;
  
}; // class icarus::opdet::SampledWaveformFunction<>


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T>
typename icarus::opdet::SampledWaveformFunction<T>::WaveformSpecs_t const
icarus::opdet::SampledWaveformFunction<T>::Waveform {
    /*
     * Post-processing by G. Petrillo from the original by A. Scarpelli:
     *  * removed introductory pedestal (68 samples, most "negative")
     *  * removed trailing samples (8 samples)
     *  * changed the sign of the samples (at that point, all were positive)
     */
    "average waveform of channel 258 of run XXXX",
    {                 // samples
      -0.002822, -0.003891, -0.004395, -0.005078, -0.005970, -0.009117, -0.019188, -0.037225,
      -0.100708, -0.204642, -0.468557, -0.949630, -1.937902, -2.274567, -2.572617, -2.614342,
      -2.893903, -2.271094, -2.054658, -1.835824, -1.870135, -1.386591, -1.199377, -1.023026,
      -1.009500, -0.737781, -0.644197, -0.568364, -0.588811, -0.444252, -0.395502, -0.360269,
      -0.382917, -0.306847, -0.294776, -0.284055, -0.309246, -0.246513, -0.234806, -0.224057,
      -0.251971, -0.210322, -0.199377, -0.185362, -0.204823, -0.163578, -0.151966, -0.146062,
      -0.159384, -0.122480, -0.115287, -0.101905, -0.108998, -0.088126, -0.087022, -0.082043,
      -0.088900, -0.068711, -0.065204, -0.067568, -0.085125, -0.070700, -0.067736, -0.061950,
      -0.072639, -0.064006, -0.069002, -0.069103, -0.077837, -0.066207, -0.069025, -0.071249,
      -0.084990, -0.075399, -0.076609, -0.077897, -0.097297, -0.089280, -0.089392, -0.095228,
      -0.121725, -0.106375, -0.114878, -0.120799, -0.142019, -0.120407, -0.115876, -0.105345,
      -0.105122, -0.077298, -0.073842, -0.073920, -0.092379, -0.076822, -0.075786, -0.076273,
      -0.092194, -0.083398, -0.087000, -0.084602, -0.094339, -0.079768, -0.082171, -0.082227,
      -0.093348, -0.077208, -0.072358, -0.069932, -0.083814, -0.073736, -0.069109, -0.064560,
      -0.078022, -0.062717, -0.064566, -0.063569, -0.069659, -0.059463, -0.058578, -0.056018,
      -0.063537, -0.054170, -0.053716, -0.056153, -0.070533, -0.055027, -0.053329, -0.052344,
      -0.059812, -0.053609, -0.055133, -0.053189, -0.062030, -0.052971, -0.053447, -0.054310,
      -0.063072, -0.052153, -0.049901, -0.051145, -0.064136, -0.055155, -0.051823, -0.050092,
      -0.061044, -0.049980, -0.052495, -0.056393, -0.064142, -0.051834, -0.057093, -0.054528,
      -0.063901, -0.059323, -0.053368, -0.050719, -0.063195, -0.051951, -0.046451, -0.049380,
      -0.057207, -0.049475, -0.050820, -0.051503, -0.057319, -0.052512, -0.054595, -0.056494,
      -0.063867, -0.050882, -0.049666, -0.051156, -0.064769, -0.057105, -0.052517, -0.048114,
      -0.057476, -0.045364, -0.046496, -0.049022, -0.053504, -0.044664, -0.045448, -0.045694,
      -0.053454, -0.044854, -0.044356, -0.043023, -0.055991, -0.044753, -0.042317, -0.041969,
      -0.047841, -0.042082, -0.042530, -0.044322, -0.048698, -0.039466, -0.038838, -0.040636,
      -0.048962, -0.039376, -0.039270, -0.036256, -0.046514, -0.037236, -0.033763, -0.033444,
      -0.041691, -0.034996, -0.035601, -0.038345, -0.041052, -0.032985, -0.034189, -0.036687,
      -0.040666, -0.032111, -0.030940, -0.030173, -0.041310, -0.035528, -0.032172, -0.033915,
      -0.039904, -0.032688, -0.037298, -0.035130, -0.039960, -0.032985, -0.039040, -0.039438,
      -0.044497, -0.034195, -0.033690, -0.033791, -0.047304, -0.038441, -0.036693, -0.036816,
      -0.043215, -0.038670, -0.040015, -0.043325, -0.046043, -0.035791, -0.041527, -0.039247,
      -0.045640, -0.039522, -0.035091, -0.034912, -0.046228, -0.039337, -0.036508, -0.037488,
      -0.043803, -0.036396, -0.039281, -0.036934, -0.042335, -0.037948, -0.038889, -0.035802,
      -0.043741, -0.034576, -0.033982, -0.032172, -0.040100, -0.033422, -0.032419, -0.029562,
      -0.036213, -0.028890, -0.030162, -0.031668, -0.033843, -0.025680, -0.028744, -0.028851,
      -0.033109, -0.027927, -0.026565, -0.025087, -0.035048, -0.027389, -0.025602, -0.024672,
      -0.027934, -0.025182, -0.027473, -0.026946, -0.028113, -0.025944, -0.027053, -0.025764,
      -0.031099, -0.026257, -0.023182, -0.025843, -0.031563, -0.026661, -0.023395, -0.023877,
      -0.027502, -0.022924, -0.028341, -0.024078, -0.026785, -0.021149, -0.023866, -0.023345,
      -0.027334, -0.022432, -0.021451, -0.021171, -0.029760, -0.027058, -0.024902, -0.024112,
      -0.029541, -0.025501, -0.027411, -0.029327, -0.028079, -0.023888, -0.026666, -0.026734,
      -0.031580, -0.024286, -0.020824, -0.021782, -0.033076, -0.025232, -0.021496, -0.023524,
      -0.027446, -0.023972, -0.025131, -0.025781, -0.030101, -0.025372, -0.026504, -0.027831,
      -0.031967, -0.022975, -0.021395, -0.023518, -0.030465, -0.023468, -0.019782, -0.019362,
      -0.025631, -0.018757, -0.020639, -0.023591, -0.019610, -0.016920, -0.019059, -0.018847,
      -0.024265, -0.017967, -0.015581, -0.016191, -0.025637, -0.019698, -0.018113, -0.018735,
      -0.021761, -0.020667, -0.021843, -0.021233, -0.025833, -0.019748, -0.020241, -0.020941,
      -0.024438, -0.016802, -0.018045, -0.017737, -0.023682, -0.021826, -0.019592, -0.016074,
      -0.021206, -0.018152, -0.019295, -0.022191, -0.021884, -0.017379, -0.021401, -0.020006,
      -0.022007, -0.018141, -0.017721, -0.018051, -0.026864, -0.019844, -0.019009, -0.020437,
      -0.025984, -0.022342, -0.023277, -0.022028, -0.022758, -0.020852, -0.022235, -0.022336,
      -0.026035, -0.019866, -0.020387, -0.021698, -0.027967, -0.023328, -0.019592, -0.021670,
      -0.022489, -0.019508, -0.020717, -0.022135, -0.022097, -0.018639, -0.020583, -0.019715,
      -0.021957, -0.018314, -0.018180, -0.016914, -0.023721, -0.019508, -0.017015, -0.017558,
      -0.019571, -0.019401, -0.023524, -0.022263, -0.025329, -0.022006, -0.024991, -0.028201,
      -0.034235, -0.032654, -0.032531, -0.034665, -0.051947, -0.048210, -0.049447, -0.049526,
      -0.061475, -0.055447, -0.056421, -0.055010, -0.060299, -0.047761, -0.050125, -0.045538,
      -0.047880, -0.034172, -0.032268, -0.027697, -0.036829, -0.029854, -0.024056, -0.022644,
      -0.025234, -0.018482, -0.020006, -0.016656, -0.014826, -0.010517, -0.010450, -0.009464,
      -0.007090, -0.006450, -0.002070, -0.000575, -0.007544, -0.006462, -0.003034, -0.003134,
      -0.007511, -0.006591, -0.009262, -0.011044, -0.008709, -0.007705, -0.010075, -0.013760,
      -0.015258, -0.014735, -0.011369, -0.008786, -0.019151, -0.015178, -0.013027, -0.014673,
      -0.015375, -0.013497, -0.014875, -0.017256, -0.015213, -0.013945, -0.014696, -0.011122,
      -0.016154, -0.012433, -0.012310, -0.011453, -0.017812, -0.012427, -0.012013, -0.010377,
      -0.014003, -0.013497, -0.013346, -0.014035, -0.013896, -0.009212, -0.011010, -0.013990,
      -0.012585, -0.009974, -0.007302, -0.007167, -0.014087, -0.010114, -0.006355, -0.007369,
      -0.011689, -0.006260, -0.012416, -0.010764, -0.011112, -0.009666, -0.010747, -0.009705,
      -0.011544, -0.009022, -0.005958, -0.008893, -0.014977, -0.012595, -0.012763, -0.010080,
      -0.012810, -0.010965, -0.014786, -0.013363, -0.013986, -0.011178, -0.013895, -0.015469,
      -0.014356, -0.013206, -0.011772, -0.012629, -0.018489, -0.016494, -0.012052, -0.012802,
      -0.014753, -0.013643, -0.014057, -0.017586, -0.014781, -0.011777, -0.016511, -0.015234,
      -0.016865, -0.013794, -0.010635, -0.010119, -0.019201, -0.014220, -0.011890, -0.012550,
      -0.015341, -0.011811, -0.014550, -0.013077, -0.012490, -0.011475, -0.013419, -0.013878,
      -0.011975, -0.010612, -0.011357, -0.008377, -0.015806, -0.015721, -0.010607, -0.010820,
      -0.013599, -0.013542, -0.016578, -0.015693, -0.014075, -0.011598, -0.014825, -0.013867,
      -0.017022, -0.013396, -0.012360, -0.011990, -0.018562, -0.015206, -0.012388, -0.012640,
      -0.015045, -0.016124, -0.016505, -0.012102, -0.013078, -0.011800, -0.012752, -0.013301,
      -0.016798, -0.011301, -0.012007, -0.011268, -0.015694, -0.015038, -0.012478, -0.011548,
      -0.016938, -0.013357, -0.015262, -0.012976, -0.013258, -0.010680, -0.013895, -0.011262,
      -0.012602, -0.012774, -0.009610, -0.011380, -0.017173, -0.014360, -0.011441, -0.013167,
      -0.014546, -0.013183, -0.014802, -0.014264, -0.013521, -0.011688, -0.014421, -0.014494,
      -0.019016, -0.012298, -0.012366, -0.012354, -0.019246, -0.016236, -0.013962, -0.011593,
      -0.017907, -0.011828, -0.013968, -0.016500, -0.016972, -0.011290, -0.016813, -0.015278,
      -0.018165, -0.015592, -0.011380, -0.013111, -0.022270, -0.017536, -0.013861, -0.014942,
      -0.017207, -0.015189, -0.016668, -0.017934, -0.016675, -0.012332, -0.015654, -0.018690,
      -0.018994, -0.015038, -0.014203, -0.013576, -0.018422, -0.014870, -0.013357, -0.012102,
      -0.017044, -0.014259, -0.014130, -0.015934, -0.016579, -0.011408, -0.016220, -0.016180,
      -0.014557, -0.013256, -0.013172, -0.012467, -0.017784, -0.014550, -0.009543, -0.011581,
      -0.013930, -0.011783, -0.013312, -0.013833, -0.011930, -0.009453, -0.012158, -0.009346,
      -0.013235, -0.009963, -0.007856, -0.009778, -0.012328, -0.012472, -0.008635, -0.009341,
      -0.010138, -0.009094, -0.008943, -0.012399, -0.008973, -0.006288, -0.011268, -0.010769,
      -0.012440, -0.009487, -0.006669, -0.005229, -0.012507, -0.011777, -0.006988, -0.008114,
      -0.010143, -0.007369, -0.011195, -0.009766, -0.007376, -0.007739, -0.010248, -0.009106,
      -0.009326, -0.008249, -0.006518, -0.006820, -0.012849, -0.009346, -0.005095, -0.008204,
      -0.008558, -0.008910, -0.009358, -0.008265, -0.007382, -0.005218, -0.009783, -0.007627,
      -0.009001, -0.008573, -0.005851, -0.005610, -0.012132, -0.008854, -0.008714, -0.007548,
      -0.008328, -0.007526, -0.008988, -0.008153, -0.004222, -0.005246, -0.007263, -0.006462,
      -0.007684, -0.003134, -0.003364, -0.001465, -0.007718, -0.006501, -0.003028, -0.000832,
      -0.001651, -0.004787, -0.003745, -0.004260, -0.005388, -0.002093, -0.005560, -0.005330,
      -0.004923, -0.004053, -0.001913, -0.001594, -0.008401, -0.006439, -0.001600, -0.001790,
      -0.003097, -0.000642, -0.003521, -0.004081
    }
  , 1_ns              // sampleDuration
  , 7e6f              // gain
}; // icarus::opdet::SampledWaveformFunction<T>::Waveform


// -----------------------------------------------------------------------------
template <typename T>
icarus::opdet::SampledWaveformFunction<T>::SampledWaveformFunction
  (Time peakTime, float gain)
  : fSource    { &Waveform }
  , fSamples   { buildSamples(*fSource, gain) }
  , fGain      { gain }
  , fPeakSample{ findPeak(fSamples) }
  , fRefTime   { peakTime - fPeakSample * sampleDuration() }
{
  assert(fSource); // way too late, if not true it would have already crashed...
  
} // icarus::opdet::SampledWaveformFunction<>::SampledWaveformFunction()


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::SampledWaveformFunction<T>::doEvaluateAt(Time time) const
  -> ADCcount
{
  std::ptrdiff_t const iSample = static_cast<std::ptrdiff_t>
    (std::floor((time - fRefTime)/sampleDuration()));
  return hasSample(iSample)? fSamples[iSample]: ADCcount{ 0 };
} // icarus::opdet::SampledWaveformFunction<>::doEvaluateAt()


// -----------------------------------------------------------------------------
template <typename T>
void icarus::opdet::SampledWaveformFunction<T>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  
  out << firstIndent
      << "Pulse: " << fSource->name
    << "\n" << indent
      << "  from " << fSamples.size() << "x " << sampleDuration() << " samples"
      << ", peak at " << Base_t::peakTime()
      << " with amplitude " << Base_t::peakAmplitude()
    << "\n" << indent
      << "  start at " << fRefTime << ", gain " << fGain
      << " (integral: " << integral() << ")"
    << '\n'
    ;
  
  // === BEGIN DEBUG FIXME DELME ===============================================
  out << "\n" << indent << "  samples:";
  for (std::size_t iSample = 0; iSample < fSamples.size(); ++iSample) {
    out << "\n" << indent << "   [" << iSample << "]  " << fSamples[iSample];
  }
  out << "\n";
  // === END DEBUG =============================================================
  
} // icarus::opdet::SampledWaveformFunction<>::doDump()


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::SampledWaveformFunction<T>::integral() const -> ADCcount
  { return std::reduce(fSamples.begin(), fSamples.end()); }


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::SampledWaveformFunction<T>::buildSamples
  (WaveformSpecs_t const& waveform, float targetGain) const
  -> std::vector<ADCcount>
{
  
  /*
   * Sample conversion
   * ------------------
   * 
   * The waveform is expected in millivolt, and it expresses the response with/
   * the photodetector set at a known gain.
   * Our target is a waveform in ADC and the target gain in argument.
   * The conversion factor is based on the full range of the digitizer,
   * that is 2 V across 14 bit.
   */
  constexpr float VoltageRange = 2'000.0; // millivolt
  constexpr unsigned short int ADCbits = 14;
  
  // 2 V in 14 bits (=> 8.192):
  constexpr float mVtoADC = (1 << ADCbits) / VoltageRange;
  
  float const factor = (targetGain / waveform.gain) * mVtoADC;
  
  auto voltageToADC = [factor](float mV)
    { return static_cast<ADCcount>(factor * mV); };
  
  std::vector<ADCcount> samples;
  samples.reserve(waveform.samples.size());
  std::transform(waveform.samples.begin(), waveform.samples.end(),
    back_inserter(samples), voltageToADC);
  
  return samples;
  
} // icarus::opdet::SampledWaveformFunction<>::buildSamples()


// -----------------------------------------------------------------------------
template <typename T>
std::size_t icarus::opdet::SampledWaveformFunction<T>::findPeak
  (std::vector<ADCcount> const& samples)
{
  assert(!samples.empty());
  auto const sbegin = samples.begin();
  auto const [ min, max ] = std::minmax_element(sbegin, samples.end());
  // assume baseline 0:
  return ((min->abs() > max->abs())? min: max) - sbegin;
} // icarus::opdet::SampledWaveformFunction<T>::findPeak()


// -----------------------------------------------------------------------------

#endif //  ICARUSCODE_PMT_ALGORITHMS_SAMPLEDWAVEFORMFUNCTION_H
