/**
 * @file   icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.tcc
 * @brief  Algorithm to produce trigger gates out of optical readout waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h`
 *         `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.cxx`
 * 
 */


#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_TCC
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_TCC

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_H
# error "ManagedTriggerGateBuilder.tcc must not be included directly."\
        " #include \"icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h\" instead."
#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_H

#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_EXTRADEBUG 0


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // icarus::trigger::ADCCounts_t
#include "icarusalg/Utilities/WaveformOperations.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/CoreUtils/counter.h"
// #include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h" // MF_LOG_TRACE()

// range library
#include "range/v3/view/chunk_by.hpp"

// C/C++ standard libraries
#include <optional>
#include <iterator> // std::next(), std::prev()
#include <cmath> // std::round()
#include <cstddef> // std::ptrdiff_t


//------------------------------------------------------------------------------
namespace icarus::trigger::details {

  /// An interval between thresholds; may be open (no upper or lower threshold).
  struct ThresholdsBand {
    
    /// Type of a sequence of thresholds, lowest to highest.
    using ThresholdCollection_t = std::vector<ADCCounts_t>;
    
    /// An iterator to one of the thresholds in a sequence.
    using ThresholdIterPtr_t
      = std::optional<ThresholdCollection_t::const_iterator>;
    
    /// An optional threshold (there may be no threshold at all).
    struct Threshold {
      
      ThresholdIterPtr_t thr = std::nullopt;
      
      /// Returns whether this threshold is set.
      bool hasThreshold() const { return thr.has_value(); }
      
      /// Returns the value of the threshold (unchecked).
      ADCCounts_t threshold() const { return **thr; }
      
      /// Rendering to string of the value of the threshold, if any.
      operator std::string() const;
      
      /// Set to no threshold.
      void remove() { thr = std::nullopt; }
      
      /// Jump to the higher threshold (unchecked).
      void goHigher() { ++*thr; }
      
      /// Jump to the lower threshold (unchecked).
      void goLower() { --*thr; }
      
      /// Returns whether `sample` is lower than the threshold.
      bool sampleLower(ADCCounts_t sample) const
        { return hasThreshold() && lower(sample); }
      
      /// Returns whether `sample` is not lower than the threshold.
      bool sampleHigher(ADCCounts_t sample) const
        { return hasThreshold() && !lower(sample); }
      
        private:
      
      /// Unchecked threshold comparison for `sampleLower()`.
      bool lower(ADCCounts_t sample) const { return sample < **thr; }
      
    }; // Threshold
    
    ThresholdIterPtr_t const bottom; ///< Iterator to the lowest threshold.
    ThresholdIterPtr_t const top; ///< Iterator to the highest threshold.
    
    Threshold lower; /// Current lower threshold.
    Threshold upper; /// Current upper threshold.
    
    /// Constructor: spans the lowest and highest `thresholds`, starts from low.
    ThresholdsBand(ThresholdCollection_t const& thresholds);
    
    /// Returns whether `sample` is lower than the current lower threshold.
    bool lowerThanLower(ADCCounts_t sample) const;
    
    /// Returns whether `sample` is higher than the current upper threshold.
    bool higherThanUpper(ADCCounts_t sample) const;
    
    /// Moves the current band to the next lower threshold available.
    bool stepLower();
    
    /// Moves the current band to the next higher threshold available.
    bool stepHigher();
    
  }; // ThresholdsBand

} // namespace icarus::trigger::details


//------------------------------------------------------------------------------
//--- icarus::trigger::ManagedTriggerGateBuilder
//------------------------------------------------------------------------------
template <typename GateMgr>
auto icarus::trigger::ManagedTriggerGateBuilder::unifiedBuild
  (GateMgr&& gateManager, std::vector<WaveformWithBaseline> const& waveforms)
  const -> std::vector<TriggerGates>
{
  using GateManager_t = GateMgr;
  using GateInfo_t = typename GateManager_t::GateInfo_t;
  
  /*
   * This is the simple algorithm where each channel is treated independently,
   * and we have as many trigger gates as we have channels.
   */
  
  // create an empty TriggerGates object for each threshold;
  // thresholds are kept relative
  std::vector<TriggerGates> allGates = prepareAllGates();
  
  raw::Channel_t channel = InvalidChannel;
  
  // now group the waveforms by channel (must be already sorted!)
  // and process waveforms channel by channel
  auto sameChannel
    = [] (WaveformWithBaseline const& a, WaveformWithBaseline const& b)
      { return a.waveform().ChannelNumber() == b.waveform().ChannelNumber(); }
    ;
  
  auto byChannel = waveforms | ranges::views::chunk_by(sameChannel);
  for (auto const& channelWaveforms: byChannel) {
    
    auto const& firstWaveform = channelWaveforms.front().waveform();
    
    // assert that the waveforms are sorted by channel and then by time
    // and not overlapping
    assert(
      !isValidChannel(channel)
      || (firstWaveform.ChannelNumber() >= channel)
      );
    if (firstWaveform.ChannelNumber() != channel)
      channel = firstWaveform.ChannelNumber();
    
    // we don't know how many... (maybe C++20 ranges will tell us)
    MF_LOG_TRACE(details::TriggerGateDebugLog)
      << "Building trigger gates from waveforms on channel " << channel;
    
    std::vector<GateInfo_t> channelGates;
    channelGates.reserve(nChannelThresholds());
    for (TriggerGates& thrGates: allGates) {
      channelGates.push_back 
        (gateManager.create(thrGates.gateFor(firstWaveform)));
    }
    
    // this method will update the channel gates referenced in `channelGates`,
    // which are owned by `allGates`
    buildChannelGates(channelGates, channelWaveforms);
    
  } // for channels
  
  return allGates;
} // icarus::trigger::ManagedTriggerGateBuilder::unifiedBuild()


//------------------------------------------------------------------------------
template <typename GateInfo, typename Waveforms>
void icarus::trigger::ManagedTriggerGateBuilder::buildChannelGates(
  std::vector<GateInfo>& channelGates,
  Waveforms const& channelWaveforms
) const
{
  using ops = icarus::waveform_operations::NegativePolarityOperations<float>;
  
  if (channelWaveforms.empty()) return;
  
  //
  // extract all thresholds from each waveform in one pass
  //
  
  raw::OpDetWaveform const& firstWaveform
    = channelWaveforms.front().waveform();
  raw::Channel_t const channel = firstWaveform.ChannelNumber();
  
  // used only in debug mode:
  optical_tick lastWaveformTick [[gnu::unused]]
    = timeStampToOpticalTick(firstWaveform.TimeStamp());
  
  /*
   * The algorithm finds gate openings and closing.
   * The actual actions on opening and closing depends on the gate info class.
   * For example, while a dynamic gate duration algorithm will perform open and
   * close operations directly, a fixed gate duration algorithm may perform
   * both opening and closing at open time, and nothing at all at closing time.
   */
  unsigned int nWaveforms = 0U;
  for (auto const& waveformData: channelWaveforms) {

    raw::OpDetWaveform const& waveform = waveformData.waveform();
    
    ops const waveOps { waveformData.baseline().baseline() };
    
    // baseline subtraction is performed in floating point,
    // but then rounding is applied again
    auto subtractBaseline = [waveOps](float sample) -> ADCCounts_t
      {
        return
          ADCCounts_t::castFrom(std::round(waveOps.subtractBaseline(sample)));
      };
    
    
    ++nWaveforms;
    assert(waveform.ChannelNumber() == channel);
    
    // start of the waveform (tick #0) in electronics time scale
    // and in optical tick units
    optical_tick const waveformTickStart
      = timeStampToOpticalTick(waveform.TimeStamp());
    optical_tick const waveformTickEnd
      = waveformTickStart + optical_time_ticks::castFrom(waveform.size());
    
    MF_LOG_TRACE(details::TriggerGateDebugLog)
      << "Waveform with " << waveform.size() << " ticks"
      <<  ": [ " << waveformTickStart << " -- " << waveformTickEnd << " ]";
    
    // assert that the waveforms are sorted by channel and then by time
    // and not overlapping
    assert(lastWaveformTick <= waveformTickStart);
    lastWaveformTick = waveformTickEnd;
    
    // register this waveform with the gates (this feature is unused here)
    for (auto& gateInfo: channelGates) gateInfo.addTrackingInfo(waveform);
    
    // all gates start closed; this gate is not necessarily closed, but the
    // waveform is not above the gate threshold any more.
    auto nextGateToOpen = channelGates.begin();
    
    // we keep track of whether we have no lower or higher thresholds available
    // to simplify the checks;
    // we name them "pp" because they behave (almost) like pointers to pointers
    
    // start at bottom with no lower threshold:
    details::ThresholdsBand thresholds{ channelThresholds() };
    
    std::ptrdiff_t const bStart = fSampleOffset;
    std::ptrdiff_t const bStep = fBlockSize;
    std::ptrdiff_t const bEnd = waveform.size();
    for (std::ptrdiff_t iBStart = bStart; iBStart < bEnd; iBStart += bStep) {
      
      // TODO check the time alignment of the waveform
      // TODO make sure we don't exceed the size of the waveform
      
      #if ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_EXTRADEBUG
      // this is too much also for regular debugging...
      MF_LOG_TRACE(details::TriggerGateDebugLog)
        << "block at sample #" << iBStart << ":";
      #endif
      
      ADCCounts_t blockRelSample = std::numeric_limits<ADCCounts_t>::min();
      for (std::ptrdiff_t sampleOffset: fPatternIndices) {
        std::ptrdiff_t iSample = iBStart + sampleOffset;
      
        // baseline subtraction is always a subtraction (as in "A minus B"),
        // regardless the polarity of the waveform
        auto const sample = waveform[iSample];
        ADCCounts_t const relSample = subtractBaseline(sample);
        blockRelSample = std::max(relSample, blockRelSample);
        
        #if ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_EXTRADEBUG
        // this is too much also for regular debugging...
        MF_LOG_TRACE(details::TriggerGateDebugLog)
          << "  sample +" << sampleOffset << ": " << sample << " [=> "
          << relSample << "]";
        #endif
        
      } // for enabled sample in block
      
      #if ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_EXTRADEBUG
      // this is too much also for regular debugging...
      MF_LOG_TRACE(details::TriggerGateDebugLog)
        << "block level: " << blockRelSample << ", thresholds lower: "
          << std::string(thresholds.lower)
        << ", upper: " << std::string(thresholds.upper) << " [*]"
        ;
      #endif
      
      auto const blockTimeTick
        = optical_time_ticks{ iBStart + fBlockTimeReference };
      
      //
      // if the level of this block is lower than the current lower threshold,
      // we are just tracking the thresholds: gate closing has already happened
      //
      if (thresholds.lowerThanLower(blockRelSample)) {
        
        MF_LOG_TRACE(details::TriggerGateDebugLog)
          << "Block " << iBStart << " (" << blockRelSample << " on "
          << waveOps.baseline() << ") leaving thresholds at "
          << waveformTickStart << " + " << blockTimeTick;
        
        do { // we keep opening gates at increasing thresholds
          
          // if there is an lower threshold, there must also be a open gate!
          assert(nextGateToOpen != channelGates.begin());
          
          
          (--nextGateToOpen)
            ->belowThresholdAt(waveformTickStart + blockTimeTick);
          
          MF_LOG_TRACE(details::TriggerGateDebugLog)
            << "  => decreasing threshold " << thresholds.lower.threshold();
          
          thresholds.stepLower();
          
        } while (thresholds.lowerThanLower(blockRelSample));
        
      } // if closing gate
      
      //
      // if this sample is greater or matching the next threshold,
      // we *are* opening gate(s)
      //
      else if (thresholds.higherThanUpper(blockRelSample)) {
        
        MF_LOG_TRACE(details::TriggerGateDebugLog)
          << "Block " << iBStart << " (" << blockRelSample << " on "
          << waveOps.baseline() << ") passing thresholds at "
          << waveformTickStart << " + " << blockTimeTick;
        
        do { // we keep opening gates at increasing thresholds
          
          // note that it is not guaranteed that gates at lower thresholds are
          // still open (that depends on the builder implementation)
          
          // if there is an upper threshold, there must also be a closed gate!
          assert(nextGateToOpen != channelGates.end());
          
          MF_LOG_TRACE(details::TriggerGateDebugLog)
            << "Opening thr=" << thresholds.upper.threshold();
          
          thresholds.stepHigher();
          
          (nextGateToOpen++)
            ->aboveThresholdAt(waveformTickStart + blockTimeTick);
          
        } while (thresholds.higherThanUpper(blockRelSample));
        
      } // if opening gate
      
    } // for block

  } // for waveforms
  
  MF_LOG_TRACE(details::TriggerGateDebugLog)
    << "Trigger gates from " << nWaveforms << " waveforms on channel "
    << channel << " completed.";
  
} // icarus::trigger::ManagedTriggerGateBuilder::buildChannelGates()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_TCC
