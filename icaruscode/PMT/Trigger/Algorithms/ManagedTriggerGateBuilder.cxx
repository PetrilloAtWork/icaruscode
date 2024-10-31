/**
 * @file   icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.cxx
 * @brief  Algorithm to produce trigger gates out of optical readout waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.tcc`
 *         `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h`
 * 
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// C/C++ standard libraries
#include <algorithm> // std::max()


// -----------------------------------------------------------------------------
// ---  icarus::trigger::details::ThresholdsBand implementation
// -----------------------------------------------------------------------------
icarus::trigger::details::ThresholdsBand::Threshold::operator std::string() const
  { return (hasThreshold()? util::to_string(threshold()): "none"); }


// TODO different direction may require different starting point (maybe not though)
icarus::trigger::details::ThresholdsBand::ThresholdsBand
  (std::vector<ADCCounts_t> const& thresholds)
  : bottom{ thresholds.begin() }
  , top{ thresholds.end() }
  , lower{ std::nullopt }
  , upper{ thresholds.empty()? std::nullopt: bottom }
{}


bool icarus::trigger::details::ThresholdsBand::lowerThanLower
  (ADCCounts_t sample) const
  { return lower.sampleLower(sample); }


bool icarus::trigger::details::ThresholdsBand::higherThanUpper
  (ADCCounts_t sample) const
  { return upper.sampleHigher(sample); }


bool icarus::trigger::details::ThresholdsBand::stepLower() {
  if (!lower.hasThreshold()) return false;
  // was this the bottom threshold?
  upper = lower;
  if (lower.thr == bottom) {
    lower.remove();
    return false;
  }
  lower.goLower(); // point to previous threshold
  return true;
} // icarus::trigger::details::ThresholdsBand::stepLower()


bool icarus::trigger::details::ThresholdsBand::stepHigher() {
  if (!upper.hasThreshold()) return false; // ?!
  lower = upper;
  upper.goHigher(); // point to next threshold
  // was it the top threshold?
  if (upper.thr != top) return true;
  upper.remove();
  return false;
} // icarus::trigger::details::ThresholdsBand::stepHigher()


//------------------------------------------------------------------------------
//--- ManagedTriggerGateBuilder implementation
//------------------------------------------------------------------------------
icarus::trigger::ManagedTriggerGateBuilder::ManagedTriggerGateBuilder
  (Config const& config)
  : Base_t(config)
  , fPatternIndices(patternToIndices(config.SamplingPattern()))
  , fBlockSize(config.SamplingPattern().size())
  , fBlockTimeReference(config.BlockTimeReference())
  , fSampleOffset(config.SampleOffset())
{
  if (fBlockSize == 0) {
    throw cet::exception{ "ManagedTriggerGateBuilder" }
      << "Configuration error: '" << config.SamplingPattern.name()
      << "' has no samples.\n";
  }
  
  if (fPatternIndices.empty()) {
    throw cet::exception{ "ManagedTriggerGateBuilder" }
      << "Configuration error: '" << config.SamplingPattern.name()
      << "' specifies a pattern with no active sample.\n";
  }
  
}


//------------------------------------------------------------------------------
void icarus::trigger::ManagedTriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  Base_t::doDumpConfiguration(out, indent, firstIndent);
  dumpLocalConfiguration(out, indent, "\n" + indent);
  
} // icarus::trigger::ManagedTriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::ManagedTriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  class IndentFactory {
    std::string const fNewLineIndent;
    std::string const& fFirstIndent;
    
    unsigned int fCount = 0;
    
      public:
    IndentFactory(std::string const& indent, std::string const& firstIndent)
      : fNewLineIndent('\n' + indent), fFirstIndent(firstIndent) {}
    std::string const& operator() ()
      { return (fCount++ == 0)? fFirstIndent: fNewLineIndent; }
  }; // IndentFactory
  
  IndentFactory newLine{ indent, firstIndent };
  
  out << newLine()
    << " * state is decided by blocks of " << fBlockSize << " samples";
  if (fBlockSize != fPatternIndices.size()) {
    out << newLine()
      << " * only " << fPatternIndices.size()
      << " samples will be tested in each block:";
    for (std::ptrdiff_t const index: fPatternIndices)
      out << " " << index;
  }
  if ((fBlockSize != 1) || (fBlockTimeReference != 0)) {
    out << newLine() << " *";
    if (fBlockTimeReference == 0) out << " the first sample";
    else out << " sample #" << fBlockTimeReference << " from the start";
    out << " of each block will be the block time reference";
  }
  
  if (fSampleOffset > 0) {
    out << newLine()
      << " * skip the first " << fSampleOffset << " samples of each waveform";
  }
  
} // icarus::trigger::ManagedTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------
std::vector<std::ptrdiff_t>
icarus::trigger::ManagedTriggerGateBuilder::patternToIndices
  (std::vector<bool> const& pattern)
{
  std::vector<std::ptrdiff_t> indices;
  for (auto const& [ index, enabled ]: util::enumerate(pattern)) {
    if (enabled) indices.push_back(index);
  }
  return indices;
} // icarus::trigger::ManagedTriggerGateBuilder::patternToIndices()


//------------------------------------------------------------------------------
