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

// C/C++ standard libraries
#include <algorithm> // std::max()


//------------------------------------------------------------------------------
//--- ManagedTriggerGateBuilder implementation
//------------------------------------------------------------------------------
icarus::trigger::ManagedTriggerGateBuilder::ManagedTriggerGateBuilder
  (Config const& config)
  : Base_t(config)
  , fSamplePrescale(std::max(config.SamplePrescale(), std::size_t{ 1 }))
  , fSampleOffset(config.SampleOffset())
{}


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
    std::string const& fIndent;
    std::string const& fFirstIndent;
    
    unsigned int fCount = 0;
    
      public:
    IndentFactory(std::string const& indent, std::string const& firstIndent)
      : fIndent(indent), fFirstIndent(firstIndent) {}
    std::string const& operator() ()
      { return (fCount++ == 0)? fFirstIndent: fIndent; }
  }; // IndentFactory
  
  IndentFactory nextIndent{ indent, firstIndent };

  if (fSamplePrescale > 1) {
    out << nextIndent()
      << " * use only one out of " << fSamplePrescale << " samples";
  }
  if (fSampleOffset > 0) {
    out << nextIndent()
      << " * skip the first " << fSampleOffset << " samples of each waveform";
  }
  
} // icarus::trigger::ManagedTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------
