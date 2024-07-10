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
