/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h
 * @brief  Applies sliding window trigger patterns.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWPATTERNALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWPATTERNALG_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <limits> // std::numeric_limits<>
#include <cstddef> // std::size_t

/**
 * @brief Applies sliding window patterns to discriminated waveforms.
 * 
 * This algorithm takes as input one trigger gate for each window.
 * Each window is identified by an index, and the input gates are assigned one
 * per window, in order.
 * 
 * The window mapping that defines the topology of the windows is passed to the
 * algorithm, as the beam gate definition is.
 * 
 * For the definition of the windows, see
 * `icarus::trigger::details::WindowChannelMap`.
 * 
 */
class icarus::trigger::SlidingWindowPatternAlg {
  
  /// Record of the trigger response.
  using TriggerInfo_t = icarus::trigger::details::TriggerInfo_t;
  
  /// Type of trigger gate extracted from the input event.
  using InputTriggerGate_t = icarus::trigger::MultiChannelOpticalTriggerGate;
  
  /// A list of trigger gates from input.
  using TriggerGates_t = std::vector<InputTriggerGate_t>;

  /// Type of gate data without channel information.
  using TriggerGateData_t = InputTriggerGate_t::GateData_t;
  
  /// Type holding information about composition and topology of all windows.
  using WindowTopology_t = icarus::trigger::details::WindowChannelMap;
  
  /// Additional information on the trigger.
  struct MoreInfo_t {
    
    /// Index of the window which led the trigger.
    std::size_t windowIndex = std::numeric_limits<std::size_t>::max();
    
  }; // struct MoreInfo_t
  
  
  /**
   * @brief Constructor: configures window topology and times.
   * @param windowTopology full composition and topology description of windows
   * @param beamGate object applying the beam gate to trigger gates
   * 
   * Window topology can be computed with TODO
   * 
   * A beam gate application class can be constructed via
   * `icarus::trigger::makeApplyBeamGate()` helper function.
   * 
   */
  icarus::trigger::SlidingWindowPatternAlg(
    WindowTopology_t windowTopology,
    icarus::trigger::ApplyBeamGateClass beamGate
    );
  
  /// Returns the trigger response from the specified `gates`.
  std::pair<TriggerInfo_t, MoreInfo_t> simulateResponse
    (TriggerGates_t const& gates) const;
  

  /// Returns a new collection of gates, set each in coincidence with beam gate.
  TriggerGates_t applyBeamGate(TriggerGates_t const& gates) const;
  
  
  
  /// Data structure to communicate internally a trigger response.
  struct WindowTriggerInfo_t {
    
    std::size_t windowIndex = std::numeric_limits<std::size_t>::max();
    TriggerInfo_t info;
    
    bool fired() const { return info.fired(); }
    
    operator bool() const { return bool(info); }
    bool operator! () const { return !info; }
    
    void emplace(std::size_t index, TriggerInfo_t info)
      { windowIndex = index; this->info = std::move(info); }
    
  }; // WindowTriggerInfo_t
  
  /// Definition of the neighborhood of each window in terms of window indices.
  WindowTopology_t const fWindowTopology;
  
  
  /// Time interval when to evaluate the trigger.
  icarus::trigger::ApplyBeamGateClass const fBeamGate;
  
  std::string const fLogCategory; ///< Message category tag.
  
  /**
   * @brief Chceks `gates` are compatible with the current window configuration.
   * @param gates the combined sliding window trigger gates, per cryostat
   * @throw cet::exception (category: `SlidingWindowTriggerEfficiencyPlots`)
   *        or derived, if an incompatibility is found
   * 
   * The method verifies that the current channel mapping is compatible with the
   * gates.
   * 
   * This currently means that the `gates` are in the expected order and have
   * the expected channel content.
   */
  void verifyInputTopology(TriggerGates_t const& gates) const;
  
}; // class icarus::trigger::SlidingWindowPatternAlg


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWPATTERNALG_H
