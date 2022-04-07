/**
 * @file   icaruscode/PMT/SampledWaveformFunctionTool_tool.cc
 * @brief  Toolization of `icarus::opdet::SampledWaveformFunction<nanosecond>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/SampledWaveformFunction.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::SinglePhotonPulseFunctionTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/SampledWaveformFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electromagnetism.h" // picocoulomb
#include "lardataalg/Utilities/quantities_fhicl.h" // nanoseconds from FHiCL

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Utilities/Exception.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::opdet { struct SampledWaveformFunctionTool; }
/**
 * @brief Creates a `SampledWaveformFunction` pulse shape.
 * @see `icarus::opdet::SinglePhotonPulseFunctionTool`
 * 
 * This tool creates a `icarus::opdet::SampledWaveformFunction<nanosecond>`
 * function to describe a R5912 PMT pulse attached to the ICARUS detector.
 * 
 * See `icarus::opdet::SampledWaveformFunction` for the details of the function.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description SampledWaveformFunctionTool` (or read `Config`
 * data structure) for a short explanation of the meaning of the parameters.
 * 
 * The amplitude of the function is scaled in terms of PMT gain
 * (`Gain` parameter).
 */
struct icarus::opdet::SampledWaveformFunctionTool
  : public icarus::opdet::SinglePhotonPulseFunctionTool
{
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<nanoseconds> TransitTime {
      Name("TransitTime"),
      Comment("peak time from the beginning of the waveform [ns]")
      // mandatory
      };
    fhicl::Atom<float> Gain {
      Name("Gain"),
      Comment("PMT amplification gain factor")
      };
    
  }; // struct Config

  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<Config>;
  
  /// Constructor: sets the configuration.
  SampledWaveformFunctionTool(Parameters const& config)
    : fPulseFunction(makePulseFunction(config())) {}
  
  
    private:
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /// Returns the function that was created at construction time.
  virtual std::unique_ptr<PulseFunction_t> doGetPulseFunction() override
    { return std::move(fPulseFunction); }
  
  // --- END -- Virtual interface ----------------------------------------------
  
  /// Function stored while waiting to be delivered.
  std::unique_ptr<PulseFunction_t> fPulseFunction;
  
  
  /// Creates and returns a pulse function with the specified configuration.
  static std::unique_ptr<PulseFunction_t> makePulseFunction
    (Config const& config);
  
}; // icarus::opdet::SampledWaveformFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::SampledWaveformFunctionTool implementation
//------------------------------------------------------------------------------
auto icarus::opdet::SampledWaveformFunctionTool::makePulseFunction
  (Config const& config) -> std::unique_ptr<PulseFunction_t>
{
  
  using MyFunction_t = icarus::opdet::SampledWaveformFunction<nanoseconds>;
  return std::make_unique<MyFunction_t>(
      config.TransitTime()  // peakTime
    , config.Gain()         // gain
    );
  
} // icarus::opdet::SampledWaveformFunctionTool::makePulseFunction()


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::SampledWaveformFunctionTool)


//------------------------------------------------------------------------------

