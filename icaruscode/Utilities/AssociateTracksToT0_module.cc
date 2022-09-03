/**
 * @file   AssociateTracksToT0_module.cc
 * @brief  Creates a direct association between tracks and time.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 31, 2022
 *
 */


// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <memory> // std::make_unique()
#include <vector>
#include <string>
#include <utility> // std::move()
#include <cassert>


//------------------------------------------------------------------------------
/**
 * @brief Creates a direct association between `recob::Track` and `anab::T0`.
 *
 * Pandora pattern recognition chose to associate time to a particle flow object
 * (PFO) rather than to the track reconstructed from it. Nevertheless, it may be
 * convenient to have a direct association, for example in the case we want to
 * use algorithms or modules that deal with tracks and not PFO. Time tagging via
 * optical or scintillation detectors may for example not go through PFO at all.
 * 
 * This module creates a direct association between `recob::Track` and
 * `anab::T0` objects that share a connection with the same `recob::PFParticle`.
 * The connection discovery is purely mechanic, relying on the common _art_
 * pointers in the two associations.
 * 
 * The module requires as input:
 *  * an association between tracks and particle flow objects
 *  * an association between time and particle flow objects
 * 
 * Note that the user algorithm needs to be made aware of the input tag for
 * the newly created associations, which is not the same as any of the inputs.
 * No new collection (neither of tracks nor of times) is created.
 * 
 *
 * Output
 * =======
 * 
 * * An association `art::Assns<recob::Track, anab::T0>` (instance name: empty
 *   by default) between the tracks and their time.
 * 
 * 
 * Input data products
 * ====================
 * 
 * This module uses as input:
 *  * An association of particle flow objects to tracks
 *    (`art::Assns<recob::PFParticle, recob::Track>`).
 *  * An association of particle flow objects to tracks
 *    (`art::Assns<recob::PFParticle, anab::T0>`).
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 *  * `PFOTrackAssnTag` (input tag, mandatory): data product with the
 *    associations between particle flow objects and their tracks.
 *  * `PFOT0AssnTag` (input tag, default: same as `PFOTrackAssnTag`): data
 *    product with the associations between particle flow objects and time.
 *  * `InstanceName` (string, default: empty): instance name of
 *    the produced associations between the tracks and the times.
 *  * `LogCategory` (string, default: `"AssociateTracksToT0"`): name of the
 *    message facility stream to use for the output to console.
 *
 */
class AssociateTracksToT0: public art::SharedProducer {
  
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> PFOTrackAssnTag {
      Name{ "PFOTrackAssnTag" },
      Comment{ "input association between particle flow objects and tracks" }
      };
    
    fhicl::OptionalAtom<art::InputTag> PFOT0AssnTag {
      Name{ "PFOT0AssnTag" },
      Comment{
        "input association between particle flow objects and times"
        " [default: as PFOTrackAssnTag]"
        }
      };
    
    fhicl::Atom<std::string> InstanceName {
      Name{ "InstanceName" },
      Comment{
        "instance name of produced associations between the tracks and times"
        },
        ""
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of message facility category for console message stream" },
      "AssociateTracksToT0"
      };
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  
  AssociateTracksToT0(Parameters const& params, const art::ProcessingFrame&);
  
  
  virtual void produce(art::Event& event, const art::ProcessingFrame&) override;
  
  
    private:
  
  // --- BEGIN --  Configuration  ----------------------------------------------
  
  art::InputTag const fPFOTrackTag; ///< The PFO-track associations for input.
  art::InputTag const fPFOT0Tag;    ///< The PFO-time associations for input.
  
  /// Instance name for the association of tracks and times.
  std::string const fInstanceName;
  
  std::string const fLogCategory; ///< Category for messagefacility stream.
  
  // --- END ----  Configuration  ----------------------------------------------
  
  
}; // class AssociateTracksToT0


//------------------------------------------------------------------------------
//---  implementation
//------------------------------------------------------------------------------
AssociateTracksToT0::AssociateTracksToT0
  (Parameters const& params, const art::ProcessingFrame&)
  : art::SharedProducer{ params }
  // configuration parameters
  , fPFOTrackTag { params().PFOTrackAssnTag() }
  , fPFOT0Tag    { params().PFOT0AssnTag().value_or(fPFOTrackTag) }
  , fInstanceName{ params().InstanceName() }
  , fLogCategory { params().LogCategory() }
{
  async<art::InEvent>();
  
  //
  // declaration of input
  //
  consumes<art::Assns<recob::PFParticle, recob::Track>>(fPFOTrackTag);
  consumes<art::Assns<recob::PFParticle, anab::T0>>(fPFOT0Tag);
  
  //
  // declaration of output
  //
  produces<art::Assns<recob::Track, anab::T0>>(fInstanceName);
  
  //
  // configuration dump
  //
  {
    mf::LogInfo log{ fLogCategory };
    log << "Joining tracks from '" << fPFOT0Tag.encode()
      << "' and times from '" << fPFOT0Tag.encode() << "' via shared particle flow objects.";
  }
  
} // AssociateTracksToT0::AssociateTracksToT0()


//------------------------------------------------------------------------------
void AssociateTracksToT0::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  /*
   * The algorithm goes through all the PFO associated with a T0, and uses
   * `art::FindManyP` to find the corresponding tracks to be associated to each
   * `T0`. Null pointers are skipped.
   */
  
  // read inputs
  auto const& PFOtoT0
    = event.getProduct<art::Assns<recob::PFParticle, anab::T0>>(fPFOT0Tag);
  
  std::vector<art::Ptr<recob::PFParticle>> PFOlist;
  PFOlist.reserve(PFOtoT0.size());
  for (auto const& [ PFOptr, T0ptr ]: PFOtoT0)
    if (PFOptr && T0ptr) PFOlist.push_back(PFOptr);
  art::FindManyP<recob::Track> const PFOtoTrack{ PFOlist, event, fPFOTrackTag };
  
  { // --- BEGIN -- DEBUG
    
    mf::LogTrace log{ fLogCategory };
    log << "'" << fPFOT0Tag.encode() << "': " << PFOtoT0.size()
      << " PFO-T0 associations";
    for (auto const& [ PFOptr, T0ptr ]: PFOtoT0)
      log << "\n - " << PFOptr << " <=> " << T0ptr;
    
    auto const& PFOtoTrack
      = event.getProduct<art::Assns<recob::PFParticle, recob::Track>>
        (fPFOTrackTag)
      ;
    log << "\n'" << fPFOTrackTag.encode() << "': " << PFOtoTrack.size()
      << " PFO-T0 associations";
    for (auto const& [ PFOptr, trackPtr ]: PFOtoTrack)
      log << "\n - " << PFOptr << " <=> " << trackPtr;
    
  } // --- END ---- DEBUG
  
  
  // join associations
  mf::LogTrace{ fLogCategory } << "Producing track-T0 associations:";
  art::Assns<recob::Track, anab::T0> trackToTime;
  for (auto const& [ iPFO, PFOptr ]: util::enumerate(PFOlist)) {
    
    auto const& thisPFOtoT0 = PFOtoT0.at(iPFO); // one T0
    auto const& T0ptr = thisPFOtoT0.second;
    assert(PFOptr == thisPFOtoT0.first);
    
    for (art::Ptr<recob::Track> const& trackPtr: PFOtoTrack.at(iPFO)) {
      if (!trackPtr) continue;
      trackToTime.addSingle(trackPtr, T0ptr);
      mf::LogTrace{ fLogCategory } << " - " << trackPtr << " <=> " << T0ptr;
    }
  } // for
  mf::LogTrace{ fLogCategory }
    << "Produced " << trackToTime.size() << " track-T0 associations.";
  
  
  // save the data products
  event.put(
    std::make_unique<art::Assns<recob::Track, anab::T0>>
      (std::move(trackToTime)),
    fInstanceName
    );
  
} // AssociateTracksToT0::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(AssociateTracksToT0)


//------------------------------------------------------------------------------
