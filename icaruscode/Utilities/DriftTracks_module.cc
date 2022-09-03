/**
 * @file   DriftTracks_module.cc
 * @brief  Creates a new set of tracks, drifted according to a time.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 31, 2022
 *
 */


// ICARUS libraries
#include "icaruscode/Utilities/TrackShifterAlg.h"
#include "icaruscode/Utilities/ReplicateAssociations.h" // replaceAssnsFirst()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::demangle()
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h" // also art::fill_ptr_vector()
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <memory> // std::make_unique()
#include <unordered_map>
#include <vector>
#include <string>
#include <utility> // std::move()
#include <type_traits> // std::remove_reference_t


//------------------------------------------------------------------------------
/**
 * @brief Creates a new set of tracks, drifted according to a time.
 * @see `icarus::TrackShifterAlg`
 *
 * This module produces a new set of tracks (`recob::Track`) with a trajectory
 * that is shifted in the drift coordinate. The amount of shift is determined
 * from the difference between the track reference time and a new track time.
 * 
 * The module requires as input:
 *  * the collection of tracks to be shifted
 *  * their association with the time used to reconstruct them
 *  * their association with the new time
 * 
 * A new set of tracks is produced, in the same order as the ones in the input
 * track collection. The shifting happens as follows:
 *  * if a track is associated to a new time:
 *      * if the track is associated also to an old time, shift the track by the
 *        time difference between the two, after porting them to the same
 *        reference time; the new track will be associated to the object with
 *        the new time;
 *      * if the track is not associated to an old time, pretend the old time
 *        of the track is the reference of the old times, and shift the track
 *        accordingly; the new track will be associated to the object with the
 *        old time;
 *  * if a track is not associated to a new time:
 *      * if the track is associated to an old time, the track is copied
 *        unshifted and it is associated to the object holding that time;
 *      * if the track is also not associated to an old time, the track is
 *        still copied unshifted, and it is not associated to any time object.
 * 
 * The actual shifting is delegated to the `icarus::TrackShifterAlg` algorithm.
 * 
 * Additional associations may be produced for the new tracks mirroring the
 * ones already existing for the input tracks. The new associations link the
 * new track to exactly the same objects as the old tracks (for example, even
 * for an association with `recob::SpacePoint`, the new association points to
 * the original space point objects, and they are not modified even if they
 * describe the positions before the shifting).
 * 
 * @todo Decide whether to support the case whether only associations are
 *       provided as input, and not the track collection. Tracks would be
 *       identified by their pointers in the two time associations, and be
 *       shifted if present among the new time associations, and just copied if
 *       only present among the old time associations.
 *       The order of the tracks would need to be decided (maybe the same order
 *       as the sorted collection of pointers).
 * 
 * @todo Decide whether to implement an option where the shift is performed only
 *       if the track is _not_ associated to an old time (implying that we trust
 *       that time more). Without this option, a list of new times should be
 *       produced excluding the ones from such tracks, and that list should be
 *       then used as input for this module.
 * 
 * 
 * 
 * Determination of the drift amount
 * ----------------------------------
 * 
 * See the documentation of `icarus::TrackShifterAlg` for an explanation and
 * full derivation of the shift principles.
 * 
 * 
 * Cathode-crossing tracks
 * ------------------------
 * 
 * See the documentation of `icarus::TrackShifterAlg` for important observations
 * about the effects on cathode-crossing tracks.
 *
 *
 * Output
 * =======
 * 
 * * A `std::vector<recob::Track>` collection with all tracks, one per input
 *   track, shifted accordingly to the algorithm criteria.
 * * An association `art::Assns<recob::Track, anab::T0>` (instance name: empty
 *   by default) between the new tracks and the times utilized for the
 *   shifting; the pointers to these time object may come from either the new
 *   time input data product or the old time input data product.
 * * An association `art::Assns<recob::Track, anab::T0>` (instance name:
 *   set by `NewTracksAssnToOldT0` configuration parameter) between the new
 *   tracks and the times utilized for the shifting; the pointers to these time
 *   object may come from either the new time input data product or the old time
 *   input data product.
 * * Additional associations to the new tracks (`recob::Track`) and to the new
 *   times (`anab::T0`) as requested by the configuration.
 * * Additional associations may be produced for the new tracks mirroring the
 *   ones already existing for the input tracks. For the associations including
 *   metadata, the metadata is copied. In order for an association to be created,
 *   the proper configuration parameter must be explicitly specified.
 *   The supported types of associations are listed below.
 *     * With TPC hits, and metadata
 *       (`art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>`)
 *       via parameter `AddAssnsToHitWithMetadataFrom`.
 *     * With particle flow objects
 *       (`art::Assns<recob::Track, recob::PFParticle>`)
 *       via parameter `AddAssnsToPFParticleFrom`.
 * 
 * Note that it's not possible to build a direct _art_ association between the
 * unshifted and the shifted tracks, because `art::Assns` does not support
 * associations between objects of the same type.
 * 
 * 
 * @todo Decide whether to produce a new `anab::T0` collection (and
 *       associations) to the new tracks. Without it, associations of the
 *       shifted tracks to `anab::T0` may span different collections (the new
 *       times when available, the old time if the new is not available), which
 *       may be a complication for the downstream modules.
 *       If we produce that new collection, we need to shift all times to the
 *       reference of the new time.
 * 
 * @todo Maybe add the option to customize the instance name of the associations
 *       to the new time.
 * 
 * @todo Additional associations are a big to-do, and a proper design should be
 *       added to deal with them.
 * 
 * 
 * Input data products
 * ====================
 * 
 * This module uses as input:
 *  * A collection of tracks (`std::vector<recob::Track>`) to be shifted
 *    (`Tracks` configuration parameter).
 *  * An association of tracks from that collection to the times used to place
 *    them in drift direction (`art::Assns<recob::Track, anab::T0>`)
 *    (`TrackTimes` configuration parameter).
 *  * An association of tracks from that collection to the times to be used to
 *    relocate them in drift direction (`art::Assns<recob::Track, anab::T0>`)
 *    (`NewTimes` configuration parameter).
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 *  * `Tracks` (input tag, mandatory): list of the tracks to shift or copy.
 *  * `TrackTimes` (input tag, optional): list of the associations between the
 *    input tracks and the time used to place them in the drift direction.
 *    If not specified, the same tag as `Tracks` is used. If it is specified
 *    empty, it is assumed that all the tracks have been placed using the
 *    "old reference time".
 *  * `NewTimes` (input tag, optional): list of the associations of the input
 *    tracks to their new time.
 *  * `OldTimeReference` (keyword; default: `TriggerTime`): the time reference
 *    for the time values in `anab::T0` side of `TrackTimes`; this is also the
 *    time assumed when a track is not associated to any time.
 *  * `NewTimeReference` (keyword; default: as in `OldTimeReference`): the time
 *    reference for the time values in `anab::T0` side of `NewTimes`.
 *  * `NewTrackAssnToOldT0` (string, default: `"unshifted"`): instance name of
 *    the produced associations between the old tracks and the new times.
 *  * `UpdateAssociations` (configuration table) encloses the options to
 *    recreate the associations of the original, unshifted tracks: they are
 *    specified by type, each as a list of elements, and each element being a
 *    configuration table containing two optional elements:
 *      * `From` (input tag): the input tag to read the original association
 *        from; if not specified, the tag of the input track (`Tracks`) is used.
 *      * `As` (string): create the new association with this instance name
 *        _instead of the one of the original association_; to create an
 *        association with no instance name from an original association that
 *        had one, `As` must be specified empty (`As: ""`).
 *    
 *    Remember that this module can create only one association of a given type
 *    for each instance name.
 *    For example:
 *        
 *        UpdateAssociations: {
 *          PFParticle: [ {}, { From: "Fitter"  As: "kalman" } ],
 *          HitWithMetadata: [ { From: "pandora:meta"  As: "" } ]
 *        }
 *        
 *    will recreate three associations. The first two are to particle flow
 *    objects: an association with the input tag from `Tracks` and write the
 *    modified one with the same instance name as in `Tracks`; and, the
 *    association from the data product called `"Fitter"` (no instance name)
 *    will be processed and written with instance name `"kalman"`.
 *    Then, the association `pandora:meta` with hits is read (see below); the
 *    resulting association will have no instance name (`As: ""`): if this had
 *    not been specified, the output association would have had instance name
 *    `"meta"`.
 *    Thus, to specify an association that was produced by the same module as
 *    the input tracks and with no instance name, the shortest specification is
 *    an empty table.
 *    The following parameters control the supported associations; they are all
 *    lists of specifications as described above, and by default none of them
 *    is enabled.
 *      * `HitWithMetadata`: associations between the new tracks and the hits
 *        with hit metadata.
 *      * `PFParticle`: association between the new tracks and the particle flow
 *        objects.
 * 
 */
class DriftTracks: public art::SharedProducer {
  
    public:
  
  /// Available time references.
  enum class TimeScales { TriggerTime, BeamGateTime };
  
  /// Specification for an association to be replicated.
  struct AssnSpec_t {
    art::InputTag from; /// Tag of the input association.
    std::string instanceName; ///< New instance name.
  };
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    /// Configuration of all associations to be replicated.
    struct Associations {
      
      // Configuration of one type of associations to be replicated.
      struct Spec {
        
        fhicl::OptionalAtom<art::InputTag> From {
          Name{ "From" },
          Comment{
            "tag of the association with the unshifted tracks"
            " [default: same as those tracks]"
            }
          };
        
        fhicl::OptionalAtom<std::string> As {
          Name{ "As" },
          Comment{
            "instance name for the new association"
            " [default: as the original one]"
            }
          };
        
      }; // struct Spec
      
      using AssociationSequence = fhicl::OptionalSequence<fhicl::Table<Spec>>;
      
      AssociationSequence HitWithMetadata {
        Name{ "HitWithMetadata" },
        Comment
          { "associations to recob::Hit with recob::TrackHitMeta metadata" }
        };
      
      AssociationSequence PFParticle {
        Name{ "PFParticle" },
        Comment{ "associations to recob::PFParticle" }
        };
      
    }; // struct Associations
    
    
    /// Selector for `OldTimeReference` and `NewTimeReference` parameters.
    static util::MultipleChoiceSelection<TimeScales> const TimeScaleSelector;
    
    
    fhicl::Atom<art::InputTag> Tracks {
      Name{ "Tracks" },
      Comment{ "input track data product" }
      };
    
    fhicl::Atom<art::InputTag> NewTimes {
      Name{ "NewTimes" },
      Comment{
        "input association of tracks and times to be used to shift the tracks"
        }
      };
    
    fhicl::OptionalAtom<art::InputTag> TrackTimes {
      Name{ "TrackTimes" },
      Comment
        { "input association of tracks and times that were used with them" }
      };
    
    fhicl::Atom<std::string> OldTimeReference {
      Name{ "OldTimeReference" },
      Comment{
        "time reference for the old T0 times: "
        + TimeScaleSelector.optionListString()
        },
      TimeScaleSelector.get(TimeScales::TriggerTime).name()
      };
    
    fhicl::Atom<std::string> NewTimeReference {
      Name{ "NewTimeReference" },
      Comment{
        "time reference for the new T0 times: "
        + TimeScaleSelector.optionListString()
        + " (default: as OldTimeReference)" },
      ""
      };
    
    fhicl::Atom<std::string> NewTrackAssnToOldT0 {
      Name{ "NewTrackAssnToOldT0" },
      Comment{
        "instance name of the produced associations between the old tracks and"
        " the new times"
        },
        "unshifted"
      };
    
    fhicl::OptionalTable<Associations> UpdateAssociations {
      Name{ "UpdateAssociations" },
      Comment
        { "Associations of unshifted tracks to replicate with the new tracks" }
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of message facility category for console message stream" },
      "DriftTracks"
      };
    
    
    template <typename DataMemberPtr>
    TimeScales getTimeScale(DataMemberPtr dataPtr) const
      {
        try {
          return TimeScaleSelector.parse((this->*dataPtr)());
        }
        catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e) {
          throw art::Exception(art::errors::Configuration)
            << "Invalid value for '" << (this->*dataPtr).name()
            << "' parameter: '" << e.label() << "'; valid options: "
            << TimeScaleSelector.optionListString() << ".\n";
        }
      } // getTimeScale()
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  
  DriftTracks(Parameters const& params, const art::ProcessingFrame&);
  
  
  virtual void produce(art::Event& event, const art::ProcessingFrame&) override;
  
  
    private:
  
  using nanoseconds = util::quantities::intervals::nanoseconds;
  
  /// Specifications about the track associations to update.
  template <typename Right, typename Meta = void>
  struct UpdateAssociationSpecs_t: std::vector<AssnSpec_t> {
    using left_t = recob::Track;
    using right_t = Right;
    using data_t = Meta;
    using Assns_t = art::Assns<left_t, right_t, data_t>;
    static constexpr bool HasMetadata = !std::is_void_v<data_t>;
    template <typename... Args> UpdateAssociationSpecs_t(Args&&... args)
      : std::vector<AssnSpec_t>{ std::forward<Args>(args)... } {}
  }; // UpdateAssociationSpecs_t
  
  
  // --- BEGIN --  Configuration  ----------------------------------------------
  
  art::InputTag const fTrackTag; ///< The tracks for input.
  art::InputTag const fOldTimeTag; ///< Track-old T0 associations in input.
  art::InputTag const fNewTimeTag; ///< Track-new T0 associations in input.
  
  /// Input tags of the associations to hits with metadata to replicate.
  std::vector<art::InputTag> const fAssnsToHits;
  /// Input tags of the associations to particele flow objects to replicate.
  std::vector<art::InputTag> const fAssnsToPFParticles;
  
  TimeScales const fOldTimeReference; ///< Time reference for old track times.
  TimeScales const fNewTimeReference; ///< Time reference for new track times.
  
  /// Instance name for the association of new tracks and old times.
  std::string const fNewTrackAssnToOldT0instanceName;
  
  /// Specifications of the associations to recreate.
  /// @{
  UpdateAssociationSpecs_t<recob::Hit, recob::TrackHitMeta> const
    fAssnHitWithMeta;
  UpdateAssociationSpecs_t<recob::PFParticle> const fAssnPFParticle;
  /// @}
  
  std::string const fLogCategory; ///< Category for messagefacility stream.
  
  // --- END ----  Configuration  ----------------------------------------------
  
  
  // --- BEGIN --  Cached services  --------------------------------------------
  
  geo::GeometryCore const& fGeom; ///< Geometry service provider.
  
  // --- END ----  Cached services  --------------------------------------------
  
  
  /// Returns the reference time for the specified time scale.
  detinfo::timescales::electronics_time getTimeReference
    (TimeScales timeRef, detinfo::DetectorTimings const& detTimings) const;
  
  
  // --- BEGIN -- Management of association replica ----------------------------
  /// @name Management of association replica
  /// @{
  
  /// Helper to convert the association specifications from the configuration.
  std::vector<AssnSpec_t> parseAssociationSpecs
    (std::optional<std::vector<Config::Associations::Spec>> const& config) const;
  
  /// Declares the consummation of the needed association data products.
  template <typename AssociationSpecs>
  void consumesReplicaAssociations(AssociationSpecs const& specs);
  
  /// Declares the production of the needed association data products.
  template <typename AssociationSpecs>
  void producesReplicaAssociations(AssociationSpecs const& specs);
  
  /// Prints the configuration for the replica of the specified associations.
  template <typename Stream, typename AssociationSpecs>
  void printAssociationConfig(
    Stream&& out, AssociationSpecs const& specs,
    std::string const& indent = "",
    std::string const& rightName = "", std::string const& metaName = ""
    ) const;
  
  /// Replicates and stores all the associations in the specifications `specs`.
  template <typename AssociationSpecs, typename ReplicaMap>
  void replicateAndPutAssociations(
    art::Event& event, AssociationSpecs const& specs, ReplicaMap const& firstMap
    ) const;
  
  /// Helper to supply a "empty default" to the association specifications.
  static Config::Associations AssociationParams(Parameters const& params);
  
  template <typename T>
  static std::string nameOf(std::string const& name = "");
  
  /// @}
  // --- END ---- Management of association replica ----------------------------
  
}; // class DriftTracks


//------------------------------------------------------------------------------
//---  implementation
//------------------------------------------------------------------------------
namespace {
  
  template <typename T>
  std::unique_ptr<T> move_to_unique_ptr(T&& data)
    { return std::make_unique<std::remove_reference_t<T>>(std::move(data)); }
  
} // local namespace


//------------------------------------------------------------------------------
//---  DriftTracks
//------------------------------------------------------------------------------
util::MultipleChoiceSelection<DriftTracks::TimeScales> const
DriftTracks::Config::TimeScaleSelector
  {
      { TimeScales::TriggerTime,  "TriggerTime",  "Trigger" }
    , { TimeScales::BeamGateTime, "BeamGateTime", "BeamTime", "BeamGate" }
  };


//------------------------------------------------------------------------------
DriftTracks::DriftTracks(Parameters const& params, const art::ProcessingFrame&)
  : art::SharedProducer{ params }
  // configuration parameters
  , fTrackTag{ params().Tracks() }
  , fOldTimeTag{ params().TrackTimes().value_or(fTrackTag) }
  , fNewTimeTag{ params().NewTimes() }
  , fOldTimeReference{ params().getTimeScale(&Config::OldTimeReference) }
  , fNewTimeReference{
      params().NewTimeReference().empty()
        ? fOldTimeReference: params().getTimeScale(&Config::NewTimeReference)
      }
  , fNewTrackAssnToOldT0instanceName{ params().NewTrackAssnToOldT0() }
  , fAssnHitWithMeta
    { parseAssociationSpecs(AssociationParams(params).HitWithMetadata()) }
  , fAssnPFParticle
    { parseAssociationSpecs(AssociationParams(params).PFParticle()) }
  , fLogCategory{ params().LogCategory() }
  // caches
  , fGeom{ *(lar::providerFrom<geo::Geometry const>()) }
{
  async<art::InEvent>();
  
  //
  // declaration of input
  //
  consumes<std::vector<recob::Track>>(fTrackTag);
  if (!fOldTimeTag.empty())
    consumes<art::Assns<recob::Track, anab::T0>>(fOldTimeTag);
  consumes<art::Assns<recob::Track, anab::T0>>(fNewTimeTag);
  consumesReplicaAssociations(fAssnHitWithMeta);
  consumesReplicaAssociations(fAssnPFParticle);
  
  //
  // declaration of output
  //
  produces<std::vector<recob::Track>>();
  produces<art::Assns<recob::Track, anab::T0>>();
  produces<art::Assns<recob::Track, anab::T0>>
    (fNewTrackAssnToOldT0instanceName);
  producesReplicaAssociations(fAssnHitWithMeta);
  producesReplicaAssociations(fAssnPFParticle);
  
  //
  // configuration dump
  //
  {
    mf::LogInfo log { fLogCategory };
    log << "Shifting tracks from '" << fTrackTag.encode()
      << "' using the times from '" << fNewTimeTag.encode() << "'"
      << " (reference: " << Config::TimeScaleSelector.get(fNewTimeReference)
      << ")";
    if (fOldTimeTag.empty())
      log << "\n - times used for input track: reference time for all tracks";
    else {
      log << "\n - times used for input track: '" << fOldTimeTag.encode()
        << "'";
    }
    log
      << " (reference: " << Config::TimeScaleSelector.get(fOldTimeReference)
        << ")"
      << "\nReplicated associations (if any):"
      ;
    printAssociationConfig(log, fAssnHitWithMeta, "  ");
    printAssociationConfig(log, fAssnPFParticle, "  ");
  }
  
} // DriftTracks::DriftTracks()


//------------------------------------------------------------------------------
void DriftTracks::produce(art::Event& event, art::ProcessingFrame const&) {
  
  using namespace util::quantities::time_literals;
  
  auto const detTimings = detinfo::makeDetectorTimings(
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    );
  icarus::TrackShifterAlg const shiftTrack {
      fLogCategory
    , fGeom
    , art::ServiceHandle<detinfo::DetectorPropertiesService const>()
      ->DataFor(event, detTimings.clockData())  // detProp
    };
  
  //
  // read inputs
  //
  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector
    (tracks, event.getHandle<std::vector<recob::Track>>(fTrackTag));
  
  std::optional<art::FindOneP<anab::T0>> trackToOldTime;
  if (!fOldTimeTag.empty())
    trackToOldTime.emplace(tracks, event, fOldTimeTag);
  art::FindOneP<anab::T0> trackToNewTime{ tracks, event, fNewTimeTag };
  
  //
  // establish reference times and time correction
  //
  auto const oldReference = getTimeReference(fOldTimeReference, detTimings);
  auto const newReference = getTimeReference(fNewTimeReference, detTimings);
  // add to old to get new:
  nanoseconds const oldToNewTimeRef = newReference - oldReference;
  
  //
  // prepare the output
  //
  art::PtrMaker<recob::Track> makeTrackPtr{ event };
  std::vector<recob::Track> shiftedTracks;
  std::unordered_map<art::Ptr<recob::Track>, art::Ptr<recob::Track>>
    shiftedTrackPtrs;
  art::Assns<recob::Track, anab::T0> newTimeToTracks;
  art::Assns<recob::Track, anab::T0> oldTimeToTracks;
  
  for (art::Ptr<recob::Track> const& track: tracks) {
    
    //
    // find the pertinent times
    //
    art::Ptr<anab::T0> const oldT0 = trackToOldTime
      ? trackToOldTime->at(track.key()): art::Ptr<anab::T0>{};
    art::Ptr<anab::T0> const newT0 = trackToNewTime.at(track.key());
    
    //
    // compute the time shift
    //
    nanoseconds trackShift { 0_ns }; // no shift by default
    if (newT0) {
      nanoseconds const newTrackTime { newT0->Time() };
      
      nanoseconds oldTrackTime = oldToNewTimeRef;
      if (oldT0) {
        mf::LogTrace{ fLogCategory }
          << "Track ID=" << track->ID() << " had time " << oldT0->Time() << " ns";
        oldTrackTime += nanoseconds{ oldT0->Time() };
      }
      else if (trackToOldTime) {
        // if there is no old time data product, don't bother with this message
        mf::LogTrace{ fLogCategory }
          << "Track ID=" << track->ID() << " had no time (using reference)";
      }
      
      trackShift = newTrackTime - oldTrackTime;
      mf::LogTrace{ fLogCategory }
        << "  new time: " << newTrackTime << ", shift: " << trackShift;
    } // if new T0
    else {
      mf::LogTrace{ fLogCategory }
        << "Track ID=" << track->ID() << " has no new time: won't be shifted";
    }
    
    //
    // shift the track
    //
    recob::Track shiftedTrack = shiftTrack(*track, trackShift);
  
    //
    // update the output collections
    //
    std::size_t const trackIndex = shiftedTracks.size(); // before adding
    
    art::Ptr<anab::T0> const targetT0 = newT0? newT0: oldT0;
    art::Ptr<recob::Track> newTrackPtr = makeTrackPtr(trackIndex);
    
    if (oldT0) oldTimeToTracks.addSingle(newTrackPtr, oldT0);
    if (targetT0) newTimeToTracks.addSingle(newTrackPtr, targetT0);
    
    // add the shifted track
    shiftedTracks.push_back(std::move(shiftedTrack));
    
    shiftedTrackPtrs[track] = std::move(newTrackPtr);
    
  } // for tracks
  
  //
  // save the data products
  //
  event.put(move_to_unique_ptr(std::move(shiftedTracks)));
  event.put(move_to_unique_ptr(std::move(newTimeToTracks)));
  event.put(
    move_to_unique_ptr(std::move(oldTimeToTracks)),
    fNewTrackAssnToOldT0instanceName
    );
  
  //
  // remap the associations
  //
  auto const replacer = [&shiftedTrackPtrs](auto const& track)
    { return shiftedTrackPtrs.at(track); };
  replicateAndPutAssociations(event, fAssnHitWithMeta, replacer);
  replicateAndPutAssociations(event, fAssnPFParticle, replacer);
  
} // DriftTracks::produce()


//------------------------------------------------------------------------------
detinfo::timescales::electronics_time DriftTracks::getTimeReference
  (TimeScales timeRef, detinfo::DetectorTimings const& detTimings) const
{
  switch (timeRef) {
    case TimeScales::TriggerTime:  return detTimings.TriggerTime();
    case TimeScales::BeamGateTime: return detTimings.BeamGateTime();
    default:
      throw art::Exception{ art::errors::LogicError }
        << "Unexpected time reference: '"
        << Config::TimeScaleSelector.get(timeRef) << "'\n"
        ;
  } // switch
  
} // DriftTracks::getTimeReference()


//------------------------------------------------------------------------------
auto DriftTracks::parseAssociationSpecs
  (std::optional<std::vector<Config::Associations::Spec>> const& config) const
  -> std::vector<AssnSpec_t>
{
  if (!config) return {};
  
  std::vector<AssnSpec_t> specs;
  
  for (Config::Associations::Spec const& oneConfig: *config) {
    art::InputTag from = oneConfig.From().value_or(fTrackTag);
    std::string instanceName = oneConfig.As().value_or(from.instance());
    specs.push_back({ std::move(from), std::move(instanceName) });
  } // for
  
  return specs;
} // DriftTracks::parseAssociationSpecs()


//------------------------------------------------------------------------------
auto DriftTracks::AssociationParams(Parameters const& params)
  -> Config::Associations
{
  return params().UpdateAssociations()
    .value_or(DriftTracks::Config::Associations{});
}


//------------------------------------------------------------------------------
template <typename AssociationSpecs>
void DriftTracks::consumesReplicaAssociations(AssociationSpecs const& specs) {
  
  for (AssnSpec_t const& spec: specs)
    consumes<typename AssociationSpecs::Assns_t>(spec.from);
  
} // DriftTracks::consumesReplicaAssociations()


//------------------------------------------------------------------------------
template <typename AssociationSpecs>
void DriftTracks::producesReplicaAssociations(AssociationSpecs const& specs) {
  
  for (AssnSpec_t const& spec: specs)
    produces<typename AssociationSpecs::Assns_t>(spec.instanceName);
  
} // DriftTracks::producesReplicaAssociations()


//------------------------------------------------------------------------------
template <typename AssociationSpecs, typename ReplicaMap>
void DriftTracks::replicateAndPutAssociations(
  art::Event& event, AssociationSpecs const& specs, ReplicaMap const& firstMap
) const {
  using Assns_t = typename AssociationSpecs::Assns_t;
  
  for (auto const& [ tag, instance ]: specs) {
    
    auto const& origAssns = event.getProduct<Assns_t>(tag);
    
    event.put(
      move_to_unique_ptr(util::replaceAssnsFirst(origAssns, firstMap)),
      instance
      );
    
  } // for associations
  
} // DriftTracks::replicateAndPutAssociations()


//------------------------------------------------------------------------------
template <typename Stream, typename AssociationSpecs>
void DriftTracks::printAssociationConfig(
  Stream&& out, AssociationSpecs const& specs,
  std::string const& indent /* = "", */,
  std::string const& rightName /* = "" */,
  std::string const& metaName /* = "" */
) const {
  
  using Right_t = typename AssociationSpecs::right_t;
  using Meta_t = typename AssociationSpecs::data_t;
  
  if (specs.empty()) return; // print nothing
  
  out << "\n - " << specs.size() << " associations: recob::Track <=> "
    << nameOf<Right_t>(rightName);
  if constexpr(AssociationSpecs::HasMetadata) {
    out << " (metadata: " << nameOf<Meta_t>(metaName) << ")";
  }
  out << ":";
  for (AssnSpec_t const& spec: specs) {
    out << indent << "\n   * replica of '" << spec.from.encode() << "'";
    if (!spec.instanceName.empty())
      out << " (instance name: '" << spec.instanceName << "')";
  } // for
  
} // DriftTracks::printAssociationConfig()


//------------------------------------------------------------------------------
template <typename T>
std::string DriftTracks::nameOf(std::string const& name /* = "" */)
  { return name.empty()? lar::debug::demangle<T>(): name; }


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(DriftTracks)


//------------------------------------------------------------------------------

