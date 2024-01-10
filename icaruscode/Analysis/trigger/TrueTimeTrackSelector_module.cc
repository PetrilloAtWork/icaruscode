/**
 * @file    icaruscode/Analysis/trigger/TrueTimeTrackSelector_module.cc
 * @date    September 17, 2021
 * @authors Animesh Chatterjee (ANC238@pitt.edu),
 *          Gianluca Petrillo (petrillo@slac.stanford.edu),
 *          Jacob Zettlemoyer (jzettle@fnal.gov)
 * 
 */

// #define MF_DEBUG // to always allow MF_LOG_DEBUG()/MF_LOG_TRACE() output
#define USE_ATOMICPASSCOUNTER 1

// SBN and ICARUS libraries
#ifdef USE_ATOMICPASSCOUNTER
#include "icarusalg/Utilities/AtomicPassCounter.h"
#else // !USE_ATOMICPASSCOUNTER
#include "icarusalg/Utilities/PassCounter.h"
#endif // USE_ATOMICPASSCOUNTER

// LArSoft libraries
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventory.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"

// framework libraries
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard library
#include <algorithm> // std::sort(), std::unique()
#include <vector>
#include <string>
#include <atomic>
#include <memory>
#include <limits>

// -----------------------------------------------------------------------------
namespace sbn { class TrueTimeTrackSelector; }
/**
 * @brief Selects all simulated tracks, adding true time information.
 * 
 * This module produces a list of "tracks" that are associated to a time.
 * Optionally, it can select only the tracks with that time in a specified
 * range.
 * 
 * This module is a filter that will return `true` for an event if at least one
 * of its tracks is selected. This threshold can be adjusted by configuration.
 * 
 * 
 * Input
 * ------
 * 
 * The input file is expected to contain the track collection(s) to be processed
 * and enough truth information to allow backtracking.
 * 
 * 
 * Output
 * -------
 * 
 * The filter _produces_ a list of pointers to the selected tracks; the list
 * will mix tracks from different data products if multiple input collections
 * are specified. A track is selected if all the following apply:
 * 
 * 1. the track can be associated to a time;
 * 2. only if a time range is specified, that track time must fall within.
 * 
 * The filter passes the event if:
 *  * the number of selected tracks is within the configured range of requested
 *    tracks per event.
 * 
 * The output includes:
 *  * `std::vector<art::Ptr<recob::Track>>` (if `SaveTracks` is set): a list of
 *    the pointers to the selected tracks.
 *  * `std::vector<anab::T0>`: a list of times for the input tracks (in the same
 *    order as the input tracks are processed). The object contains:
 *      * `Time()`: the backtracked time of the track, in Pandora's T0
 *        convention of @ref DetectorClocksTriggerTime "trigger time scale"
 *        (for simulation it is actually beam gate opening scale) but in
 *        nanoseconds, or a large negative number if time is not valid;
 *      * `TriggerType()`: `2` ("Monte Carlo truth" according to documentation);
 *      * `TriggerBits()`: `1` if time is valid, `0` otherwise;
 *      * `ID()`: same ID as the track being tagged;
 *      * `TriggerConfidence()`: `1.0` if the time is valid, `-1.0` (default)
 *        otherwise.
 *    Therefore it is possible to identify the tracks which couldn't be
 *    associated to a time by testing either `TriggerBits()` or
 *    `TriggerConfidence()` value.
 *  * `art::Assns<recob::Track, anab::T0>`: the list of selected tracks and
 *    the `anab::T0` they are associated to. Only the track with a valid time
 *    are included in this association.
 * 
 * Configuration options
 * ----------------------
 * 
 * * `TrackTags` (list of data product tags, required): data product of the
 *    tracks to be selected.
 * * `MinT0` (real, optional): if specified, tracks are selected only if their
 *     associated time is not earlier than this value. Time is in the same time
 *     scale as the associated track time, which is expected to be the
 *     @ref DetectorClocksTriggerTime "trigger time scale" but in nanoseconds.
 * * `MaxT0` (real, optional): if specified, tracks are selected only if their
 *     associated time is earlier than this value. This time is in the same
 *     time scale as `MinT0`.
 * * `MinTracks` (integer, default: `1`): the filter "passes" the event only if
 *     at least these many tracks are selected; disable this by setting it to
 *     `0`.
 * * `MaxTracks` (integer, default: a large number): if specified, filter
 *     "passes" the event only if at most these many tracks are selected.
 * * `SaveTracks` (flag, default: `true`): produces a list of _art_ pointers
 *     to the selected particles. This is the default behaviour of the module.
 *     On the other end, the module can be used solely for its filtering
 *     capability, in which case saving the list is not necessary and this
 *     option allows omitting it.
 * * `LogCategory` (string, default: `TrueTimeTrackSelector`): name of the
 *     message facility stream where to send the module consoed output to.
 * 
 */
class sbn::TrueTimeTrackSelector: public art::SharedFilter {
  
public:
  
  static constexpr double NoMinTime { std::numeric_limits<double>::lowest() };
  static constexpr double NoMaxTime { std::numeric_limits<double>::max() };
  static constexpr double NoAngleCut {std::numeric_limits<double>::lowest() };
  static constexpr unsigned int NoMinTracks { 0U };
  static constexpr unsigned int NoMaxTracks
    { std::numeric_limits<unsigned int>::max() };
  
  /// Information about one input track data product.
  struct TrackTags_t {
    art::InputTag tracks;
    art::InputTag trackToHits;
  };

  /// Module configuration parameters.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    struct TrackInputConfig {
      
      fhicl::Atom<art::InputTag> TrackTag {
        Name{ "TrackTag" },
        Comment{ "data product with tracks to be selected" }
      };
      
      fhicl::OptionalAtom<art::InputTag> AssociatedHitsTag {
        Name{ "AssociatedHitsTag" },
        Comment{
          "associations between input tracks and their hits (default: same tag as tracks)"
        }
      };
      
    }; // TrackInputConfig
    
    
    fhicl::Sequence<fhicl::TableAs<TrackTags_t, TrackInputConfig>> InputTracks {
      Name{ "InputTracks" },
      Comment{ "tags of tracks to be selected and additional information" }
    };
    
    fhicl::Atom<double> MinT0 {
      Name{ "MinT0" },
      Comment{ "Select only tracks not earlier than this time [ns]" },
      NoMinTime
    };
    
    fhicl::Atom<double> MaxT0 {
      Name{ "MaxT0" },
      Comment{ "Select only tracks earlier than this time [ns]" },
      NoMaxTime
    };
    
    fhicl::Atom<unsigned int> MinTracks {
      Name{ "MinTracks" },
      Comment{ "Pass only events with at least these many selected tracks" },
      1U
    };
    
    fhicl::Atom<unsigned int> MaxTracks {
      Name{ "MaxTracks" },
      Comment{ "Pass only events with at most these many selected tracks" },
      NoMaxTracks
    };
    
    fhicl::Atom<bool> SaveTracks {
      Name{ "SaveTracks" },
      Comment{ "Whether to write the list of selected tracks" },
      true
    };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of a message facility stream for this module" },
      "TrueTimeTrackSelector"
    };
    
  }; // Config
  
  
  using Parameters = art::SharedFilter::Table<Config>;
  
  explicit TrueTimeTrackSelector
    (Parameters const& params, art::ProcessingFrame const&);
  
  bool filter(art::Event& event, art::ProcessingFrame const&) override;
  
  /// Prints end-job summary.
  void endJob(art::ProcessingFrame const&) override;
  

private:
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  std::vector<TrackTags_t> const fTrackTags; ///< List of track input info.
  
  double const fMinT0;  ///< Minimum track time for track selection.
  double const fMaxT0;  ///< Maximum track time for track selection.

  /// Minimum selected tracks for event selection.
  unsigned int const fMinTracks;
  /// Maximum selected tracks for event selection.
  unsigned int const fMaxTracks;
  
  bool const fSaveTracks;  ///< Whether to save selected tracks into the event.
  
  std::string const fLogCategory; ///< Message facility stream name.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  // --- BEGIN -- Caches -------------------------------------------------------
  
  detinfo::DetectorClocksService const& fDetClocks;
  cheat::BackTracker const& fBackTracker;
  cheat::ParticleInventory const& fParticleInv;
  
  // --- END ---- Caches -------------------------------------------------------
  
  /// Counter of passed events (not thread-safe).
#ifdef USE_ATOMICPASSCOUNTER
  icarus::ns::util::AtomicPassCounter<> fPassRate;
#else // !USE_ATOMICPASSCOUNTER
  icarus::ns::util::PassCounter<> fPassRate;
#endif // USE_ATOMICPASSCOUNTER
  
  
  /**
   * @brief Adds to `selectedTracks` qualifying tracks from `tracks`.
   * @param trackHandle _art_ handle to the tracks to be selected
   * @param makeT0ptr pointer maker for the module T0 data product
   * @param[out] times collection to expand with track times (one per input)
   * @param[out] selectedTracksAndT0 collection to expand with qualifying tracks
   * @return the number of qualifying tracks found in `tracks` and added
   */
  template <typename TrackHandle>
  unsigned int selectTracks(
    TrackHandle const& trackHandle,
    art::FindMany<recob::Hit> const& trackHitAssns,
    art::PtrMaker<anab::T0> const& makeT0ptr,
    detinfo::DetectorClocksData const& clockData,
    std::vector<anab::T0>& times,
    art::Assns<recob::Track, anab::T0>& selectedTracksAndT0
    ) const;
  
  /**
   * @brief Returns the `anab::T0` to associate to the specified `track`.
   * @param track the track to be associated to a time
   * @param trackHits sequence of pointers to the hits associated to `track`
   * @param clockData detector clocks data (needed for backtracking)
   * @return the time object, with special values if no time could be associated
   * 
   * In case the track can't be associated to a time, the returned object
   * follows the prescription described in the module class documentation.
   */
  anab::T0 findTrackTime(
    recob::Track const& track,
    std::vector<recob::Hit const*> const& trackHits,
    detinfo::DetectorClocksData const& clockData
    ) const;
  
  /// Returns whether the specified track (with specified time) qualifies.
  bool isTrackSelected
    (recob::Track const& track, anab::T0 const& time) const;

  /// Returns if the number of tracks `nTracks` satisfies filter requirements.
  bool selectedTracksRequirement(unsigned int nTracks) const;
  
  /// Unpacks input track data product configuration.
  friend TrackTags_t convert(Config::TrackInputConfig const& config);
  
}; // sbn::TrueTimeTrackSelector

// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  /// Moves the content of `obj` into a new `std::unique_ptr` and returns it.
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& obj)
    { return std::make_unique<T>(std::move(obj)); }
  
  // ---------------------------------------------------------------------------
  /// Appends copies of all the elements of `src` to the end of `dest`.
  /// @return `dest`
  template <typename DestColl, typename SrcColl>
  DestColl& append(DestColl& dest, SrcColl const& src) {
    using std::begin, std::end;
    dest.insert(end(dest), begin(src), end(src));
    return dest;
  } // append()

  // ---------------------------------------------------------------------------
  /// Sorts and removes duplicate elements of `coll` in place.
  template <typename Coll>
  void makeUniqueCollection(Coll& coll) {
    using std::begin, std::end;
    std::sort(begin(coll), end(coll));
    auto const newEnd = std::unique(begin(coll), end(coll));
    coll.erase(newEnd, end(coll));
  } // makeUniqueCollection()

  // ---------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
namespace sbn {

  sbn::TrueTimeTrackSelector::TrackTags_t convert
    (sbn::TrueTimeTrackSelector::Config::TrackInputConfig const& config)
  {
    art::InputTag const trackTag = config.TrackTag();
    return {
        trackTag                                       // tracks
      , config.AssociatedHitsTag().value_or(trackTag)  // trackToHits
    };
  } // convert(sbn::TrueTimeTrackSelector::Config::TrackInputConfig)
  
} // namespace sbn

// -----------------------------------------------------------------------------
// ---  sbn::TrueTimeTrackSelector
// -----------------------------------------------------------------------------
sbn::TrueTimeTrackSelector::TrueTimeTrackSelector
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedFilter{ params }
    // configuration
  , fTrackTags{ params().InputTracks() }
  , fMinT0{ params().MinT0() }
  , fMaxT0{ params().MaxT0() }
  , fMinTracks{ params().MinTracks() }
  , fMaxTracks{ params().MaxTracks() }
  , fSaveTracks{ params().SaveTracks() }
  , fLogCategory{ params().LogCategory() }
    // caches
  , fDetClocks{ *art::ServiceHandle<detinfo::DetectorClocksService const>() }
  , fBackTracker{ *lar::providerFrom<cheat::BackTrackerService const>() }
  , fParticleInv{ *lar::providerFrom<cheat::ParticleInventoryService const>() }
{
  async<art::InEvent>();
  
  //
  // consume declarations
  //
  for (auto const& [ trackTag, trackHitAssnTag ]: fTrackTags) {
    consumes<std::vector<recob::Track>>(trackTag);
    consumes<art::Assns<recob::Track, recob::Hit>>(trackTag);
  } // for
  
  //
  // produce declarations
  //
  if (fSaveTracks) {
    produces<std::vector<art::Ptr<recob::Track>>>();
  }
  produces<std::vector<anab::T0>>();
  produces<art::Assns<recob::Track, anab::T0>>();
  
  //
  // configuration dump
  //
  {
    mf::LogInfo log{ fLogCategory };
    log << "Configuration:"
      << "\n  * tracks required to be back-tracked to a time"
      ;
    log << "\n  * " << fTrackTags.size() << " input track collections:";
    for (auto const& [ trackTag, trackHitAssnTag ]: fTrackTags) {
      log << " " << trackTag.encode();
      if (trackTag != trackHitAssnTag)
        log << " (+" << trackHitAssnTag.encode() << ")";
    } // for
    log << "\n  * track time:";
    if (fMinT0 == NoMinTime) {
      if (fMaxT0 == NoMaxTime) log << " any";
      else                     log << " before " << fMaxT0;
    }
    else {
      log << " from " << fMinT0;
      if (fMaxT0 == NoMaxTime) log << " on";
      else                     log << " to " << fMaxT0;
    }
    
    log << "\n  * selected tracks per event: ";
    if (fMinTracks == NoMinTracks) {
      if (fMaxTracks == NoMaxTracks) log << "any";
      else log << fMaxTracks << " or less";
    }
    else {
      if (fMaxTracks == NoMaxTracks) log << fMinTracks << " or more";
      else log << "between " << fMinTracks << " and " << fMaxTracks;
    }
    log << "\n  * selected tracks will" << (fSaveTracks? "": " not")
      << " be saved";
  } // end local scope
  
} // sbn::TrueTimeTrackSelector::TrueTimeTrackSelector()


// -----------------------------------------------------------------------------
bool sbn::TrueTimeTrackSelector::filter
(art::Event& event, art::ProcessingFrame const&)
{
  
  mf::LogDebug(fLogCategory) << "Processing " << event.id();
  
  //
  // select tracks from each of the input tags
  //
  unsigned int nTracks = 0;
  art::PtrMaker<anab::T0> makeT0ptr{ event };
  detinfo::DetectorClocksData const clockData = fDetClocks.DataFor(event);
  
  std::vector<anab::T0> Times;
  art::Assns<recob::Track, anab::T0> selectedTracksAndT0;
  for (auto const& [ trackTag, trackHitAssnTag ]: fTrackTags) {
    auto const& trackHandle
      = event.getValidHandle<std::vector<recob::Track>>(trackTag);
    
    art::FindMany<recob::Hit> tracksToHits
      { trackHandle, event, trackHitAssnTag};
    
    nTracks += trackHandle->size();
    unsigned int const newTracks = selectTracks(
      trackHandle, tracksToHits,
      makeT0ptr, clockData,
      Times, selectedTracksAndT0
      );
    
    mf::LogTrace(fLogCategory)
      << "From '" << trackTag.encode() << "': "
      << newTracks << " tracks selected"
      ;
    
  } // for

  std::vector<art::Ptr<recob::Track>> selectedTracks;
  for (auto const& TrackAndT0: selectedTracksAndT0)
    selectedTracks.push_back(TrackAndT0.first);
  unsigned int const nSelectedTracks = selectedTracks.size();
  
  
  //
  // filter logic
  //
  
  bool const passed = selectedTracksRequirement(nSelectedTracks);
  fPassRate.add(passed);
  mf::LogTrace(fLogCategory) << event.id()
    << ' ' << (passed? "passed": "rejected") << " (" << nSelectedTracks
    << "/" << nTracks << " selected tracks).";
  
  //
  // save track list in the event
  //
  
  // after this, selectedTracks may be empty.
  if (fSaveTracks) {
    if (!selectedTracks.empty()) {
      mf::LogTrace(fLogCategory)
        << "InputTag for this product is: "
        << event.getProductDescription(selectedTracks[0].id())->inputTag();
    }
    event.put(moveToUniquePtr(selectedTracks));
  }
  event.put(moveToUniquePtr(Times));
  event.put(moveToUniquePtr(selectedTracksAndT0));
  
  return passed;
  
} // sbn::TrueTimeTrackSelector::filter()


// -----------------------------------------------------------------------------
void sbn::TrueTimeTrackSelector::endJob(art::ProcessingFrame const&) {
  
  mf::LogInfo(fLogCategory) << "Selected " << fPassRate.passed()
    << '/' << fPassRate.total() << " events with qualifying tracks.";
  
} // sbn::TrueTimeTrackSelector::endJob()


// -----------------------------------------------------------------------------
template <typename TrackHandle>
unsigned int sbn::TrueTimeTrackSelector::selectTracks(
  TrackHandle const& trackHandle,
  art::FindMany<recob::Hit> const& tracksToHits,
  art::PtrMaker<anab::T0> const& makeT0ptr,
  detinfo::DetectorClocksData const& clockData,
  std::vector<anab::T0>& times,
  art::Assns<recob::Track, anab::T0>& selectedTracksAndT0
) const {
  static_assert
    (std::is_same_v<typename TrackHandle::element_type::value_type, recob::Track>);
  
  unsigned int nSelectedTracks { 0U };
  for (auto const& [ iTrack, track ]: util::enumerate(*trackHandle)) {
    
    std::vector<recob::Hit const*> const& trackHits = tracksToHits.at(iTrack);
    
    anab::T0 const t0 = findTrackTime(track, trackHits, clockData);
    times.push_back(t0); // unconditional
    
    if (!isTrackSelected(track, t0)) continue;
    
    art::Ptr<recob::Track> trackPtr{ trackHandle, iTrack };
    
    MF_LOG_TRACE(fLogCategory) << "Track #" << trackPtr.key() << " selected.";
    
    selectedTracksAndT0.addSingle(trackPtr, makeT0ptr(times.size() - 1));
    ++nSelectedTracks;
  } // for
  
  return nSelectedTracks;
  
} // sbn::TrueTimeTrackSelector::selectTracks()


// -----------------------------------------------------------------------------

anab::T0 sbn::TrueTimeTrackSelector::findTrackTime(
  recob::Track const& track,
  std::vector<recob::Hit const*> const& trackHits,
  detinfo::DetectorClocksData const& clockData
) const {
  
  /*
   * 1. finds all particles (via their G4 ID) contributing to a track
   * 2. finds all particle ancestors of those particles
   * 3. success if there is only one ancestor
   * 
   * The ancestor selection may be refined by returning the one responsible of
   * the most energy in the hits (should go via `TrackIDE` in that case).
   */
  
  std::vector<int> particleIDs;
  for (recob::Hit const* hit: trackHits) {
    if (!hit) continue; // not sure if it's possible
    append(particleIDs, fBackTracker.HitToTrackIds(clockData, *hit));
  } // for hits
  
  makeUniqueCollection(particleIDs);
  
  std::vector<int> eveIDs;
  for (int ID: particleIDs)
    eveIDs.push_back(fParticleInv.TrackIdToEveTrackId(ID));
  
  makeUniqueCollection(eveIDs);
  
  simb::MCParticle const* particle = (eveIDs.size() == 1)
    ? fParticleInv.TrackIdToParticle_P(eveIDs.front()): nullptr;
  if (!particle) {
    
    if (eveIDs.empty()) { // backtrack failed
      mf::LogDebug{ fLogCategory } << "Track ID=" << track.ID()
        << " could not be tracked back to any true particle.";
    }
    else if (eveIDs.size() == 1) { // backtracking succeeded but...
      mf::LogDebug{ fLogCategory }
        << "Track ID=" << track.ID()
        << " was tracked back to particle ID=" << eveIDs.front()
        << ", which is not available: skipped.";
    }
    else { // backtracking pulled in too many particles
      mf::LogDebug log{ fLogCategory };
      auto it = eveIDs.cbegin(), eend = eveIDs.cend();
      log << "Track ID=" << track.ID()
        << " was tracked back to " << eveIDs.size() << " particles: " << *it;
      while (++it != eend) log << ", " << *it;
      log << "; skipped.";
    }
    
    return anab::T0{
        std::numeric_limits<double>::lowest()  // Time
      , 2                                      // TriggerType
      , 0                                      // TriggerBits
      , track.ID()                             // ID
      , -1.0                                   // TriggerConfidence
      };
  } // if backtracking failed in any way
  
  // simb::MCParticle stores times in nanoseconds from generation time,
  // which is assumed to be the same as beam gate time
  // (also, `detinfo::DetectorTimings` can't correctly convert between the two)
  double const time = particle->T();
  
  return anab::T0{
      time        // Time
    , 2           // TriggerType
    , 1           // TriggerBits
    , track.ID()  // ID
    , +1.0        // TriggerConfidence
    };
  
} // sbn::TrueTimeTrackSelector::findTrackTime()


// -----------------------------------------------------------------------------
bool sbn::TrueTimeTrackSelector::isTrackSelected
(recob::Track const& track, anab::T0 const& time) const
{
  
  double const T0 = time.Time();
  MF_LOG_TRACE(fLogCategory) << "Track time: " << T0;

  if ((T0 < fMinT0) || (T0 >= fMaxT0)) {
    MF_LOG_TRACE(fLogCategory) << "Time out of range [ " << fMinT0 << "; "
      << fMaxT0 << " ] => discarded!";
    return false;
  }

  return true;
} // sbn::TrueTimeTrackSelector::isTrackSelected()


// -----------------------------------------------------------------------------
bool sbn::TrueTimeTrackSelector::selectedTracksRequirement
(unsigned int nTracks) const
{
  return (nTracks >= fMinTracks) && (nTracks <= fMaxTracks);
} // sbn::TrueTimeTrackSelector::selectedTracks()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::TrueTimeTrackSelector)


// -----------------------------------------------------------------------------
