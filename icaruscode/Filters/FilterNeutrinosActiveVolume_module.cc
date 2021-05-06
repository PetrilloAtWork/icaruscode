/**
 * @file   icaruscode/Filters/FilterNeutrinosActiveVolume_module.cc
 * @brief  Discards events with no neutrino interaction in the specified volume.
 * @authors Christian Farnese (farnese@pd.infn.it),
 *          Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 6, 2020
 * 
 * Run
 *     
 *     lar --print-description FilterNeutrinosActiveVolume
 *     
 * for configuration directions.
 * 
 */


// ICARUS libraries
#include "icarusalg/Utilities/WeakCurrentType.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/MCDumpers/MCDumperUtils.h" // sim::TruthXXXName(), ...
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ROOTGeometryNavigator.h"
#include "larcorealg/Geometry/GeoNodePath.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::...::makeFromCoords()
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::begin(), util::end()
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT libraries
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

// C++ standard libraries
#include <ostream>
#include <regex>
#include <algorithm> // std::sort(), std::binary_search()
#include <vector>
#include <string>
#include <atomic>
#include <utility> // std::move()
#include <cmath> // std::abs()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::simfilter::details {

  /// A particle, defined by its belonging to one of the specified flavors.
  struct InteractingParticle_t {
    std::vector<int> flavors; ///< List of flavors as PDG ID.
    auto begin() const noexcept { return flavors.begin(); }
    auto end() const noexcept { return flavors.end(); }
  }; // InteractingParticle_t

  std::ostream& operator<< (std::ostream& out, InteractingParticle_t const&);

  /// An interaction between particles.
  struct InteractionComposition_t {
    using Particles_t = std::vector<InteractingParticle_t>;
    /// List of particles involved in the interaction.
    Particles_t particles;
    auto begin() const noexcept { return particles.begin(); }
    auto end() const noexcept { return particles.end(); }
    std::size_t size() const noexcept { return particles.size(); }
    bool empty() const noexcept { return particles.empty(); }
  }; // InteractionComposition_t

  std::ostream& operator<< (std::ostream& out, InteractionComposition_t const&);
  
} // namespace icarus::simfilter::details


// -----------------------------------------------------------------------------
namespace icarus::simfilter { class FilterNeutrinosActiveVolume; }

/**
 * @brief Accepts only neutrino-like events with vertices in a specified volume.
 * 
 * The module identifies a list of volumes to filter on.
 * The filter selects "qualifying" neutrino interactions that:
 * 
 *  * have interaction vertex in one of a selected list of volumes
 *  * optionally, has a specific interaction type
 *    (like in `simb::MCNeutrino::InteractionType()`)
 *  * optionally, has a specific current type (see `simb::MCNeutrino::CCNC()`)
 * 
 * The event is kept if there is _at least_ one qualifying interaction.
 * In that case, the whole event is passed (including any other interaction).
 * 
 * 
 * Input
 * ======
 * 
 * Truth information from all data products in the form
 * `std::vector<simb::MCTruth>` will be evaluated.
 * It is an error for an event to have no such data product (but they can all be
 * empty, in which case the event does _not_ qualify).
 * 
 * 
 * Configuration
 * ==============
 * 
 * Parameters
 * -----------
 * 
 * The module can be operated in one of two ways: active volumes or generic.
 * In the active volume mode, all TPC active volumes are used for filtering.
 * In the generic mode, all the specified volumes are used for filtering.
 * 
 * * `inActive` (flag, default: `false`): if set to `true`, the "active volume
 *   mode" is activated;
 * * `volumeBoxes` (list of coordinate sets; default: empty): a list of selected
 *   volumes is specified as a collection of box coordinates, each one in the
 *   form of a table:
 *     * `Xmin`: lower _x_ coordinate
 *     * `Xmax`: upper _x_ coordinate
 *     * `Ymin`: lower _y_ coordinate
 *     * `Ymax`: upper _y_ coordinate
 *     * `Zmin`: lower _z_ coordinate
 *     * `Zmax`: upper _z_ coordinate
 *   each coordinate is in centimeters, and in the "world" reference frame
 * * `volumeNames` (list of strings; default: empty): a list of selected volume
 *   patterns is specified as volume names; all volumes in the detector geometry
 *   with the specified names are included in the list (for example, specifying
 *   `volCryostat` _all_ cryostats will be included); each pattern is a regular
 *   expression like parsed by `std::regex` with default flags;
 * * `interactionTypes` (list of integers, optional): if specified, interactions
 *   only qualify if they belong to _any_ of these interaction types;
 *   see `simb::MCNeutrino::InteractionType()`, and `simb::int_type_` definition
 *   for the values and their meanings;
 * * `weakCurrent` (`CC` or `NC`, optional): if specified, interactions
 *   qualify only if they are tagged as charged or neutral currents;
 *   see `icarus::WeakCurrentTypes::parse()`;
 * * `initialState` (list of particle configurations, default: any):
 *   consider only interactions with a composition matching one of the ones
 *   listed here; this parameter is a list of configurations: only if an
 *   interaction does not match _any_ of the configurations, it is excluded;
 *   a configuration is a list of interacting particles, and each interacting
 *   particle is defined by a list of particle identifier (PDG ID). Therefore,
 *   this parameter is effectively a triply nested collection; for example:
 *       
 *       initialState: [
 *         [ # first accepted process:
 *           [ "nu_e", "nu_mu" ], # either a muon or an electron neutrino
 *           [ "e-" ]             # and an electron
 *         ],
 *         [ # second accepted process:
 *           [ "nu_e_bar", "nu_mu_bar" ], # either a muon or an electron antineutrino
 *           [ 1000180400 ]               # and a 40Ar nucleus
 *         ],
 *         [
 *           [ "nu_Re", "nu_Rmu", "nu_Re_bar", "nu_Rmu_bar" ] # any right-handed neutrino
 *         ]
 *       ] # initialState
 *       
 *   will accept only interactions of neutrinos with electrons and antineutrinos
 *   with argon nuclei, or any right-handed neutrino whatever it hits. Particle
 *   names are parsed via ROOT (`TDatabasePDG`). Warnings are emitted if a
 *   _particle ID_ is not known to ROOT; using a _name_ that is not known to
 *   ROOT, instead, is a configuration error.
 *   The initial particles in an event are defined by their status code `0` (see
 *   also `sim::ParticleStatusName()`).
 * * `logCategory` (string, default: `FilterNeutrinosActiveVolume`): name of the
 *   category this module uses to send messages to the message facility.
 * 
 * It is required that at least one among `inActive`, `volumeNames`,
 * `volumeBoxes`, `weakCurrent` and `interactionTypes` are specified.
 * 
 * 
 * Filtering
 * ----------
 * 
 * This is a filter module: its configuration must appear in the `filters`
 * table and to be included in a path in `trigger_path`, and the output module
 * must be told to include that path in `SelectEvents` configuration key; e.g.:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * physics: {
 *   filters: {
 *     eventActive: @local::icarus_FilterNeutrinoActive
 *     # ...
 *   } # filters
 *   # ...
 *   
 *   appliedFilters: [ eventActive ]
 *   
 *   trigger_paths: [ appliedFilters ]
 *   
 * } # physics
 * 
 * 
 * outputs: {
 *   rootoutput: {
 *     module_type:  RootOutput
 *     fileName:    "%ifb_%tc-%p.root"
 *     SelectEvents: [ appliedFilters ]
 *   }
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */
class icarus::simfilter::FilterNeutrinosActiveVolume: public art::EDFilter {
  
  public:
    
    /// Configuration of box volume geometry.
    struct BoxCoordConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<double> Xmin {
        Name("Xmin"),
        Comment("minimum x coordinate of the volume [cm]")
        };
      fhicl::Atom<double> Ymin {
        Name("Ymin"),
        Comment("minimum y coordinate of the volume [cm]")
        };
      fhicl::Atom<double> Zmin {
        Name("Zmin"),
        Comment("minimum z coordinate of the volume [cm]")
        };
      fhicl::Atom<double> Xmax {
        Name("Xmax"),
        Comment("maximum x coordinate of the volume [cm]")
        };
      fhicl::Atom<double> Ymax {
        Name("Ymax"),
        Comment("maximum y coordinate of the volume [cm]")
        };
      fhicl::Atom<double> Zmax {
        Name("Zmax"),
        Comment("maximum z coordinate of the volume [cm]")
        };
      
    }; // struct BoxCoordConfig
    
    
    /// Configuration parameter structure.
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<bool> inActive {
        Name("inActive"),
        Comment("selects only events with interactions in TPC active volumes"),
        false // default
        };
      
      fhicl::OptionalSequence<fhicl::Table<BoxCoordConfig>> volumeBoxes {
        Name("volumeBoxes"),
        Comment("interactions in box volumes by world coordinates [cm]")
        };
      
      fhicl::Sequence<std::string> volumeNames {
        Name("volumeNames"),
        Comment("interactions in box volumes by name (std::regex pattern)"),
        std::vector<std::string>{} // default value
        };
      
      fhicl::Sequence<int> interactionTypes {
        Name("interactionTypes"),
        Comment(
          "require an interaction of one of these types (`simb::MCNeutrino::InteractionType()`)"
          ),
        std::vector<int>{}
        };
      
      fhicl::Sequence<fhicl::Sequence<fhicl::Sequence<std::string>>> initialState {
        Name("initialState"),
        Comment("list of accepted initial states; each entry is an interaction,"
          " as a list of interacting particles; each interacting particle is"
          " a list of particle types, as a PDG ID or particle name; if empty"
          " (default), no requirement on initial state particle flavour is made"
          ),
        std::vector<std::vector<std::vector<std::string>>>{}
        };
      
      fhicl::Atom<std::string> weakCurrent {
        Name("weakCurrent"),
        Comment(
          "require an interaction of charged or neutral current (\"CC\", \"NC\", \"any\")"
          ),
        "any" // default
        };
      
      fhicl::Atom<std::string> logCategory {
        Name("logCategory"),
        Comment("category name for message facility message stream"),
        "FilterNeutrinosActiveVolume" // default
        };
      
    }; // Config
    
    using Parameters = art::EDFilter::Table<Config>;
    
    
    /// Constructor: reads configuration and extracts information from geometry.
    explicit FilterNeutrinosActiveVolume(Parameters const& config);
    
    /// Framework hook: applies the filter.
    virtual bool filter(art::Event& event) override;
    
    /// Framework hook: prints the summary of the passed events.
    virtual void endJob() override;
    
  private:
    
    // --- BEGIN -- Configuration parameters -----------------------------------
    
    /// Volumes for qualifying interactions.
    std::vector<geo::BoxBoundedGeo> fVolumes;
    
    /// List of qualifying interaction types.
    std::vector<int> const fInteractions;
    
    /// Allowed initial states.
    std::vector<details::InteractionComposition_t> const fInitialStates;
    
    icarus::WeakCurrentType const fWeakCurrentType; ///< Selected weak current.
    
    std::string const fLogCategory; ///< Category name for the output stream.
    
    // --- END -- Configuration parameters -------------------------------------
    
    
    // --- BEGIN -- Counters ---------------------------------------------------
    
    std::atomic<unsigned int> fNObserved { 0U }; ///< Number of observed events.
    std::atomic<unsigned int> fNPassed { 0U }; ///< Number of passed events.
    
    // --- END -- Counters -----------------------------------------------------
    
    /// Adds all active volumes of detector into the qualifying volume list.
    void addActiveVolumes();
    
    /// Adds the specified volumes into the qualifying volume list.
    void addVolumeBoxes
      (fhicl::OptionalSequence<fhicl::Table<BoxCoordConfig>> const& boxConfig);
    
    /// Adds all `volName` from geometry into the qualifying volume list.
    unsigned addVolumeByName(std::string const& volumeName);
    
    
    /// Returns whether the interaction described in `truth` qualifies.
    bool qualifying(simb::MCTruth const& truth) const;
    
    /// Returns whether the interaction type is qualifying.
    bool qualifyingInteractionType(int const interactionType) const;

    /// Returns whether the initial state particles qualify.
    bool qualifyingInitialState(std::vector<int> const& initialState) const;

    /// Returns whether the weak current type is qualifying.
    bool qualifyingWeakCurrent(int const CCNC) const;
    
    /// Returns whether the location is among the accepted ones.
    bool qualifyingLocation(geo::Point_t const& location) const;
    
    /// Returns whether `initialState` is compatible with `interaction`.
    bool interactionMatches(
      details::InteractionComposition_t const& interaction,
      std::vector<int> const& initialState
      ) const;

    /// Returns a sorted copy of the specified collection.
    template <typename Coll>
    static Coll sorted(Coll const& coll);
    
    /// Converts `initialState` parameter value into a list of interactions.
    static std::vector<details::InteractionComposition_t> parseInitialState
      (std::vector<std::vector<std::vector<std::string>>> const& interSpecs);
    
    
    /// Converts an interaction from configuration into an object.
    static details::InteractionComposition_t parseInteraction
      (std::vector<std::vector<std::string>> const& interSpec);
    
    /// Converts an interacting particle from configuration into an object.
    static details::InteractingParticle_t parseInteractingParticle
      (std::vector<std::string> const& partSpec);

    /// Returns the PDG ID of a particle with the specified name... or ID.
    static int parseParticleID(std::string const& IDstr);

    /// Returns all PDG ID of initial state particles in `truth` (status `0`).
    static std::vector<int> extractInitialStateParticles
      (simb::MCTruth const& truth);
    
}; // icarus::simfilter::FilterNeutrinosActiveVolume



// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
std::ostream& icarus::simfilter::details::operator<<
  (std::ostream& out, InteractingParticle_t const& particle)
{
  auto const& flavors = particle.flavors;
  if (flavors.size() == 1) {
    out << sim::ParticleName(flavors.front());
  }
  else {
    auto iFlav = begin(flavors);
    auto const fend = end(flavors);
    out << "{";
    if (iFlav != fend) out << " " << sim::ParticleName(*(iFlav++));
    while (iFlav != fend) out << " / " << sim::ParticleName(*(iFlav++));
    out << " }";
  }
  return out;
} // operator<< (InteractingParticle_t)


// -----------------------------------------------------------------------------
std::ostream& icarus::simfilter::details::operator<<
  (std::ostream& out, InteractionComposition_t const& interaction)
{
  auto const& particles = interaction.particles;
  if (particles.empty()) out << "any";
  else {
    auto iPart = begin(particles);
    auto const pend = end(particles);
    out << *iPart;
    while (++iPart != pend) out << " + " << *iPart;
  }
  return out;
} // operator<< (InteractionComposition_t)


// -----------------------------------------------------------------------------
namespace icarus::simfilter::details { class InteractionMatcher; }

/// One-shot algorithm to match interaction specification to an initial state.
class icarus::simfilter::details::InteractionMatcher {
  
    private:
  
  static constexpr int NoParticle = std::numeric_limits<int>::min();
  
  /// Status of the current matching.
  struct MatchStatus {
    
    /// On creation, assigns a particle match; on destruction, releases it.
    struct Assigner {
      
      MatchStatus& matchStatus; ///< Global status of the match.
      
      /// Matched initial state particle.
      int particle = NoParticle;
      int& target; ///< Slot of matched initial state particle.
      
      Assigner(int& part, MatchStatus& status) noexcept
        : matchStatus(status), target(part)
        {
          std::swap(target, particle);
          matchStatus.matched.push_back(particle);
        }
      
      ~Assigner()
        {
          assert(!matchStatus.matched.empty());
          assert(matchStatus.matched.back() == particle);
          matchStatus.matched.pop_back();
          std::swap(particle, target);
        }
      
    }; // Assigner
    
    /// List of available initial state particles. Size is fixed.
    std::vector<int> available;
    /// List of matched initial state particles, in matched order.
    std::vector<int> matched;
    
    MatchStatus(std::vector<int> initialState)
      : available(std::move(initialState))
      {}
    
    Assigner assign [[nodiscard]] (int& particle)
      { return Assigner{ particle, *this }; }
    
    bool completeMatch() const noexcept
      { return matched.size() == available.size(); }
    
  }; // struct MatchStatus
  
  
  MatchStatus matchStatus; ///< Status of the current matching.
  
    public:
  using InitialState_t = std::vector<int>; ///< List of particle flavours.
  
  InteractionMatcher(InitialState_t initialState)
    : matchStatus(std::move(initialState))
    {}
  
  bool match(InteractionComposition_t const& interaction)
    {
      return interaction.empty()
        ? true: matchRemaining(interaction.begin(), interaction.end());
    }
  bool operator() (InteractionComposition_t const& interaction)
    { return match(interaction); }
  
    private:
  using ParticleIter_t = InteractionComposition_t::Particles_t::const_iterator;
  
  /// Match all remaining initial state particles (recursive).
  bool matchRemaining(ParticleIter_t begin, ParticleIter_t end)
    {
      assert(begin != end);
      
      // here we try to match this one
      InteractingParticle_t const& targetParticle = *begin;
      for (int& particle: matchStatus.available) {
        if (particle == NoParticle) continue; // not available any more
        if (!particleFlavourMatches(targetParticle, particle)) continue;
        
        // we matched this role!
        // let's assign the matched role and remove it from the available ones
        auto const assigned = matchStatus.assign(particle); // this is a guard
        auto const next = std::next(begin);
        if (next == end) return true; // happy ending
        
        // try to match the rest of the particles now
        if (matchRemaining(next, end)) return true;
        
      } // for all available particles
      return false;
    } // matchRemaining()
  
  static bool particleFlavourMatches
    (InteractingParticle_t const& particle, int flavor)
    {
      auto const& flavors = particle.flavors;
      return std::find(begin(flavors), end(flavors), flavor) != end(flavors);
    }
  
}; // class icarus::simfilter::details::InteractionMatcher


// -----------------------------------------------------------------------------
// ---  icarus::simfilter::FilterNeutrinosActiveVolume
// -----------------------------------------------------------------------------
icarus::simfilter::FilterNeutrinosActiveVolume::FilterNeutrinosActiveVolume
  (Parameters const& config)
  : art::EDFilter(config)
  , fInteractions(sorted(config().interactionTypes()))
  , fInitialStates(parseInitialState(config().initialState()))
  , fWeakCurrentType(config().weakCurrent())
  , fLogCategory(config().logCategory())
{
  
  { // local scope
    mf::LogInfo log(fLogCategory);
    
    log << "Filter configuration:";
    
    if (!fInteractions.empty()) {
      log << "\n * required one of these " << size(fInteractions)
        << " interaction types:";
      for (int intType: fInteractions)
        log << "\n    - " << sim::TruthInteractionTypeName(intType);
    } // if
    
    if (!fInitialStates.empty()) {
      log << "\n * required one of " << fInitialStates.size()
        << " initial state categories:";
      for (details::InteractionComposition_t const& interaction: fInitialStates)
        log << "\n     * " << interaction;
    } // if
    
    log << "\n * weak current type: " << std::string(fWeakCurrentType);
    
    log << "\nConfiguration of qualifying volumes:";
    
  }
  
  //
  // load volumes from the different sources
  //
  if (config().inActive()) addActiveVolumes();
  
  addVolumeBoxes(config().volumeBoxes);
  
  for (std::string const& volName: config().volumeNames())
    addVolumeByName(volName);
  
  //
  // check that we are doing at least something
  //
  if (fVolumes.empty()
    && fInteractions.empty()
    && fInitialStates.empty()
    && (fWeakCurrentType == icarus::AnyWeakCurrentType)
  ) {
    
    throw art::Exception(art::errors::Configuration)
      << "No filtering action specified"
        " (volume, current, initial state nor interaction type).\n"
      ;
    
  } // if no filter
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::FilterNeutrinosActiveVolume()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterNeutrinosActiveVolume::filter(art::Event& event) {
  
  ++fNObserved;
  
  /*
   * Consider all truth information available in the event.
   * Any record of any truth data product will be enough to pass the event.
   */
  std::vector<art::Handle<std::vector<simb::MCTruth>>> allTruth;
  event.getManyByType(allTruth);
  
  if (allTruth.empty()) { // is this real data?
    throw art::Exception(art::errors::ProductNotFound)
      << event.id() << " has no truth information!\n";
  } // if no truth
  
  mf::LogDebug(fLogCategory)
    << "Event " << event.id() << " (#" << fNObserved << ") has "
    << allTruth.size() << " truth data products.";
  
  for (auto const& handle: allTruth) {
    
    art::InputTag const& tag [[gnu::unused]] = handle.provenance()->inputTag();
    
    std::vector<simb::MCTruth> const& truths = *handle;
    if (truths.empty()) {
      mf::LogTrace(fLogCategory)
        << "No truth records from " << tag.encode() << ": skipped.";
      continue;
    } // if no truth
    
    for (auto const& [ iTruth, truth ]: util::enumerate(truths)) {
      
      mf::LogTrace(fLogCategory)
        << "Processing record [" << (iTruth + 1U) << "/" << truths.size()
        << "] from " << tag.encode();
      
      if (!qualifying(truth)) continue; // neeext!
      
      ++fNPassed;
      mf::LogTrace(fLogCategory) << "Event " << event.id()
        << " (#" << fNObserved << ") passes the filter in virtue of "
        << tag.encode() << " record " << (iTruth + 1U)
        << " (" << fNPassed << "/" << fNObserved << " passed so far)."
        ;
      return true;
      
    } // for truth record
    
  } // for truth data product
  
  mf::LogTrace(fLogCategory) << "Event " << event.id() << " (#" << fNObserved
    << ")  does not pass the filter (" << fNPassed << "/" << fNObserved
    << " passed so far).";
  
  return false;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::filter()


// -----------------------------------------------------------------------------
void icarus::simfilter::FilterNeutrinosActiveVolume::endJob() {
  
  mf::LogInfo log(fLogCategory);
  log
    << "FilterNeutrinosActiveVolume: passed " << fNPassed << " / " << fNObserved
    << " events";
  if (fNObserved > 0U)
    log << " (" << (float(fNPassed) * 100. / fNObserved) << "%)";
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::endJob()


// -----------------------------------------------------------------------------
void icarus::simfilter::FilterNeutrinosActiveVolume::addActiveVolumes() {
  
  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
  
  for (geo::TPCGeo const& TPC: geom.IterateTPCs()) {
    
    geo::BoxBoundedGeo const& box = TPC.ActiveBoundingBox();
    
    mf::LogVerbatim(fLogCategory)
      << "[volume #" << fVolumes.size() << "] active volume from " << TPC.ID()
      << ": [ " << box.Min() << " -- " << box.Max() << " ]";
    
    fVolumes.push_back(box);
    
  } // for all TPCs
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::addActiveVolumes()


// -----------------------------------------------------------------------------
void icarus::simfilter::FilterNeutrinosActiveVolume::addVolumeBoxes
  (fhicl::OptionalSequence<fhicl::Table<BoxCoordConfig>> const& boxConfig)
{
  std::vector<BoxCoordConfig> boxParams;
  
  if (!boxConfig(boxParams)) return;
  
  for (auto const& [ iBox, boxParam ]: util::enumerate(boxParams)) {
    
    geo::BoxBoundedGeo box {
      boxParam.Xmin(), boxParam.Xmax(),
      boxParam.Ymin(), boxParam.Ymax(),
      boxParam.Zmin(), boxParam.Zmax()
      };
    
    mf::LogVerbatim(fLogCategory)
      << "[volume #" << fVolumes.size() << "] box coordinates #" << iBox
      << ": [ " << box.Min() << " -- " << box.Max() << " ]";
    
    fVolumes.push_back(std::move(box));
    
  } // for boxes
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::addVolumeBoxes()


// -----------------------------------------------------------------------------
unsigned int icarus::simfilter::FilterNeutrinosActiveVolume::addVolumeByName
  (std::string const& volumePattern)
{
  
  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
  
  //
  // find the path of all volumes matching the given pattern
  //
  std::regex const namePattern { volumePattern };
  std::vector<geo::GeoNodePath> volumePaths;
  auto const matchMe = [&pattern=namePattern](std::string const& s)
    { std::smatch match; return (std::regex_match(s, match, pattern)); };
  auto const findVolume = [&volumePaths, &patternMatcher=matchMe](auto& path)
    {
      if (patternMatcher(path.current().GetVolume()->GetName()))
        volumePaths.push_back(path);
      return true;
    };
  
  geo::ROOTGeometryNavigator navigator { *(geom.ROOTGeoManager()) };
  
  navigator.apply(findVolume);
  
  if (volumePaths.empty()) {
    throw art::Exception(art::errors::Configuration)
      << "No volume matching '" << volumePattern
      << "' has been found in the detector '" << geom.DetectorName()
      << "'.\n";
  }
  
  //
  // convert each volume into world coordinates and add it to the list
  //
  for (auto const& [ iVolume, path ] : util::enumerate(volumePaths)) {
    
    //
    // find the coordinates of the volume in local coordinates
    //
    TGeoShape const* pShape = path.current().GetVolume()->GetShape();
    auto pBox = dynamic_cast<TGeoBBox const*>(pShape);
    if (!pBox) {
      throw cet::exception("FilterNeutrinosActiveVolume") << "Volume '"
        << path.current().GetName() << "' is a " << pShape->IsA()->GetName()
        << ", not a TGeoBBox.\n";
    }
    
    geo::Point_t const origin
      = geo::vect::makeFromCoords<geo::Point_t>(pBox->GetOrigin());
    geo::Vector_t const diag = {
      std::abs(pBox->GetDX()), std::abs(pBox->GetDY()), std::abs(pBox->GetDZ())
      };
    
    //
    // convert to world coordinates
    //
    
    auto const trans
      = path.currentTransformation<geo::TransformationMatrix>();
    
    geo::Point_t min, max;
    trans.Transform(origin - diag, min);
    trans.Transform(origin + diag, max);
    
    //
    // add to the coordinates
    //
    geo::BoxBoundedGeo box { min, max };
    
    mf::LogVerbatim(fLogCategory)
      << " c* [volume #" << fVolumes.size() << "] volume box '"
      << path.current().GetVolume()->GetName()
      << "' [(" << (iVolume + 1U) << "/" << volumePaths.size()
      << "): [ " << box.Min() << " -- " << box.Max() << " ]";
    
    fVolumes.push_back(std::move(box));
    
  } // for all volume paths
  
  return volumePaths.size();
} // icarus::simfilter::FilterNeutrinosActiveVolume::addVolumeByName()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterNeutrinosActiveVolume::qualifying
  (simb::MCTruth const& truth) const
{
  /*
   * Apply all the needed cuts:
   * * interaction type
   * * current type
   * * location
   * 
   */
  
  //
  // only neutrino record types may qualify:
  //
  if (!truth.NeutrinoSet()) {
    mf::LogTrace(fLogCategory)
      << "Interaction does not qualify because it is not tagged as neutrino.";
    return false;
  }
  
  simb::MCNeutrino const& nuInfo = truth.GetNeutrino();
  simb::MCParticle const& nu = nuInfo.Nu();
  
  //
  // interaction type
  //
  if (!fInteractions.empty() && !qualifyingInteractionType(nuInfo.InteractionType()))
    return false;
  
  //
  // current type
  //
  if ((fWeakCurrentType != icarus::AnyWeakCurrentType) && !qualifyingWeakCurrent(nuInfo.CCNC()))
    return false;
  
  //
  // location
  //
  if (!fVolumes.empty() && !qualifyingLocation({ nu.Vx(), nu.Vy(), nu.Vz() }))
    return false;
  
  //
  // initial state
  //
  if (!fInitialStates.empty() && !qualifyingInitialState(extractInitialStateParticles(truth)))
    return false;
  
  // success, after all
  return true;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::qualifying()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingInteractionType
  (int const interactionType) const
{
  mf::LogTrace log(fLogCategory);
  log
    << "Interaction type: " << sim::TruthInteractionTypeName(interactionType)
    << " (" << interactionType << ")";
  
  bool const pass = std::binary_search
    (begin(fInteractions), end(fInteractions), interactionType);
  log << " => :-" << (pass? ')': '(');
  return pass;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingInteractionType()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingInitialState
  (std::vector<int> const& initialState) const
{
  {
    mf::LogTrace log(fLogCategory);
    log << "Initial state: " << initialState.size() << " particles:";
    for (int const part: initialState) log << " " << sim::ParticleName(part);
  }
  for (details::InteractionComposition_t const& interaction: fInitialStates)
    if (interactionMatches(interaction, initialState)) return true;
  
  return false;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingInitialState()


// -----------------------------------------------------------------------------
/// Returns whether the weak current type is qualifying.
bool icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingWeakCurrent
  (int const CCNC) const
{
  using icarus::WeakCurrentType;
  
  mf::LogTrace log(fLogCategory);
  log << "Interaction current: " << sim::TruthCCNCname(CCNC)
    << " (" << CCNC << ")";
  
  bool pass = false;
  switch (fWeakCurrentType) {
    case WeakCurrentType::CC:  pass = (CCNC == simb::kCC); break;
    case WeakCurrentType::NC:  pass = (CCNC == simb::kNC); break;
    case WeakCurrentType::any: pass = true; break;
  } // switch
  
  log << " => :-" << (pass? ')': '(');
  return pass;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingWeakCurrent()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingLocation
  (geo::Point_t const& location) const
{
  
  mf::LogTrace log(fLogCategory);
  log << "Interaction location: " << location << " cm";
  
  for (auto const& [ iBox, box ]: util::enumerate(fVolumes)) {
    if (!box.ContainsPosition(location)) continue;
    
    log << " => in volume #" << iBox
      << " [ " << box.Min() << " -- " << box.Max() << " ] => :-)";
    
    return true;
  } // for
  
  log << " => :-(";
  return false;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::qualifyingLocation()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterNeutrinosActiveVolume::interactionMatches(
  details::InteractionComposition_t const& interaction,
  std::vector<int> const& initialState
) const {
  
  mf::LogTrace log(fLogCategory);
  log << "  - initial state " << interaction << "?";
  
  // sanity checks first:
  if (interaction.empty()) {
    log << " yes (no constraint)"; // ??
    return false;
  }
  if (interaction.size() > initialState.size()) {
    log << " no: not enough particles";
    return false;
  }
  
  bool const matched = details::InteractionMatcher{ initialState }(interaction);
  log << " => :-" << (matched? ')': '(');
  return matched;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::interactionMatches()


// -----------------------------------------------------------------------------
template <typename Coll>
Coll icarus::simfilter::FilterNeutrinosActiveVolume::sorted(Coll const& coll) {
  
  // copy, sort, return
  auto sortedColl { coll };
  std::sort(util::begin(sortedColl), util::end(sortedColl));
  return sortedColl;
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::sorted()



// -----------------------------------------------------------------------------
auto icarus::simfilter::FilterNeutrinosActiveVolume::parseInitialState
  (std::vector<std::vector<std::vector<std::string>>> const& interSpecs)
  -> std::vector<details::InteractionComposition_t>
{
  std::vector<details::InteractionComposition_t> interactions;
  for (auto const& spec: interSpecs)
    interactions.push_back(parseInteraction(spec));
  return interactions;
} // icarus::simfilter::FilterNeutrinosActiveVolume::parseInitialState()


// -----------------------------------------------------------------------------
auto icarus::simfilter::FilterNeutrinosActiveVolume::parseInteraction
  (std::vector<std::vector<std::string>> const& interSpec)
  -> details::InteractionComposition_t
{
  details::InteractionComposition_t interaction;
  for (auto const& spec: interSpec)
    interaction.particles.push_back(parseInteractingParticle(spec));
  return interaction;
} // icarus::simfilter::FilterNeutrinosActiveVolume::parseInteraction()


// -----------------------------------------------------------------------------
auto icarus::simfilter::FilterNeutrinosActiveVolume::parseInteractingParticle
  (std::vector<std::string> const& partSpec)
  -> details::InteractingParticle_t
{
  details::InteractingParticle_t particle;
  for (std::string const& partIDstr: partSpec) {
    try {
      particle.flavors.push_back(parseParticleID(partIDstr));
    }
    catch (cet::exception const& e) {
      throw art::Exception{art::errors::Configuration, "", e }
        << "Error parsing initial state configuration particles.\n";
    }
  }
  return particle;
} // icarus::simfilter::FilterNeutrinosActiveVolume::parseInteractingParticle()


// -----------------------------------------------------------------------------
int icarus::simfilter::FilterNeutrinosActiveVolume::parseParticleID
  (std::string const& IDstr)
{
  static std::regex const IntegerPattern
    { R"([[:space:]]*[+-]?[[:digit:]]+[[:space:]]*)" };
  
  if (std::regex_match(IDstr, IntegerPattern)) {
    try {
      return std::stoi(IDstr);
    }
    catch (std::exception const& e) {
      throw cet::exception("parseParticleID")
        << "Logic error: '" << IDstr
        << "' looked like a number, but `std::stoi()` can't convert it (\""
        << e.what() << "\").\n";
    }
  }
  else { // not an integer... must be a particle name; ROOT?
    TDatabasePDG const& Particles = *(TDatabasePDG::Instance());
    
    TParticlePDG const* part = Particles.GetParticle(IDstr.c_str());
    if (!part) {
      throw cet::exception("parseParticleID")
        << "No particle with name '" << IDstr << "' is known to ROOT.\n";
    }
    return part->PdgCode();
  }
  
} // icarus::simfilter::FilterNeutrinosActiveVolume::parseParticleID()


// -----------------------------------------------------------------------------
std::vector<int>
icarus::simfilter::FilterNeutrinosActiveVolume::extractInitialStateParticles
  (simb::MCTruth const& truth)
{
  static constexpr int stInitialParticle = 0; // status: initial state particle
  
  std::vector<int> particles;
  for (auto const iPart: util::counter(truth.NParticles())) {
    simb::MCParticle const& part = truth.GetParticle(iPart);
    if (part.StatusCode() != stInitialParticle) continue;
    particles.push_back(part.PdgCode());
  } // for
  
  return particles;
} // icarus::simfilter::FilterNeutrinosActiveVolume::extractInitialStateParticles()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::simfilter::FilterNeutrinosActiveVolume)


// -----------------------------------------------------------------------------

