/**
 * @file   AnodeCathodeCrossingClusters_module.cc
 * @brief  Module: selects clusters spanning a full anode/cathode range.
 * @author Christian Farnese, Weskey Ketchum, Olivia Bitter, Gianluca Petrillo
 * @date   January 12, 2021
 * 
 * The algorithm of this module is taken from LAr purity analysis code
 * (`ICARUSPurityDQMnew`).
 * 
 */

// ICARUS libraries
#include "icaruscode/Analysis/Algorithms/AnodeToCathodeClusterAlg.h"

// LArSoft libraries
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#ifdef USE_DETECTORTIMINGS
# include "lardataalg/DetectorInfo/DetectorTimings.h"
# include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // detinfo::timescales
#else
# include "lardataalg/DetectorInfo/DetectorClocksData.h"
#endif // USE_DETECTORTIMINGS
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// framework libraries
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// C/C++ standard libraries
#include <vector>
#include <array>
#include <memory> // std::make_unique()
#include <algorithm> // std::accumulate(), ...
#include <mutex> // std::mutex, std::lock_guard
#include <atomic>
#include <utility> // std::move()



// -----------------------------------------------------------------------------
namespace icarus {
  class AnodeCathodeCrossingClusters;
  class ReasonTracker;
}

/// Keeps track of rejection reasons. Thread-safe.
class icarus::ReasonTracker {
  /*
   * Thread-safety note
   * -------------------
   * 
   * Access to the single counter values is built-in thread-safe as each counter
   * is statically allocated, unmoveable and atomic.
   * Access to the whole set of counters is made consistent via locking.
   */
    public:
  
  using reason = icarus::AnodeToCathodeClusterAlg::reason; ///< Shortcut.
  
  /// Number of reasons tracked by this class.
  static constexpr std::size_t NReasons
    = static_cast<std::size_t>(reason::NReasons);
  
  using Counters_t = std::array<unsigned int, NReasons>; /// Type of counters.
  
  // --- BEGIN -- Constructors -------------------------------------------------
  
  /// Default constructor: all counters are set to 0.
  ReasonTracker();
  
  /// Copy constructor; thread-safe.
  ReasonTracker(ReasonTracker const& from);
  
  /// Move constructor; actaually copies.
  ReasonTracker(ReasonTracker&& from);
  
  /// Copy assignment; thread-safe.
  ReasonTracker& operator= (ReasonTracker const& from);
  
  /// Move operator: actually copies.
  ReasonTracker& operator= (ReasonTracker&& from);
  
  // --- END -- Constructors ---------------------------------------------------
  
  
  /// Records `n` entries for the specified reason `r`.
  void record(reason r, unsigned int n = 1U);
  
  
  // --- BEGIN -- Access -------------------------------------------------------
  /// @name Access to counters
  /// @{
  
  /// Returns the number of reasons being tracked.
  std::size_t nReasons() const;
  
  /// Returns the current value of all counters.
  Counters_t counters() const;
  
  /// Returns the current value of the specified counter.
  unsigned int counter(reason r) const;
  
  /// Returns the name of the counter.
  static const char* name(reason r);
  
  /// Returns the currently total recorded reasons.
  unsigned int total() const;
  
  /// @}
  // --- END -- Access ---------------------------------------------------------
  
    private:
  
  /// All reason counters.
  mutable std::array<std::atomic<unsigned int>, NReasons> fCounters;
  
  mutable std::mutex fCounterAccess; ///< Access control for the whole counters.
  
  /// Locks the counter access and returns a guard for it.
  std::lock_guard<std::mutex> lockCounters() const;
  
  /// Copies values from `counters`. Use `lockCounters()` for thread-safety.
  void copyCounters(Counters_t const& counters);
  
}; // icarus::ReasonTracker


// -----------------------------------------------------------------------------
/**
 * @brief Selects clusters spanning a full anode/cathode range.
 * 
 * This module selects clusters based on their spanning a full range of drift
 * time.
 * The clusters undergo a hit selection; the cluster in output are recomputed
 * based only on the hits that are selected.
 * 
 * 
 * Algorithms
 * ===========
 * 
 * This module is a _art_ wrapper for the algorithm
 * `icraus::AnodeToCathodeClusterAlg`.
 * The selection of the clusters and the filtering of the hits are performed
 * by that algorithm.
 * Clusters are considered one by one and independently.
 * This module converts the each of the selected input clusters into a new
 * cluster object, created from the hits returned by the aforementioned
 * algorithm. Cluster creation algorithm is the "standard"
 * `cluster::StandardClusterParamsAlg`, and in addition the cluster parameters
 * are chosen as follows:
 * 
 * * **cluster ID**, **plane** and view are the same as the original cluster
 * * start point is associated to the hit with the earliest peak time
 *   (presumably, anode), the end point to the hit with the latest peak time
 *   (supposedly, cathode)
 * * wire extremes: wire number of the earliest and latest hits, projected on
 *   the cluster plane; uncertainty parameter is set to 0
 * * time extremes: time (in TDC ticks) of the peak of the earliest and latest
 *   hits; uncertainty is the uncertainty on those hit peak times;
 * 
 * 
 * Input
 * ======
 * 
 * * clusters (`recob::Cluster`) and their associated hits (`recob::Hit`);
 *   effectively, we access the cluster/hits association data products
 * 
 * 
 * Output
 * =======
 * 
 * Clusters and their associations are produced in a single data product which
 * includes all the input collections. To have one output collection per input
 * collection, multiple instances of this module will have to be run instead.
 * 
 * * `std::vector<recob::Cluster>` collection of the selected clusters;
 *   the cluster relative order is the same as in the original collection;
 * * `art::Assns<recob::Cluster, recob::Hit>>` associations of the selected
 *     clusters to the hits used to make them; the associations are promised to
 *     fulfil the
 *     @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association".
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * * **Clusters** (list of input tags, mandatory): tags of the cluster data
 *     products to be processed
 * * **MinHits** (integral, default: `30`): do not consider clusters of fewer
 *     than this number of hits
 * * **MinWires** (integral, default: `100`): do not consider clusters covering
 *     fewer than this number of wires (smallest to largest hit channel ID)
 * * **MinTicks** (integral, default: `1100`): do not consider clusters covering
 *     fewer than this number of TPC ticks (smallest/earliest to largest/latest
 *     TPC hit TDC)
 * * **MaxHitClusterDistance** (real, default: `3.0`) hits further than this
 *     much from cluster axis are ignored [cm]
 * * **MinCloseHits** (integral, default: `2`) after purging hits that are too
 *     far, keep only clusters with at least this many hits
 * * **MinDriftDistance** (real, default: `0`) clusters spanning shorter than
 *     this drift distance [cm] are discarded; this is the projection on drift
 *     direction of the size of the cluster
 * * **MaxDriftDistance** (real, optional) if specified, clusters spanning
 *     longer than this drift distance [cm] are discarded; the definition is the
 *     same as for `MinDriftDistance`
 * * **OnlyPlanes** (list of integral numbers): if specified, only
 *     clusters lying on a wire plane whose number is in this list will be
 *     considered, and all others will be discarded; if the parameter is omitted
 *     or the specified list is empty, clusters from all planes are considered.
 *     Note that cluster from all planes with matching plane number, i.e. from
 *     all TPC and all cryostats in the detector, are considered
 * * **ClusterIDoffset** (integer, default: `100000`): output clusters are given
 *     the same ID (`recob::Cluster::ID()`) as the original ones, to allow easy
 *     matching; if multiple input data products are specified, though, this may
 *     create ambiguity; this parameter will add `ClusterIDoffset` to the ID of
 *     all clusters from the _second_ configured input cluster collection
 *     (`Clusters`), twice `ClusterIDoffset` to the third collection, and so
 *     forth; this number should be chosen high enough to avoid conflicts
 * * **LogCategory** (text, default: "AnodeToCathodeClusterAlg"):
 *     category tag used to output the messages from the clustering algorithm
 *     (not the module itself) into the message facility service;
 *     see also `ModuleLogCategory`
 * * **ModuleLogCategory** (text, default: "AnodeCathodeCrossingClusters"):
 *     category tag used to output the messages from this module into the
 *     message facility service; see also `LogCategory`
 * 
 * 
 * Service requirements
 * =====================
 * 
 * * `DetectorPropertiesService`
 * * `DetectorClocksService`
 * * `Geometry`
 * 
 * 
 * Technical note: multithreading readiness
 * =========================================
 * 
 * This module is **not** supporting multithread.
 * Steps known to be taken to enable it:
 * 
 * * algorithm needs to become local to `produce()` instead of class member;
 * * access to `icarus::ReasonTracker` needs to be protected (via mutex)
 * 
 */
class icarus::AnodeCathodeCrossingClusters: public art::EDProducer {
  
    public:
  
  /// Algorithm configuration (see module documentation for details).
  struct AlgoConfig {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<unsigned int> MinHits {
      Name("MinHits"),
      Comment("do not consider clusters of fewer than this number of hits"),
      30U
      };
    
    fhicl::Atom<unsigned int> MinWires {
      Name("MinWires"),
      Comment("cluster must cover at least these many wires/channels end to end"),
      100U
      };
    
    fhicl::Atom<unsigned int> MinTicks {
      Name("MinTicks"),
      Comment("cluster must cover at least these many time ticks end to end"),
      1100U
      };
    
    fhicl::Atom<double> MaxHitClusterDistance {
      Name("MaxHitClusterDistance"),
      Comment("hits further than this much from cluster axis are ignored [cm]"),
      3.0
      };
    
    fhicl::Atom<unsigned int> MinCloseHits {
      Name("MinCloseHits"),
      Comment("minimum number of hits closer than MaxHitClusterDistance"),
      2U
      };
    
    fhicl::Atom<double> MinDriftDistance {
      Name("MinDriftDistance"),
      Comment("clusters spanning less than this space in drift are discarded [cm]"),
      0.0
      };
    
    fhicl::OptionalAtom<double> MaxDriftDistance {
      Name("MaxDriftDistance"),
      Comment("clusters spanning more than this space in drift are discarded [cm]")
      };
    
    fhicl::Sequence<geo::PlaneID::PlaneID_t> OnlyPlanes {
      Name("OnlyPlanes"),
      Comment
        ("Consider only clusters on wire planes with the specified number(s)"),
      std::vector<geo::PlaneID::PlaneID_t>{}
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("tag for logging into message facility"),
      "AnodeToCathodeClusterAlg"
      };
    
  }; // struct AlgoConfig
  
  /// Module configuration (see module documentation for details).
  struct Config: public AlgoConfig {
    
    fhicl::Sequence<art::InputTag> Clusters {
      Name("Clusters"),
      Comment("list of tags of the input cluster data products")
      };
    
    fhicl::Atom<int> ClusterIDoffset {
      Name("ClusterIDoffset"),
      Comment(
        "offset to add to ID of clusters from input collections after the first one"
        ),
        100000U
      };
    
    fhicl::Atom<std::string> ModuleLogCategory {
      Name("ModuleLogCategory"),
      Comment("tag for logging into message facility (module messages only)"),
      "AnodeCathodeCrossingClusters"
      };
    
  }; // struct Config
  
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  /// Constructor: extracts module configuration, and little more.
  AnodeCathodeCrossingClusters(Parameters const& config);
  
  
  // --- BEGIN -- Framework hooks ----------------------------------------------
  /// @name Framework hooks
  /// @{
  
  void beginJob() override;
  
  void produce(art::Event& evt) override;
  
  void endJob() override;
  
  /// @}
  // --- END -- Framework hooks ------------------------------------------------
  
  
  
    private:
  
  /// Wrapper record for `processClusterCollection()` configuration.
  struct ClusterProcessParams_t {
    util::GeometryUtilities const gser;
    cluster::ClusterParamsImportWrapper<cluster::StandardClusterParamsAlg>
      clusterParamAlg;
    art::PtrMaker<recob::Cluster> const makeClusterPtr;
    icarus::ReasonTracker& reasonTracker;
    int clusterIDoffset = 0;
  }; // struct ClusterProcessParams_t

  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  std::vector<art::InputTag> const fClusterTags; ///< Cluster input tags.
  
  int const fClusterIDoffset; ///< Base offset for ID of clusters.
  
  std::string const fLogCategory; ///< Log category for module messages.
  
  // --- END -- Configuration parameters ---------------------------------------
  
  
  // --- BEGIN -- Cached services ----------------------------------------------
  
  geo::GeometryCore const& fGeom;
#ifdef USE_DETECTORTIMINGS
  detinfo::DetectorTimings const fDetTimings;
#else // !USE_DETECTORTIMINGS
  detinfo::DetectorClocksData const fDetClocks;
#endif // ?USE_DETECTORTIMINGS
  detinfo::DetectorPropertiesData const fDetProp;
  
  // --- END -- Cached services ------------------------------------------------
  
  /// Algorithm being applied.
  icarus::AnodeToCathodeClusterAlg fAnodeToCathodeFilter;
  
  icarus::ReasonTracker fReasons; ///< Tracking of rejection reasons.
  
  
  /// Processes a single cluster data product.
  /// @return number of selected and of processed clusters, in a `pair`
  std::pair<unsigned int, unsigned int> processClusterCollection(
    art::Assns<recob::Cluster, recob::Hit> const& clustersAndHits,
    std::vector<recob::Cluster>& selectedClusters,
    art::Assns<recob::Cluster, recob::Hit>& selectedClusterHits,
    ClusterProcessParams_t& params
    );

  /**
   * @brief Returns which wire number/coordinate `wire` would have on `plane`.
   *
   * Wire orientation compatibility between `wire` and `plane` is assumed
   * and not checked.
   */
  double equivalentWireOn
    (geo::WireID const& wire, geo::PlaneID const& plane) const;
  
  /// Fills and returns dynamic parameters.
  icarus::AnodeToCathodeClusterAlg::SetupParams_t makeSetupParams
    (geo::PlaneID const& plane) const;
  
  /// Fills and returns a `ConfigParams_t` from the module configuration.
  static icarus::AnodeToCathodeClusterAlg::ConfigParams_t makeConfigParams
    (Config const& config);
  
}; // icarus::AnodeCathodeCrossingClusters



// -----------------------------------------------------------------------------
// ---  icarus::ReasonTracker implementation
// -----------------------------------------------------------------------------
inline icarus::ReasonTracker::ReasonTracker()
  { for (auto& counter: fCounters) counter = 0U; }


inline icarus::ReasonTracker::ReasonTracker(ReasonTracker const& from)
  { copyCounters(from.counters()); }


inline icarus::ReasonTracker::ReasonTracker(ReasonTracker&& from)
  : ReasonTracker(from) {}


inline icarus::ReasonTracker& icarus::ReasonTracker::operator=
  (ReasonTracker const& from)
{
  auto const counters = from.counters(); // copy source before locking
  auto const guard { lockCounters() }; copyCounters(counters);
  return *this;
}


inline icarus::ReasonTracker& icarus::ReasonTracker::operator=
  (ReasonTracker&& from)
  { return this->operator=(from); }


inline void icarus::ReasonTracker::record(reason r, unsigned int n /* = 1U */)
  { fCounters[static_cast<std::size_t>(r)] += n; }


inline std::size_t icarus::ReasonTracker::nReasons() const
  { return fCounters.size(); }


inline auto icarus::ReasonTracker::counters() const -> Counters_t {
  Counters_t counters;
  auto const guard { lockCounters() };
  std::copy(fCounters.begin(), fCounters.end(), counters.begin());
  return counters;
} // icarus::ReasonTracker::counters()


inline unsigned int icarus::ReasonTracker::counter(reason r) const
  { return static_cast<unsigned int>(fCounters[static_cast<std::size_t>(r)]); }


inline const char* icarus::ReasonTracker::name(reason r)
  { return icarus::AnodeToCathodeClusterAlg::explain(r); }


inline unsigned int icarus::ReasonTracker::total() const {
  auto const guard { lockCounters() };
  return std::accumulate(fCounters.cbegin(), fCounters.cend(), 0U); 
} // icarus::ReasonTracker::total()


inline std::lock_guard<std::mutex> icarus::ReasonTracker::lockCounters() const
  { return std::lock_guard{ fCounterAccess }; }


inline void icarus::ReasonTracker::copyCounters(Counters_t const& counters)
  { std::copy(counters.begin(), counters.end(), fCounters.begin()); }


// -----------------------------------------------------------------------------
// ---  Module implementation
// -----------------------------------------------------------------------------
icarus::AnodeCathodeCrossingClusters::AnodeCathodeCrossingClusters
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration parameters
  , fClusterTags(config().Clusters())
  , fClusterIDoffset(config().ClusterIDoffset())
  , fLogCategory(config().ModuleLogCategory())
  // services
  , fGeom(*lar::providerFrom<geo::Geometry>())
#ifdef USE_DETECTORTIMINGS
  , fDetTimings{detinfo::makeDetectorTimings
      (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
    }
  , fDetProp(
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()
      ->DataForJob(fDetTimings.clockData())
    )
#else // !USE_DETECTORTIMINGS
  , fDetClocks
      (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
  , fDetProp(
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()
      ->DataForJob(fDetClocks)
    )
#endif // ?USE_DETECTORTIMINGS
  , fAnodeToCathodeFilter(makeConfigParams(config()))
{
  
  // consume declarations
  for (art::InputTag const& tag: fClusterTags)
    consumes<art::Assns<recob::Cluster, recob::Hit>>(tag);
  
  // produce declarations
  produces<std::vector<recob::Cluster>>();
  produces<art::Assns<recob::Cluster, recob::Hit>>();
  
} // icarus::AnodeCathodeCrossingClusters::AnodeCathodeCrossingClusters()


// -----------------------------------------------------------------------------
void icarus::AnodeCathodeCrossingClusters::beginJob() {
  
  mf::LogInfo(fLogCategory)
    << "Algorithm configuration:\n" << fAnodeToCathodeFilter.config();
  
} // icarus::AnodeCathodeCrossingClusters::beginJob()


// -----------------------------------------------------------------------------
void icarus::AnodeCathodeCrossingClusters::produce(art::Event& event) {
  
  // so far, we do not distinguish between wire planes;
  // different parameters might be employed for each plane type though
  
  // this is just a header to locate the event ID in the log...
  mf::LogDebug(fLogCategory)
    << "Processing clusters from " << fClusterTags.size()
    << " cluster collections";
  
  ClusterProcessParams_t clusterProcessParams {
      util::GeometryUtilities{ fGeom, fDetClocks, fDetProp } // gser
    , cluster::ClusterParamsImportWrapper<cluster::StandardClusterParamsAlg>{}
                                                             // clusterParamAlg
    , art::PtrMaker<recob::Cluster>{ event }                 // makeClusterPtr
    , fReasons                                               // reasonTracker
    , 0                                                      // clusterIDoffset
    };

  std::vector<recob::Cluster> selectedClusters;
  art::Assns<recob::Cluster, recob::Hit> selectedClusterHits;
  
  unsigned int totalSelectedClusters = 0U;
  unsigned int totalClusters = 0U;
  for (art::InputTag const& clusterTag: fClusterTags) {
    //
    // read clusters
    //
    /*
     * We pick the clusters via their associations to hits.
     * We assume that the associations follow LArSoft recommended standard
     * and use `util::associated_groups_with_left()` to access cluster and hits
     * together.
     * In addition, we merge the associations from all the cluster data products
     * that we process together.
     * Each cluster is processed independently anyway.
     */
    
    auto const& inputClustersAndHits
       = event.getByLabel<art::Assns<recob::Cluster, recob::Hit>>(clusterTag);
    
    //
    // delegate processing and filling of output
    //
    auto [ nSelectedClusters, nClusters ] = processClusterCollection(
      inputClustersAndHits,
      selectedClusters, selectedClusterHits, clusterProcessParams
      );
    
    totalSelectedClusters += nSelectedClusters;
    totalClusters += nClusters;
    
    if (fClusterTags.size() > 1U) {
      mf::LogVerbatim(fLogCategory)
        << "Selected " << nSelectedClusters << "/" << nClusters
        << " '" << clusterTag.encode() << "' clusters";
    }
    
    clusterProcessParams.clusterIDoffset += fClusterIDoffset; // for next input
    
  } // for cluster collections
  
  mf::LogInfo(fLogCategory)
    << "Selected " << selectedClusters.size() << "/" << totalClusters
    << "' clusters";
  
  event.put(std::make_unique<std::vector<recob::Cluster>>
    (std::move(selectedClusters)));
  event.put(std::make_unique<art::Assns<recob::Cluster, recob::Hit>>
    (std::move(selectedClusterHits)));
  
} // icarus::AnodeCathodeCrossingClusters::produce()


// -----------------------------------------------------------------------------
auto icarus::AnodeCathodeCrossingClusters::processClusterCollection(
  art::Assns<recob::Cluster, recob::Hit> const& clustersAndHits,
  std::vector<recob::Cluster>& selectedClusters,
  art::Assns<recob::Cluster, recob::Hit>& selectedClusterHits,
  ClusterProcessParams_t& params
) -> std::pair<unsigned int, unsigned int> {
  
  // so far, we do not distinguish between wire planes;
  // different parameters might be employed for each plane type though
  
  // this is a sequence of pairs: cluster/sequence of hits (all via art::Ptr)
  auto const clustersWithHits
     = util::associated_groups_with_left(clustersAndHits);
  
  //
  // prepare output data structures
  //
  // used for creating clusters:
  using reason = icarus::AnodeToCathodeClusterAlg::reason;
  
  // to support offsets we need to be a bit creative... track cluster product ID
  unsigned int nSelectedClusters = 0U;
  unsigned int totalClusters = 0U;
  for (auto const& [ cluster, hits ]: clustersWithHits) {
    ++totalClusters;
    
    if (!cluster->hasPlane()) {
      mf::LogWarning(fLogCategory)
        << "Cluster " << cluster->ID()
        << " skipped because its plane is not known.";
      continue;
    }
    fAnodeToCathodeFilter.setup(makeSetupParams(cluster->Plane()));
    
    // run the algorithm here
    auto const [ response, selectedHits ]
      = fAnodeToCathodeFilter(cluster, hits);
    params.reasonTracker.record(response);
    if (response != reason::Accepted) {
      mf::LogTrace log(fLogCategory);
      log << "Cluster " << cluster->ID();
      if (cluster->hasPlane()) log << " on " << cluster->Plane();
      log << " rejected: " << fAnodeToCathodeFilter.explain(response);
      continue;
    }
    
    /*
     * The cluster is assigned a start and stop coordinate matching the ones
     * of the earliest and latest hit. Somehow arbitrary, but consistent.
     */
    ++nSelectedClusters;
    params.clusterParamAlg.ImportHits(params.gser, selectedHits);
    auto const [ earliestHit, latestHit ] = std::minmax_element(
      selectedHits.begin(), selectedHits.end(),
      [](auto const& A, auto const& B){ return A->PeakTime() < B->PeakTime(); }
      );
    geo::PlaneID const plane = cluster->Plane();
    double const earliestWire
      = equivalentWireOn((*earliestHit)->WireID(), plane);
    double const latestWire = equivalentWireOn((*latestHit)->WireID(), plane);
    cluster::ClusterCreator selectedCluster {
        params.gser
      , params.clusterParamAlg                  // algo
      , static_cast<float>(earliestWire)        // start_wire
      , 0.                                      // sigma_start_wire
      , (*earliestHit)->PeakTime()              // start_tick
      , (*earliestHit)->SigmaPeakTime()         // sigma_start_tick
      , static_cast<float>(latestWire)          // end_wire
      , 0.                                      // sigma_end_wire,
      , (*latestHit)->PeakTime()                // end_tick
      , (*latestHit)->SigmaPeakTime()           // sigma_end_tick
      , cluster->ID() + params.clusterIDoffset  // ID
      , cluster->View()                         // view
      , plane                                   // plane
      , recob::Cluster::Sentry                  // sentry
      };
    
    // save the result!
    selectedClusterHits.addMany
      (params.makeClusterPtr(selectedClusters.size()), selectedHits);
    selectedClusters.push_back(selectedCluster.move());
    
  } // for
  
  return { nSelectedClusters, totalClusters };
} // icarus::AnodeCathodeCrossingClusters::processClusterCollection()


// -----------------------------------------------------------------------------
void icarus::AnodeCathodeCrossingClusters::endJob() {
  
  using reason = icarus::AnodeToCathodeClusterAlg::reason;
  
  mf::LogInfo log(fLogCategory);
  log << "Summary of " << fReasons.total() << " evaluated clusters:";
  for (auto [ r, count ]: util::enumerate(fReasons.counters())) {
    log << "\n" << std::setw(9) << count << "  "
      << fReasons.name(static_cast<reason>(r));
  } // for
  
} // icarus::AnodeCathodeCrossingClusters::endJob()


// -----------------------------------------------------------------------------
icarus::AnodeToCathodeClusterAlg::SetupParams_t
icarus::AnodeCathodeCrossingClusters::makeSetupParams
  (geo::PlaneID const& plane) const
{
  
  return {
      fGeom.WirePitch(plane)    // wirePitch
#ifdef USE_DETECTORTIMINGS
    , fDetTimings.ClockPeriodFor<detinfo::timescales::TPCelectronics_time>().value()
                     // tickDuration (TPCelectronics_time scale is microseconds)
#else // !USE_DETECTORTIMINGS
    , fDetClocks.TPCClock().TickPeriod() // tickDuration 
#endif // ?USE_DETECTORTIMINGS
    , fDetProp.DriftVelocity()  // driftVelocity
    };
  
} // icarus::AnodeCathodeCrossingClusters::makeSetupParams()


// -----------------------------------------------------------------------------
icarus::AnodeToCathodeClusterAlg::ConfigParams_t
icarus::AnodeCathodeCrossingClusters::makeConfigParams
  (icarus::AnodeCathodeCrossingClusters::Config const& config)
{
  icarus::AnodeToCathodeClusterAlg::ConfigParams_t params;
  
  params.minHits = config.MinHits();
  params.minChannelSpan = config.MinWires();
  params.minTimeSpan = config.MinTicks();
  params.maxHitClusterDistance = config.MaxHitClusterDistance();
  params.minCloseHits = config.MinCloseHits();
  params.minDriftSpan = config.MinDriftDistance();
  config.MaxDriftDistance(params.maxDriftSpan); // stays huge if it was omitted
  params.logCategory = config.LogCategory();
  params.acceptedPlanes = config.OnlyPlanes();
  
  std::sort(params.acceptedPlanes.begin(), params.acceptedPlanes.end());
  
  return params;
} // icarus::AnodeCathodeCrossingClusters::makeConfigParams()


// -----------------------------------------------------------------------------
double icarus::AnodeCathodeCrossingClusters::equivalentWireOn
  (geo::WireID const& wire, geo::PlaneID const& plane) const
{
  // invalid something, do nothing
  if (!plane || !wire) return static_cast<double>(wire.Wire);
  
  geo::PlaneID const& wirePlane = wire.asPlaneID();
  
  // wire is on the same plane, nothing do be done
  if (wirePlane == plane) return static_cast<double>(wire.Wire);
  
  auto const& targetWire = fGeom.Wire(wire);
  auto const& refPlane = fGeom.Plane(plane);
  
  return refPlane.WireCoordinate(targetWire.GetCenter<geo::Point_t>());
  
} // icarus::AnodeCathodeCrossingClusters::equivalentWireOn()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::AnodeCathodeCrossingClusters)


// -----------------------------------------------------------------------------
