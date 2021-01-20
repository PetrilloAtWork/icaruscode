/**
 * @file   icaruscode/Analysis/Algorithms/AnodeToCathodeClusterAlg.h
 * @brief  Algorithm: selects clusters spanning a full anode/cathode range.
 * @author Christian Farnese, Weskey Ketchum, Olivia Bitter, Gianluca Petrillo
 * @date   January 12, 2021
 * @see    icaruscode/Analysis/Algorithms/AnodeToCathodeClusterAlg.cxx
 * 
 * The algorithm of this module is taken from LAr purity analysis code
 * (`ICARUSPurityDQMnew`).
 * 
 */

#ifndef ICARUSCODE_ANALYSIS_ALGORITHMS_ANODETOCATHODECLUSTERALG_H
#define ICARUSCODE_ANALYSIS_ALGORITHMS_ANODETOCATHODECLUSTERALG_H


// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::PlaneID

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <tuple>
#include <utility> // std::move(), std::pair<>
#include <limits> // std::numeric_limits<>



// -----------------------------------------------------------------------------
namespace icarus {
  //
  // definitions in this header
  //
  class AnodeToCathodeClusterAlg;
  
  //
  // forward declarations
  //
  class HitCoordinateBuilder;
  
} // namespace icarus

/**
 * @brief Anode-to-cathode cluster selection algorithm.
 * 
 * The algorithm tests a cluster to determine if it matches the criteria for
 * anode-to-cathode crossing.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::vector<std::pair<recob::Cluster const*, std::vector<recob::Hit const*>>>
 *   anodeToCathode; // ...
 * 
 * AnodeToCathodeClusterAlg anodeToCathode{{}}; // default configuration
 * unsigned int nCrossingClusters = 0U;
 * for (auto const& [ cluster, hits ]: ClustersAndHits) {
 *   
 *   if (anodeToCathode(cluster, hits)) ++nCrossingClusters;
 *   
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 * Algorithm details
 * ------------------
 * 
 * This has started as an implementation of the selection algorithm originally
 * from `ICARUSPurityDQMnew` module (Christian Farnese, Weskey Ketchum,
 * Olivia Bitter) for purity studies.
 * 
 * The algorithm does not reconstruct or discriminate the direction of the
 * cluster (anode-to-cathode or cathode-to-anode).
 * 
 *  * clusters are required to have minimum extension in both wire (channel)
 *    space (`minChannelSpan`) and time tick (drift) space (`minTimeSpan`)
 *  * points are removed from the cluster when too far (`MaxHitClusterDistance`)
 *    from the cluster axis (defined as the unweighted linear fit of the
 *    wire/drift positions its hits) one by one;
 *  * clusters with less than `minCloseHits` surviving hits are rejected;
 * 
 * Configuration parameters
 * -------------------------
 * 
 * Algorithm configuration is passed on construction and never changed.
 * The configuration data, `ConfigParams_t`, includes:
 * 
 * * `minHits` (default: `0`): minimum number of hits for a cluster to be
 *   considered
 * * `minChannelSpan` (default: `0`): do not consider clusters covering
 *   fewer than this number of channels (smallest to largest hit channel ID)
 * * `minTimeSpan` (default: `0`): do not consider clusters covering
 *   fewer than this number of TPC ticks (smallest/earliest to largest/latest
 *   TPC hit TDC)
 * * `maxHitClusterDistance` (default: `0.0`): hits further than this much [cm]
 *   from the cluster axis are ignored
 * * `minCloseHits` (default: `2`): after purging hits that are too far, keep
 *   only clusters with at least this many hits
 * 
 * 
 */
class icarus::AnodeToCathodeClusterAlg {
  
    public:
  
  /// Algorithm configuration parameters.
  struct ConfigParams_t {
    
    /// Number of the wire planes to consider (increasing).
    std::vector<geo::PlaneID::PlaneID_t> acceptedPlanes;
    
    unsigned int minHits = 0U; ///< Minimum number of hits in the cluster.
    unsigned int minChannelSpan = 0U; ///< Minimum range of channels to cover.
    float minTimeSpan = 0.0f; ///< Minimum range of TDC ticks to cover.
    
    /// Maximum distance between a hit and its cluster axis [cm]
    double maxHitClusterDistance = 0.0;
    
    /// Minimum span of cluster in drift direction [cm]
    double minDriftSpan = 0.0;
    
    /// Maximum span of cluster in drift direction [cm]
    double maxDriftSpan = std::numeric_limits<double>::max();
    
    /// Minimum number of hits close to cluster axis.
    unsigned int minCloseHits = 2U;
    
    /// Log category for algorithm messages.
    std::string logCategory = "AnodeToCathodeClusterAlg";
    
  }; // ConfigParams_t
  
  /// Parameters that may need to be updated.
  struct SetupParams_t {
    
    double wirePitch  = 0.0;        ///< [cm]
    double tickDuration = -1.0;     ///< [&micro;s]
    double driftVelocity = 0.0;     ///< [cm/&micro;s]
    double tickDriftDistance = 0.0; ///< [cm]
    
    SetupParams_t() = default; // needed, but bad
    SetupParams_t(double wirePitch, double tickDuration, double driftVelocity);
    
    bool isSet() const { return tickDuration > 0.0; }
    
  }; // SetupParams_t
  
  
  /// Rejection reason (see `explain()` for description).
  enum class reason {
      Accepted = 0     ///< Cluster is anode-to-cathode.
    , Invalid          ///< Cluster object is invalid.
    , NoPlaneInfo      ///< Cluster object is not associated to a readout plane.
    , PlaneIgnored     ///< Cluster object is on a vetoed wire plane.
    , TooFewHits       ///< Not enough hits in the cluster.
    , TooNarrowOnWires ///< Hits cross too few wires/channels.
    , TooNarrowOnTime  ///< Hits cross too few TDC ticks.
    , LinearFitFailure ///< Time/wire distribution not compatible with a line.
    , TooFewCloseHits  ///< Too few hits close to cluster axis.
    , TooShortDrift    ///< Too short span in drift direction
    , TooLongDrift     ///< Too long span in drift direction
    , NReasons         ///< Number of suppo
  }; // reason
  
  
  AnodeToCathodeClusterAlg(ConfigParams_t config)
    : fParams(std::move(config))
    {}
  
  AnodeToCathodeClusterAlg(ConfigParams_t config, SetupParams_t setup)
    : fParams(std::move(config)), fSetup(setup)
    {}
  
  
  /// Sets up the runtime parameters.
  void setup(SetupParams_t setup) { fSetup = std::move(setup); }
  
  /// Returns the current configuration parameters.
  ConfigParams_t const& config() const { return fParams; }
  
  /// Returns the runtime (set up) parameters.
  SetupParams_t const& setupParams() { return fSetup; }
  
  
  // --- BEGIN -- Operations ---------------------------------------------------
  /// @name Operations
  /// @{
  /**
   * @brief Returns whether the specified cluster qualifies as anode-to-cathode
   * @tparam Cluster type of pointer to `recob::Cluster`
   * @tparam Hits type of collection of pointers to `recob::Hit`
   * @param cluster the cluster object to be tested
   * @param hits the list of hits associated to the `cluster` being tested
   * @return `reason::accepted` if anode-to-cathode, or the `reason` why not;
   *         if accepted, a list of accepted hits is returned too
   * 
   * The `cluster` object must present the interface of a pointer to
   * `recob::Cluster` (e.g. `recob::Cluster const*`, `art::Ptr<recob::Cluster>`,
   * `std::unique_ptr<recob::Cluster>`,
   * `std::vector<recob::Cluster>::const_iterator`...).
   * The `hits` object must be a sequence of elements behaving like pointers to
   * `recob::Hit` objects. It is recommended that `art::Ptr<recob::Hit>` be used
   * if the result is used to create a new cluster or to match the input hits.
   * 
   * The return value is a pair. The first element is the `reason` of the
   * rejection. If that reason is `Accepted`, `TooShortDrift` or `TooLongDrift`,
   * the list of hits in the second element contains the hits that have been
   * considered for this cluster (e.g. excluding the ones too far from the
   * cluster axis).
   * In addition, if the reason is `Accepted`, the cluster is not rejected at
   * all.
   * 
   */
  template <typename Cluster, typename Hits>
  std::pair<reason, std::vector<typename Hits::value_type>> accept
    (Cluster const& cluster, Hits const& hits) const;
  
  /// Alias for `accept()`.
  template <typename Cluster, typename Hits>
  auto operator() (Cluster const& cluster, Hits const& hits) const
    { return accept(cluster, hits); }
  
  /// @}
  // --- END -- Operations -----------------------------------------------------
  
  /// Returns a C-string explaining rejection reason `r`.
  static const char* explain(reason r);
  
  
  /**
   * @brief Returns the span of full channel range covered by the non-null hits.
   * @tparam Hits type of collection of pointers to hits
   * @tparam hits collection of pointers to hits
   * @return minimum and maximum peak time among all the `hits`
   * 
   * Hit pointers which are null are skipped.
   * The same value of minimum and maximum are returned if the number of valid
   * hits is less than two.
   * 
   * The `Hits` collection must be forward-traversable and each element must
   * have a `recob::Hit const*` interface (member access `operator->` and
   * conversion to `bool`).
   */
  /**
   * Returns the span of the full range of time covered by the non-null hits.
   */
  template <typename Hits>
  static std::pair<float, float> hitTimeSpan(Hits const& hits);
  
  /**
   * @brief Returns the span of full channel range covered by the non-null hits.
   * @tparam Hits type of collection of pointers to hits
   * @tparam hits collection of pointers to hits
   * @return minimum and after-maximum channel number among all the `hits`
   *
   * Hits are assumed to be sorted by increasing channel number.
   * So this function is not that hard after all.
   * The second returned element is the largest ID among the hits, increased by
   * `1`.
   * 
   * The `Hits` collection must be forward-traversable and each element must
   * have a `recob::Hit const*` interface (member access `operator->` and
   * conversion to `bool`).
   */
  template <typename Hits>
  static std::pair<raw::ChannelID_t, raw::ChannelID_t> hitChannelSpan
    (Hits const& hits);
  
  
  /// Does `A` hit peak earlier than `B`?
  template <typename HitPtrA, typename HitPtrB>
  static bool earlierHitPeak(HitPtrA const& A, HitPtrB const& B)
    { return A->PeakTime() < B->PeakTime(); }
  
  
    private:
  
  /// -- BEGIN -- Data members -------------------------------------------------
  
  ConfigParams_t fParams; ///< Algorithm configuration parameters.
  
  SetupParams_t fSetup; ///< Setup parameters;
  
  /// -- END -- Data members ---------------------------------------------------
  
  
  /// Performs the analysis of the cluster (see `accept()`).
  std::pair<reason, std::vector<recob::Hit const*>> clusterAnalysis
    (recob::Cluster const& cluster, std::vector<recob::Hit const*> hits) const;

  /**
   * @brief Selects hits not too far from the line of the cluster.
   * @param hits (input/output) collection of hits to be considered
   * @return the number of non-outlier hits, and slope and intercept of the line
   * 
   * The algorithm converts the hit time and wire into points with homogeneous
   * coordinates, and fits them with linear interpolation.
   * At each iteration, the resulting line is compared with the points that it
   * was extracted from. If no point is farther than the configured maximum
   * distance the fit is accepted. Otherwise, the farthest point is discarded
   * and a new fit is considered with all previous points except this one.
   * The procedure ultimately fails when there are only 2 points left.
   * 
   * The allowed distance is set by the `maxHitClusterDistance` parameter.
   * 
   * The returned values include the number of hits in the final fit, its
   * intercept and slope. Also, `hits` collection is updated, with the pointers
   * of all hits that were found too far set to `nullptr`.
   */
  std::tuple<unsigned int, double, double> removeOutlyingHits
    (std::vector<recob::Hit const*>& hits) const;
  
  
  /// Returns a `HitCoordinateBuilder` set up with the current parameters.
  icarus::HitCoordinateBuilder makeHitCoordinates() const;
  
  /// Transforms a hit collection into a collection of pointers to `recob::Hit`.
  template <typename Hits>
  static std::vector<recob::Hit const*> makeHitVector(Hits const& hits);
  
  /// Returns whether `plane` is a plane number configured for processing.
  bool isAcceptedPlane(geo::PlaneID const& plane) const;
  
  
}; // class icarus::AnodeToCathodeClusterAlg


namespace icarus {
  
  /// Prints the specified configuration into a output stream.
  /// 
  /// Does not start a new line; it ends with a new line.
  std::ostream& operator<< (
    std::ostream& out,
    icarus::AnodeToCathodeClusterAlg::ConfigParams_t const& config
    );
  
} // namespace icarus


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <typename Cluster, typename Hits>
auto icarus::AnodeToCathodeClusterAlg::accept
  (Cluster const& cluster, Hits const& hits) const
  -> std::pair<reason, std::vector<typename Hits::value_type>>
{
  /*
   * We have a bit of metaprogramming here.
   * 
   * In general, we expect a collection of pointers to `recob::Hit` as input,
   * but we don't know which type of pointers.
   * To keep it simple, the algorithm only accepts C pointers
   * (`recob::Hit const*`). So we translate the input hit list for the algorithm
   * and the algorithm returns a list of the selected hits, also simple
   * pointers.
   * But we have promised to return hit pointers of the same type as our input.
   * The actual pointers to recob::Hit objects are still the same, but the
   * pointers and their order may be different. We need to map back from the
   * selected hit list to the input list, and we do that using the pointer to
   * the actual hit as a connecting key.
   * This is a bit cumbersome but it's reasonably fine.
   * 
   * As a bonus, though, if the input is already `recob::Hit const*` type of
   * pointers (`HitCpointers` below), we skip the mapping part completely
   * and directly use the hit list returned by the algorithm.
   * 
   */
  using SrcHit_t = typename Hits::value_type;
  static constexpr bool HitCpointers
    = std::is_same_v<std::decay_t<std::remove_pointer_t<SrcHit_t>>, recob::Hit>;
  
  if (!cluster) return { reason::Invalid, {} };
  
  std::vector<recob::Hit const*> myHits;
  std::unordered_map<recob::Hit const*, SrcHit_t> hitMap;
  
  if constexpr (HitCpointers) myHits = hits; // copy
  else {
    for (auto const& hit: hits) {
      recob::Hit const* hitAddress = &*hit;
      myHits.push_back(hitAddress);
      hitMap[hitAddress] = hit;
    } // for
  }
  
  auto [ reason, selectedHits ] = clusterAnalysis(*cluster, std::move(myHits));
  
  // remap the result into the input hits
  std::vector<SrcHit_t> mappedSelectedHits;
  for (recob::Hit const* hit: selectedHits) {
    if (!static_cast<bool>(hit)) continue;
    if constexpr (HitCpointers) mappedSelectedHits.push_back(hit);
    else mappedSelectedHits.push_back(hitMap.at(hit));
  } // for
  
  
  return { reason, mappedSelectedHits };
} // icarus::AnodeToCathodeClusterAlg::accept()


// -----------------------------------------------------------------------------
template <typename Hits>
std::pair<float, float> icarus::AnodeToCathodeClusterAlg::hitTimeSpan
  (Hits const& hits)
{
  if (empty(hits)) return { 0.0f, 0.0f };
  
  float min = std::numeric_limits<float>::max();
  float max = std::numeric_limits<float>::lowest();
  for (auto const& hit: hits) {
    if (!static_cast<bool>(hit)) continue;
    float const time = hit->PeakTime();
    if (min > time) min = time;
    if (max < time) max = time;
  } // for
  
  return { min, max };
} // icarus::AnodeToCathodeClusterAlg::hitTimeSpan()


// -----------------------------------------------------------------------------
template <typename Hits>
std::pair<raw::ChannelID_t, raw::ChannelID_t>
icarus::AnodeToCathodeClusterAlg::hitChannelSpan(Hits const& hits) {
  // assumed sorted
  static auto nonNull = [](auto* ptr){ return static_cast<bool>(ptr); };
  auto const first = std::find_if(hits.begin(), hits.end(), nonNull);
  if (first == hits.end()) return { 0, 0 }; // meh
  auto const last = std::find_if(hits.rbegin(), hits.rend(), nonNull);
  return { (*first)->Channel(), (*last)->Channel() + 1 };
} // icarus::AnodeToCathodeClusterAlg::hitChannelSpan()


// -----------------------------------------------------------------------------
template <typename Hits>
std::vector<recob::Hit const*> makeHitVector(Hits const& hits) {
  std::vector<recob::Hit const*> myHits;
  for (auto const& hit: hits) if (hit) myHits.push_back(&*hit);
  return myHits;
} // icarus::AnodeToCathodeClusterAlg::makeHitVector(Hits)


// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
icarus::AnodeToCathodeClusterAlg::SetupParams_t::SetupParams_t
  (double wirePitch, double tickDuration, double driftVelocity)
  : wirePitch(wirePitch), tickDuration(tickDuration)
  , driftVelocity(driftVelocity)
  , tickDriftDistance(driftVelocity*tickDuration)
  {}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_ANALYSIS_ALGORITHMS_ANODETOCATHODECLUSTERALG_H
