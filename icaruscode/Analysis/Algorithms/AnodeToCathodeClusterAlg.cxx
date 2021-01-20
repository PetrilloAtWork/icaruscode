/**
 * @file   icaruscode/Analysis/Algorithms/AnodeToCathodeClusterAlg.cxx
 * @brief  Algorithm: selects clusters spanning a full anode/cathode range.
 * @author Christian Farnese, Weskey Ketchum, Olivia Bitter, Gianluca Petrillo
 * @date   January 12, 2021
 * @see    icaruscode/Analysis/Algorithms/AnodeToCathodeClusterAlg.h
 * 
 * The algorithm of this module is taken from LAr purity analysis code
 * (`ICARUSPurityDQMnew`).
 * 
 */

// library header
#include "icaruscode/Analysis/Algorithms/AnodeToCathodeClusterAlg.h"
#include "icaruscode/TPC/SignalProcessing/HitFinder/HitCoordinateBuilder.h"

// LArSoft libraries
#include "lardata/Utilities/SimpleFits.h" // lar::util::LinearFit
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"

// framework libraries
#include "cetlib/pow.h" // cet::square()
#include "cetlib_except/exception.h"

// C/C++ standard library
#include <ostream>
#include <algorithm> // std::sort(), std::find_if(), ...
#include <cassert>


// -----------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  //
  // Hit comparison functions
  //
  
  /// Does `A` peak earlier than `B`?
  bool earlier(recob::Hit const* A, recob::Hit const* B)
    { return A->PeakTime() < B->PeakTime(); }
  
  /// Has `A` lower channel than `B`? if the same, does `A` peak earlier?
  bool lowerChannelThenTime(recob::Hit const* A, recob::Hit const* B)
    {
      if (A->Channel() < B->Channel()) return true;
      if (A->Channel() > B->Channel()) return false;
      return earlier(A, B);
    } // lowerChannelThenTime()
  
  // ---------------------------------------------------------------------------
  /// Returns whether value `key` is in an element between `b` and `e`.
  template <typename BIter, typename EIter, typename T>
  bool contains(BIter b, EIter e, T const& key)
    { return std::find(b, e, key) != e; }
  
  /// Returns whether value `key` is in an element of container `cont`.
  template <typename Cont, typename T>
  bool contains(Cont const& cont, T const& key)
    { return contains(begin(cont), end(cont), key); }
  
  
  // ---------------------------------------------------------------------------
  template <typename T>
  constexpr bool isMax(T value)
    { return value == std::numeric_limits<T>::max(); }
  
  
  // ---------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
auto icarus::AnodeToCathodeClusterAlg::clusterAnalysis
  (recob::Cluster const& cluster, std::vector<recob::Hit const*> hits) const
  -> std::pair<reason, std::vector<recob::Hit const*>>
{
  static auto const failure = [](reason r)
    { return std::pair<reason, std::vector<recob::Hit const*>>{ r, {} }; };
  
  if (!fSetup.isSet()) {
    throw cet::exception(fParams.logCategory)
      << "AnodeToCathodeClusterAlg algorithm executed before setup().";
  }
  
  if (!cluster.isValid()) return failure(reason::Invalid);
  if (!cluster.hasPlane()) return failure(reason::NoPlaneInfo); // could work around this...
  
  geo::PlaneID const plane = cluster.Plane();
  if (!isAcceptedPlane(plane)) return failure(reason::PlaneIgnored);
  
  
  /*
   * 6. selection:
   *    * time span > 1100, wire span > 100
   * 7. analysis
   *   1. at least 30 hits (wasn't it already required? TODO check with assert)
   *   2. linear fit of time vs. wire
   *      * one by one, progressively exclude hits that lie too far from the fit
   *        (delta^2/(slope^2+1) > 3, inhomogeneous)
   *      * stick to a list of hits without such outliers
   *   3. preserve hits with distance from fit of less than 3 mm
   *   4. keep only clusters with at least 100 good hits
   *   5. find earliest and latest hit
   */
  
  if (size(hits) < std::max(fParams.minHits, 1U))
    return failure(reason::TooFewHits);
  
  // sort the hits by increasing channel, then by increasing peak time
  std::sort(hits.begin(), hits.end(), &lowerChannelThenTime);
  
  //
  // check span in time and channel
  //
  auto const channelRange = hitChannelSpan(hits);
  if ((channelRange.second - channelRange.first) < fParams.minChannelSpan)
    return failure(reason::TooNarrowOnWires);
  
  auto const timeRange = hitTimeSpan(hits);
  if ((timeRange.second - timeRange.first) < fParams.minTimeSpan)
    return failure(reason::TooNarrowOnTime);
  
  //
  // hit selection; excluded hits are going to be set to nullptr
  //
  auto const [ nSelectedHits, slope, intcp ] = removeOutlyingHits(hits);
  if (nSelectedHits != hits.size()) {
    mf::LogTrace(fParams.logCategory)
      << "Cluster " << cluster.ID() << " lost " << (hits.size() - nSelectedHits)
      << " / " << hits.size() << " hits";
  }
  if (nSelectedHits <= 2U)
    return failure(reason::LinearFitFailure); // hard to imagine...
  
  //
  // minimum hits, extension
  //
  if (nSelectedHits < fParams.minCloseHits)
    return failure(reason::TooFewCloseHits);
  
  auto const cleanChannelRange [[maybe_unused]] = hitChannelSpan(hits);
  auto const cleanTimeRange = hitTimeSpan(hits);
  float const cleanTimeSpan = cleanTimeRange.second - cleanTimeRange.first;
  float const cleanDriftSpan = cleanTimeSpan * fSetup.tickDriftDistance;
  MF_LOG_TRACE(fParams.logCategory)
    << "Cluster " << cluster.ID() << " spans "
    << (cleanChannelRange.second - cleanChannelRange.first)
    << " channels and " << cleanTimeSpan << " TPC ticks (" << cleanDriftSpan
    << " cm).";
  
  if (cleanDriftSpan < fParams.minDriftSpan)
    return { reason::TooShortDrift, std::move(hits) };
  if (cleanDriftSpan > fParams.maxDriftSpan)
    return { reason::TooLongDrift, std::move(hits) };
  
  //
  // then it's ok
  //
  
  return { reason::Accepted, std::move(hits) };
  
} // icarus::AnodeToCathodeClusterAlg::clusterAnalysis()


// -----------------------------------------------------------------------------
std::tuple<unsigned int, double, double>
icarus::AnodeToCathodeClusterAlg::removeOutlyingHits
  (std::vector<recob::Hit const*>& hits) const
{
  
  unsigned int nSelectedHits = hits.size();
  
  // we assume all hits are on the same wire plane
  auto const getHitCoords = makeHitCoordinates();
  
  using HitCoords_t = icarus::HitCoordinateBuilder::HitCoords_t;
  
  std::vector<HitCoords_t> hitCoords;
  hitCoords.reserve(nSelectedHits);
  for (auto const& hit: hits)
    hitCoords.push_back(getHitCoords.coordinatesFromChannel(hit));
  
  lar::util::LinearFit<double> fitter;
  while (nSelectedHits > 2U) {
    
    // this is a waste, the linear fitter should support removal of points...
    fitter.clear();
    for (auto const& [ hit, coords ]: util::zip(hits, hitCoords)) {
      if (!hit) continue;
      fitter.add(coords.wire, coords.drift);
    } // for
    
    // find the farthest hit; unless it's close, and then it's ok;
    // the residual is the (square of) 2D distance of hit from fit
    // (computed as difference of drift at hit time, projected normal to fit)
    auto const getResidual2 = [
      intcp=fitter.Intercept(),
      slope=fitter.Slope(),
      projAngle2=1.0 / cet::sum_of_squares(fitter.Slope(), 1.0)
      ](double channel, double time){
        return cet::square((intcp + slope * channel) - time) * projAngle2;
      };
    
    MF_LOG_TRACE(fParams.logCategory)
      << "Cluster fit #" << (hits.size() - nSelectedHits + 1) << " on "
      << hits.size() << " hits: "
      << fitter.Intercept() << " + channel x " << fitter.Slope()
      ;
    
    double max_d2 = cet::square(fParams.maxHitClusterDistance);
    std::size_t iFarthestHit = std::numeric_limits<std::size_t>::max();
    for (auto const& [ iHit, hit, coords ]: util::enumerate(hits, hitCoords)) {
      if (!hit) continue;
      
      double const residual2 = getResidual2(coords.wire, coords.drift);
      /*
      mf::LogTrace log(fParams.logCategory);
      log 
        << "  hit #" << iHit
        << " channel " << hit->Channel() << " peak " << hit->PeakTime()
        << " at ( " << coords.wire << ", " << coords.drift << " )"
        << " => " << residual2
        ;
      */
      if (residual2 <= max_d2) continue;
      // log << " (too much and worst so far)";
      max_d2 = residual2;
      iFarthestHit = iHit;
      
    } // for
    
    // if the maximum distance does not exceed the threshold, we are good to go
    if (max_d2 == cet::square(fParams.maxHitClusterDistance)) break;
    
    // mask the hit and try again;
    // if the fitter supported it, we would just remove the hit
    MF_LOG_TRACE(fParams.logCategory)
      << "  => discarding hit #" << iFarthestHit
      << " channel " << hits[iFarthestHit]->Channel()
      << " peak " << hits[iFarthestHit]->PeakTime()
      ;
    
    assert(iFarthestHit < hits.size());
    hits[iFarthestHit] = nullptr;
    --nSelectedHits;
    
  } // while
  
  if (nSelectedHits <= 2U) return { 0U, 0.0, 0.0 }; // failure
  
  return { nSelectedHits, fitter.Intercept(), fitter.Slope() };
  
} // icarus::AnodeToCathodeClusterAlg::removeOutlyingHits()


// -----------------------------------------------------------------------------
icarus::HitCoordinateBuilder
icarus::AnodeToCathodeClusterAlg::makeHitCoordinates() const {
  assert(fSetup.isSet());
  return icarus::HitCoordinateBuilder
    { fSetup.driftVelocity, fSetup.tickDuration, fSetup.wirePitch };
} // icarus::AnodeToCathodeClusterAlg::makeHitCoordinates()


// -----------------------------------------------------------------------------
const char* icarus::AnodeToCathodeClusterAlg::explain(reason r) {
  
  switch (r) {
    
    case reason::Accepted        : return "cluster is anode-to-cathode";
    case reason::Invalid         : return "cluster object is invalid";
    case reason::NoPlaneInfo     : return "cluster object is not associated to a readout plane";
    case reason::PlaneIgnored    : return "cluster object is on a vetoed wire plane";
    case reason::TooFewHits      : return "not enough hits in the cluster";
    case reason::TooNarrowOnWires: return "hits cross too few wires/channels";
    case reason::TooNarrowOnTime : return "hits cross too few TDC ticks";
    case reason::LinearFitFailure: return "time/wire distribution not compatible with a line";
    case reason::TooFewCloseHits : return "too few hits close to cluster axis";
    case reason::TooShortDrift   : return "too short span in drift direction";
    case reason::TooLongDrift    : return "too long span in drift direction";
    case reason::NReasons        : return "number of supported reasons";
    
  } // switch
  
  return "unknown reason";
  
} // icarus::AnodeToCathodeClusterAlg::explain()


// -----------------------------------------------------------------------------
bool icarus::AnodeToCathodeClusterAlg::isAcceptedPlane
  (geo::PlaneID const& plane) const
{
  return fParams.acceptedPlanes.empty()
    || contains(fParams.acceptedPlanes, plane.Plane);
} // icarus::AnodeToCathodeClusterAlg::isAcceptedPlane()


// -----------------------------------------------------------------------------
std::ostream& icarus::operator<< (
  std::ostream& out,
  icarus::AnodeToCathodeClusterAlg::ConfigParams_t const& config
) {
  
  out <<   " * minimum number of hits in the original cluster:     ";
  if (config.minHits > 0) out << config.minHits;
  else                    out << "none";
  
  out << "\n * minimum range of channels in the original cluster:  ";
  if (config.minChannelSpan > 0) out << config.minChannelSpan << " channels";
  else                           out << "none";
  
  out << "\n * minimum range of TPC ticks in the original cluster: ";
  if (config.minTimeSpan > 0) out << config.minTimeSpan << " ticks";
  else                        out << "none";
  
  out << "\n * maximum distance of good hits from cluster axis:    ";
  if (config.maxHitClusterDistance > 0.0)
    out << config.maxHitClusterDistance << " cm";
  else out << "none";
  
  out << "\n * minumum number of good (close) hits:                ";
  if (config.minCloseHits > 0.0) out << config.minCloseHits;
  else                           out << "none";
  
  out << "\n * drift distance covered by the cluster:              ";
  if (config.minDriftSpan != 0.0) {
    if (isMax(config.maxDriftSpan))
      out << "> " << config.minDriftSpan << " cm";
    else
      out << config.minDriftSpan << " -- " << config.maxDriftSpan << " cm";
  }
  else {
    if (isMax(config.maxDriftSpan))
      out << "any";
    else
      out << ">" << config.maxDriftSpan << " cm";
  }
  
  if (!config.acceptedPlanes.empty()) {
    auto iPlane = config.acceptedPlanes.begin();
    out << "\n * includes only clusters from " << config.acceptedPlanes.size()
      << " planes:              #" << *iPlane;
    while (++iPlane != config.acceptedPlanes.end()) out << ", #" << *iPlane;
  } // if planes selected
  
  return out << "\n";
} // operator<< (icarus::AnodeToCathodeClusterAlg::ConfigParams_t)


// -----------------------------------------------------------------------------
