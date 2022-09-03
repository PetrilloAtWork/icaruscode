/**
 * @file   icaruscode/Utilities/TrackShifterAlg.cxx
 * @brief  Creates a direct association between tracks and time.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 31, 2022
 * @see    icaruscode/Utilities/TrackShifterAlg.h
 *
 */


// ICARUS libraries
#include "icaruscode/Utilities/TrackShifterAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <utility> // std::move()
#include <limits>
#include <cassert>


//------------------------------------------------------------------------------
icarus::TrackShifterAlg::TrackShifterAlg
  (std::string const& logCategory /* = "TrackShifterAlg" */)
  : icarus::ns::util::mfLoggingClass{ logCategory }
  {}


//------------------------------------------------------------------------------
void icarus::TrackShifterAlg::setup(
  geo::GeometryCore const& geom,
  detinfo::DetectorPropertiesData detProp
  )
{
  fGeom = &geom;
  fDetProp.emplace(std::move(detProp));
} // icarus::TrackShifterAlg::setup()


//------------------------------------------------------------------------------
void icarus::TrackShifterAlg::unsetup() {
  // mostly symbolic
  fGeom = nullptr;
  fDetProp.reset();
} // icarus::TrackShifterAlg::unsetup()


//------------------------------------------------------------------------------
bool icarus::TrackShifterAlg::isSetup() const
  { return fGeom && fDetProp; }


//------------------------------------------------------------------------------
recob::Track icarus::TrackShifterAlg::shiftTrack
  (recob::Track const& track, nanoseconds shift) const
{
  /*
   * For each point, we apply a shift by -v(d) * shift.
   * One problem is to determine the drift velocity v(d) (which is a vector).
   * 
   * So for each point in the trajectory we find which TPC the point is closest
   * to, and use its geometric drift direction as v(d).
   * If a point is invalid, though, detection of the TPC may be unreliable:
   * we use instead the drift direction from the previous point, no matter what
   * (at the very beginning, instead, that direction is set to the one for the
   * first valid point).
   * 
   * We limit all searches to the cryostat the track belongs (that is the one
   * closest to the median track valid point -- more or less).
   */
  using util::quantities::intervals::microseconds;
  using namespace util::quantities::time_literals;
  
  if (!isSetup())
    throw cet::exception{ "TrackShifterAlg" } << "Algorithm not set up!\n";
  
  if (shift == 0_ns) return track; // no shift, no party
  
  recob::TrackTrajectory const& traj = track.Trajectory();
  
  recob::TrackTrajectory::Positions_t shiftedPositions
    { traj.Trajectory().Positions() };
  
  double const driftVelocity = fDetProp->DriftVelocity(); // cm/us
  double const shiftSize = driftVelocity * microseconds{ shift }.value();
  
  geo::CryostatGeo const& cryo = *(trackCryostat(track.Trajectory()));
  assert((track.CountValidPoints() == 0) || (&cryo != nullptr));
  
  geo::TPCGeo const* lastTPC = nullptr;
  geo::Vector_t driftDirection;
  if (std::size_t const iFirst = track.FirstValidPoint();
      iFirst != recob::TrackTrajectory::InvalidIndex
  ) {
    lastTPC = findNearestTPC(track.LocationAtPoint(iFirst), cryo);
    assert(lastTPC);
    driftDirection = lastTPC->DriftDir<geo::Vector_t>();
    mfLogTrace() << "In " << lastTPC->ID() << ": shifting by "
      << (shiftSize*driftDirection) << " cm";
  }
  
  std::size_t iPoint = 0;
  for (recob::tracking::Point_t& pos: shiftedPositions) {
    
    // we don't touch invalid points
    if (!track.HasValidPoint(iPoint++)) continue;
    
    assert(lastTPC);
    
    // if `pos` is not inside the last TPC, we look for a new one
    if (!lastTPC->ContainsPosition(pos)) {
      lastTPC = findNearestTPC(pos, cryo); // may still be the same as before
      driftDirection = lastTPC->DriftDir<geo::Vector_t>();
      mfLogTrace() << "Point [" << (iPoint-1) << "] now in "
        << lastTPC->ID() << ": shifting by " << (shiftSize*driftDirection)
        << " cm";
    }
    
    // FIXME: this is a major issue: we may be moving points across cathode,
    //        and shifting by +t and then -t makes a mess instead of reverting
    
    recob::tracking::Vector_t const shift = driftDirection * shiftSize;
    
    pos += shift;
  } // for points
  
  // everything is a copy of the original track, except for the positions
  auto copyOf = [](auto v){ return v; };
  return {
      std::move(shiftedPositions)              // positions
    , copyOf(traj.Trajectory().Momenta())      // momenta
    , copyOf(traj.Flags())                     // flags
    , track.HasMomentum()                      // hasMomenta
    , track.ParticleId()                       // PId
    , track.Chi2()                             // Chi2
    , track.Ndof()                             // Ndof
    , copyOf(track.VertexCovarianceLocal5D())  // CovVertex
    , copyOf(track.EndCovarianceLocal5D())     // CovEnd
    , track.ID()                               // tkID
    };
  
} // DriftTracks::shiftTrack()


//------------------------------------------------------------------------------
geo::CryostatGeo const* icarus::TrackShifterAlg::trackCryostat
  (recob::TrackTrajectory const& track) const
{
  constexpr auto InvalidIndex = recob::TrackTrajectory::InvalidIndex;
  
  // pick the valid point "closest" to the median point of the track;
  // a track is required to have at least two points.
  std::size_t const iFirst = track.FirstValidPoint();
  if (iFirst == InvalidIndex) return nullptr;
  
  std::size_t const iLast = track.LastValidPoint();
  assert(iLast != InvalidIndex);
  
  std::size_t const iMiddle = (iFirst + iLast) / 2;
  
  std::size_t iChosen = iMiddle;
  if (!track.HasValidPoint(iMiddle)) {
    
    std::size_t const iAfter = track.NextValidPoint(iMiddle);
    assert((iAfter == InvalidIndex) || (iAfter > iMiddle));
  
    std::size_t const iBefore = track.PreviousValidPoint(iMiddle);
    assert((iBefore == InvalidIndex) || (iBefore < iMiddle));
    
    if (iAfter == InvalidIndex) iChosen = iBefore;
    else if (iBefore == InvalidIndex) iChosen = iAfter;
    else {
      iChosen = ((iMiddle - iBefore) <= (iAfter - iMiddle))? iBefore: iAfter;
    }
  }
  assert(iChosen != InvalidIndex);
  
  return findNearestCryostat(track.LocationAtPoint(iChosen));
} // icarus::TrackShifterAlg::trackCryostat()


//------------------------------------------------------------------------------
geo::CryostatGeo const* icarus::TrackShifterAlg::findNearestCryostat
  (geo::Point_t const& point) const
{
  return findNearestBox(point, fGeom->IterateCryostats());
} // icarus::TrackShifterAlg::findNearestCryostat()



//------------------------------------------------------------------------------
geo::TPCGeo const* icarus::TrackShifterAlg::findNearestTPC
  (geo::Point_t const& point, geo::CryostatGeo const& cryo) const
{
  return findNearestBox(point, cryo.IterateTPCs());
} // icarus::TrackShifterAlg::findNearestTPC()


//------------------------------------------------------------------------------
template <typename Boxes>
typename Boxes::value_type const* icarus::TrackShifterAlg::findNearestBox
  (geo::Point_t const& point, Boxes const& boxes)
{
  using std::empty;
  using Box_t = typename Boxes::value_type;
  
  if (empty(boxes)) return nullptr;
  
  double minDistance2 = std::numeric_limits<double>::max();
  Box_t const* closestBox = nullptr;
  for (Box_t const& box: boxes) {
    double const d2 = (box.Center() - point).Mag2();
    if (closestBox && (minDistance2 <= d2)) continue;
    minDistance2 = d2;
    closestBox = &box;
  } // for all boxes
  
  return closestBox;
} // icarus::TrackShifterAlg::findNearestBox()


//------------------------------------------------------------------------------

