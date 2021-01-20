/**
 * @file   icaruscode/TPC/SignalProcessing/HitFinder/HitCoordinateBuilder.cxx
 * @brief  Simple utility for converting a hit into 2D homogeneous coordinates.
 * @author Gianluca Petrillo
 * @date   January 15, 2021
 * @see    `icaruscode/TPC/SignalProcessing/HitFinder/HitCoordinateBuilder.h`
 * 
 * This library is header-only.
 */


// library header
#include "icaruscode/TPC/SignalProcessing/HitFinder/HitCoordinateBuilder.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"


// -----------------------------------------------------------------------------
icarus::HitCoordinateBuilder icarus::makeHitCoordinates(
  detinfo::DetectorClocksData const& detClocks,
  detinfo::DetectorPropertiesData const& detProp,
  geo::PlaneGeo const& plane
) {
  return {
      detProp.DriftVelocity()           // driftVelocity
    , detClocks.TPCClock().TickPeriod() // tickDuration 
    , plane.WirePitch()                 // wirePitch
    };
} // icarus::makeHitCoordinates(PlaneGeo)


// -----------------------------------------------------------------------------
icarus::HitCoordinateBuilder icarus::makeHitCoordinates(
  geo::GeometryCore const& geom,
  detinfo::DetectorClocksData const& detClocks,
  detinfo::DetectorPropertiesData const& detProp,
  geo::PlaneID const& plane
) {
  return makeHitCoordinates(detClocks, detProp, geom.Plane(plane));
} // icarus::makeHitCoordinates(GeometryCore)


// -----------------------------------------------------------------------------
