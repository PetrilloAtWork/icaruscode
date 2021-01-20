/**
 * @file   icaruscode/TPC/SignalProcessing/HitFinder/HitCoordinateBuilder.h
 * @brief  Simple utility for converting a hit into 2D homogeneous coordinates.
 * @author Gianluca Petrillo
 * @date   January 15, 2021
 * @see    `icaruscode/TPC/SignalProcessing/HitFinder/HitCoordinateBuilder.h`
 * 
 * The header is sufficient for the use of `icarus::HitCoordinateBuilder`,
 * while to use `icarus::makeHitCoordinates()` the library
 * `icaruscode_TPC_SignalProcessing_HitFinder` (nice name!) needs to be linked.
 */

#ifndef ICARUSCODE_TPC_SIGNALPROCESSING_HITFINDER_HITCOORDINATEBUILDER_H
#define ICARUSCODE_TPC_SIGNALPROCESSING_HITFINDER_HITCOORDINATEBUILDER_H


// LArSoft libraries
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::WireID

// C/C++ standard libraries
#include <utility> // std::get() (for convenience)


// -----------------------------------------------------------------------------
namespace geo { class GeometryCore; class PlaneGeo; }
namespace detinfo { class DetectorClocksData; class DetectorPropertiesData; }
namespace icarus {
  class HitCoordinateBuilder; 
  
  /**
   * @brief Creates and returns a `icarus::HitCoordinateBuilder` object.
   * @param geom geometry service provider
   * @param detClocks data from detector clocks provider
   * @param detProp data from detector properties provider
   * @param plane ID of the plane to make the hit coordinate builder for
   * @return an instance of `icarus::HitCoordinateBuilder` for hits on `plane`
   * 
   * This wrapper function extracts the information from the service providers
   * necessary to initialize a `icarus::HitCoordinateBuilder` for hits in the
   * specified plane.
   */
  HitCoordinateBuilder makeHitCoordinates(
    geo::GeometryCore const& geom,
    detinfo::DetectorClocksData const& detClocks,
    detinfo::DetectorPropertiesData const& detProp,
    geo::PlaneID const& plane
    );
  
  /**
   * @brief Creates and returns a `icarus::HitCoordinateBuilder` object.
   * @param detClocks data from detector clocks provider
   * @param detProp data from detector properties provider
   * @param plane plane hosting the hits the hit coordinate builder is for
   * @return an instance of `icarus::HitCoordinateBuilder` for hits on `plane`
   * 
   * This wrapper function extracts the information from the service providers
   * necessary to initialize a `icarus::HitCoordinateBuilder` for hits in the
   * specified plane.
   */
  HitCoordinateBuilder makeHitCoordinates(
    detinfo::DetectorClocksData const& detClocks,
    detinfo::DetectorPropertiesData const& detProp,
    geo::PlaneGeo const& plane
    );
  
} // namespace icarus

/**
 * @brief Returns a 2D homogeneous coordinate for a hit
 * @see `icarus::makeHitCoordinates()`
 * 
 * The coordinates represent wire position and hit peak time.
 * Conversion is performed via drift velocity and wire pitch, but the absolute
 * values of the coordinates, i.e. their reference point, is not defined.
 * 
 * Example of usage in LArSoft:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
 * auto const detTimings {
 *   detinfo::makeDetectorTimings
 *     (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
 *   };
 * detinfo::DetectorPropertiesData const detProp {
 *   art::ServiceHandle<detinfo::DetectorPropertiesService const>()
 *     ->DataForJob(detTimings.clockData())
 *   };
 * 
 * icarus::HitCoordinateBuilder const hitCoords {
 *   geom.WirePitch(plane)
 *   , detTimings.ClockPeriodFor<detinfo::timescale::TPCelectronics_time>().value()
 *   , detProp.DriftVelocity()
 *   };
 * 
 * auto [ wire, time ] = hitCoords(hit);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * where `hit` is a `recob::Hit` object
 * 
 * Different options are available:
 * * `coordinates()` explicitly requires the two coordinates, and it accepts a
 *     reference value; this function uses a wire pitch step on the first
 *     coordinate and a tick distance step on the second one;
 * * `coordinatesOf()` extracts the coordinates from the specified hit, using
 *     wire number for wire and peak time for tick;
 * * `coordinatesFromChannel()` extracts the coordinates from the specified hit,
 *     using channel number for wire and peak time for tick;
 * 
 * 
 */
class icarus::HitCoordinateBuilder {
  
  double fDrift;    ///< driftVelocity ionization charge drift velocity [cm/&micro;s]
  double fTickTime; ///< tickDuration time duration of the TDC tick [&micro;s]
  double fPitch;    ///< wirePitch distance between wires [cm]
  
    public:
  
  /// Hit homogeneous coordinates.
  struct HitCoords_t {
    double wire;  ///< Coordinate on wire plane.
    double drift; ///< Coordinate on drift direction.
  }; // HitCoords_t
  
  /**
   * @brief Sets this conversion object up.
   * @param driftVelocity ionization charge drift velocity [cm/&micro;s]
   * @param tickDuration time duration of the TDC tick [&micro;s]
   * @param wirePitch distance between wires [cm]
   * @see `icarus::makeHitCoordinates()`
   */
  HitCoordinateBuilder
    (double driftVelocity, double tickDuration, double wirePitch)
    : fDrift(driftVelocity), fTickTime(tickDuration), fPitch(wirePitch)
    {}
  
  
  //@{
  /// Converts the specified `wire` and `tick` coordinates.
  HitCoords_t coordinates(double wire, double tick) const;
  HitCoords_t coordinates(double wire, double tick, double refWire) const
    { return coordinates(wire - refWire, tick); }
  HitCoords_t coordinates
    (geo::WireID const& wire, double tick, double refWire = 0.0) const
    { return coordinates(static_cast<double>(wire.Wire), tick, refWire); }

  HitCoords_t operator() (double wire, double tick, double refWire = 0.0) const
    { return coordinates(wire, tick, refWire); }
  HitCoords_t operator()
    (geo::WireID const& wire, double tick, double refWire = 0.0) const
    { return coordinates(wire, tick, refWire); }
  //@}
  
  
  //@{
  /// Returns the wire and time coordinate of the specified `hit`.
  HitCoords_t coordinatesOf(recob::Hit const& hit, double refWire = 0.0) const
    { return coordinates(hit.WireID(), hit.PeakTime(), refWire); }
  HitCoords_t coordinatesOf(recob::Hit const* hit, double refWire = 0.0) const
    { return coordinatesOf(*hit, refWire); }
  
  HitCoords_t operator() (recob::Hit const& hit, double refWire = 0.0) const
    { return coordinatesOf(hit, refWire); }
  HitCoords_t operator() (recob::Hit const* hit, double refWire = 0.0) const
    { return coordinatesOf(hit, refWire); }
  //@}
  
  
  //@{
  /// Returns the channel and time coordinate of the specified `hit`.
  HitCoords_t coordinatesFromChannel
    (recob::Hit const& hit, raw::ChannelID_t refChannel = 0) const
    {
      return coordinates
        (static_cast<double>(hit.Channel() - refChannel), hit.PeakTime());
    }
  HitCoords_t coordinatesFromChannel
    (recob::Hit const* hit, raw::ChannelID_t refChannel = 0) const
    { return coordinatesFromChannel(*hit, refChannel); }
  //@}
  
}; // icarus::HitCoordinateBuilder


// -----------------------------------------------------------------------------
inline auto icarus::HitCoordinateBuilder::coordinates
  (double wire, double tick) const -> HitCoords_t
{
  return {
    wire * fPitch,             // wire
    tick * fTickTime * fDrift  // time
    };
} // icarus::HitCoordinateBuilder::coordinatesOf()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_TPC_SIGNALPROCESSING_HITFINDER_HITCOORDINATEBUILDER_H
