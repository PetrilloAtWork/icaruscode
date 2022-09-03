/**
 * @file   icaruscode/Utilities/TrackShifterAlg.h
 * @brief  Creates a direct association between tracks and time.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 31, 2022
 * @see    icaruscode/Utilities/TrackShifterAlg.cxx
 *
 */

#ifndef ICARUSCODE_UTILITIES_TRACKSHIFTERALG_H
#define ICARUSCODE_UTILITIES_TRACKSHIFTERALG_H


// ICARUS libraries
#include "icarusalg/Utilities/mfLoggingClass.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t...

// C/C++ standard libraries
#include <string>
#include <optional>
#include <utility> // std::forward()


//------------------------------------------------------------------------------
namespace geo {
  class GeometryCore;
  class CryostatGeo;
  class TPCGeo;
}

//------------------------------------------------------------------------------
namespace icarus { class TrackShifterAlg; }
/**
 * @brief Moves the track along the drift direction using time shift.
 * 
 * This algorithm applies a shift to the trajectory of the specified tracks.
 * The shift is a relative time shift, which is converted into a spatial shift
 * according to the location of the trajectory in the drift volumes.
 * 
 * The algorithm applies a shift independently to each point of the trajectory
 * of the track, as follows:
 * 
 *  * invalid points are not shifted at all but they are preserved;
 *  * points inside a TPC are shifted according to the nominal drift direction
 *    in that TPC;
 *  * points outside a TPC active volume are shifted according to the nominal
 *    drift direction of the closest TPC.
 * 
 * Note that this simple algorithm has effects on tracks crossing the cathode
 * that may be undesired, tearing apart or adding wrinkles to them (more below).
 * 
 * 
 * Determination of the drift amount
 * ----------------------------------
 * 
 * A track position in the drift direction is determined by its TPC charge
 * deposition time and the actual time the particle crosses the TPC. The former
 * is ultimately derived from the TPC digitized waveforms, in TPC readout clock
 * ticks from the start of the readout, usually via the peak time of a
 * reconstructed hit (`recob::Hit`). The track time is an "absolute" time
 * usually called @f$ t_{0} @f$. The position in the drift direction is obtained
 * using the drift velocity _v_, currently assumed constant, and an additional
 * assumption relating the track time and the hit time.
 * 
 * The common assumption in Pandora pattern recognition algorithms is that a
 * track crossing the anode plane exactly at trigger time is associated with a
 * known hit time set up in hardware. In ICARUS, this time is physically
 * determined by adding a delay of about 1.3 ms to the global trigger signal,
 * and using the delayed signal to stop the 4096-digit TPC readout window.
 * Given that the TPC readout clock has a period of 400 ns, this means that the
 * trigger time is 3250 ticks before the end of the buffer, i.e. at tick `846`
 * (all these values can be retrieved from the `detinfo::DetectorClocks` and
 * `detinfo::DetectorProperties` services). Pandora assigns a `anab::T0` time
 * value of `0` to such tracks.
 * 
 * This algorithm shifts the tracks according to the difference between the
 * reference track time, i.e. the time that was assumed when putting together
 * the input track, and a new time.
 *
 * The logic steps of the timing follow.
 * 
 * 1. Let's have a particle creating charge at position @f$ \vec{x} @f$ and
 *    @f$ t @f$, both absolute values
 * 2. The time the charge takes to get to the anode, "drift time", is
 *    @f$ \Delta t_{d} = -(\vec{x} - \vec{x}_{A})/\vec{v}_{d} @f$,
 *    with @f$ \vec{x}_{A} @f$ the anode position and @f$ \vec{v}_{d} @f$ the
 *    drift velocity. Here vector notation is used loosely: more properly, all
 *    vectors and positions should be projected to the same axis with the drift
 *    direction (in LArSoft SBN it's @f$ \hat{x} @f$).
 * 3. the time that charge arrives at the anode is
 *    @f$ t_{c} = t + \Delta t_{d} = t - (\vec{x} - \vec{x}_{A})/\vec{v}_{d} @f$.
 * 4. If we measure the particle time _relative_ to a time reference
 *    @f$ t_{R} @f$, the charge is measured on the anode
 *    @f$ \Delta t_{c} = t_{c} - t_{R} = (t - t_{R}) - (\vec{x} - \vec{x}_{A})/\vec{v}_{d} @f$
 *    after that time reference (typical reference is the trigger time).
 *    It's important to note that _charge from a particle crossing the cathode
 *    (at @f$ \vec{x} = \vec{x}_{A} @f$) exactly at the reference time
 *    @f$ t = t_{R} @f$ will be observed with @f$ \Delta t_{c} = 0 @f$.
 * 5. the number of TPC readout tick at which the charge is observed is
 *    @f$ T_{c} = T_{R} + \Delta t_{c}/\Delta t_{T} @f$,
 *    where @f$ \Delta t_{T} @f$ is ADC the sampling period (400 ns in ICARUS,
 *    500 ns in most other experiments) and @f$ T_{R} @f$ is a reference tick
 *    that is defined by the very same relation: when @f$ \Delta t_{c} = 0 @f$,
 *    @f$ T_{c} = T_{R} @f$, which spells out that when, as above,
 *    "a particle crossed the cathode
 *    (at @f$ \vec{x} = \vec{x}_{A} @f$) exactly at the reference time
 *    @f$ t = t_{R} @f$", we record its charge at the TPC readout reference
 *    time. Note that _the reference time we use when measuring the crossing
 *    time of the particle, @f$ t_{R} @f$, and the TPC readout reference tick
 *    number @f$ T_{R} @f$ are not the same concept nor value (although
 *    eventually they may be tied together by the hardware).
 * 6. Summarising, the reconstruction which is given the tick number
 *    @f$ T_{c} @f$ at which the charge from particle crossing at time @f$ t @f$
 *    was observed, will conclude that the particle was at a position
 *    @f$ \vec{x} = \vec{x}_{A} + \left[ (t - t_{R}) - (T_{c} - T_{R}) \Delta t_{T} \right] \vec{v}_{d} @f$
 *
 * This is the "master formula" to pin down the position of the track in the
 * drift direction.
 * 
 * Suppose now that we have a new measurement for the absolute particle time,
 * @f$ t^{*} @f$, which replaces @f$ t @f$. We have already a track shifted
 * to @f$ \vec{x} @f$ according to @f$ t @f$ using the prescription above, and
 * we know the new position @f$ \vec{x}^{*} @f$ using the same prescription with
 * the new time @f$ t^{*} @f$. Assuming that we know both the time @f$ t @f$
 * used in computing @f$ \vec{x} @f$ and the new time @f$ t^{*} @f$, the shift
 * that we need to add to the position @f$ \vec{x} @f$ of the track to relocate
 * it to @f$ \vec{x}^{*} @f$ is
 * @f$ \Delta \vec{x} = \vec{x}^{*} - \vec{x} = \left[ (t^{*} - t^{*}_{R}) - (t - t_{R}) \right] \vec{v}_{d} @f$
 * where @f$ t^{*}_{R} @f$ is the reference time for @f$ t^{*} @f$, and if we
 * can obtain the new time in the same reference as the old one, we get
 * @f$ \Delta \vec{x} = (t^{*} - t) \vec{v}_{d} @f$.
 * This was a long trip to get to a known place. The most important take-away is
 * still not completely apparent from the latter formula: the drift velocity has
 * its sign opposite to the distance of @f$ \vec{x} @f$ from the anode, meaning
 * that if the new time is later than the old one, @f$ \Delta \vec{x} @f$ is
 * "negative" and the track is moved closer to the anode. The other important
 * take-away is the relevance of the reference for the two times: this
 * prescription spells out the admittedly unsurprising way to deal with them,
 * and more importantly stresses that we do need to deal with them.
 * 
 * 
 * Cathode-crossing tracks
 * ------------------------
 * 
 * The parts at the two sides of the cathode of a trajectory will be pulled to
 * opposite absolute directions by a shift in time (which is the principle on
 * which the timing of cathode-crossing tracks is based on). For instance,
 * postponing the time of the track (@f$ t^{*} > t @f$ in the notation in the
 * previous section) will move each point closer to the anode of its TPC,
 * effectively tearing the track apart.
 * 
 * This algorithm does not try to reason with this feature, it just applies it.
 * It is assumed that the choice of the new time is sensible.
 * Note however that if for example the original time was biassed by a
 * misestimation of the drift velocity, a shift is not enough to fix the
 * reconstruction, which uses the drift velocity as a key ingredient.
 * 
 * 
 * Multithreading notes
 * =====================
 * 
 * This algorithm is not designed for general multithreading, but it can be used
 * in multithreading as long as the state (`setup()`) does not need to be
 * changed.
 * 
 */
class icarus::TrackShifterAlg: private icarus::ns::util::mfLoggingClass {
  
  // --- BEGIN --  Setup  ------------------------------------------------------
  
  geo::GeometryCore const* fGeom = nullptr; ///< Geometry service provider.
  
  /// Detector properties information.
  std::optional<detinfo::DetectorPropertiesData> fDetProp;
  // "optional" works around the absence of default constructor
  
  // --- END ----  Setup  ------------------------------------------------------
  
    public:
  
  using nanoseconds = util::quantities::intervals::nanoseconds; // alias
  
  
  /**
   * @brief Algorithm constructor: configuration.
   * @param logCategory (default: `"TrackShifterAlg"`) name of the category
   *                    used for console messages
   * 
   * The algorithm will need to be set up before being usable.
   */
  TrackShifterAlg(std::string const& logCategory = "TrackShifterAlg");
  
  /**
   * @brief Algorithm constructor: configures and sets context information.
   * @param logCategory (default: `"TrackShifterAlg"`) name of the category
   *                    used for console messages
   * @param setupArg first argument for `setup()` call
   * @param moreSetupArgs the rest of the arguments for `setup()` call
   * @see `TrackShifterAlg(std::string const&)`, `setup()`
   * 
   * This constructor is a one-stop configure + setup call.
   * 
   * The first argument of the constructor is the same as the one of the other,
   * configuration-only constructor (which is called).
   * The strange form of the two last arguments is just a C++ cheat.
   * The actual arguments for `setup()` are documented in that function,
   * and all the arguments specified here will be forwarded to `setup()`.
   */
  template <typename SetupArg, typename... MoreSetupArgs>
  TrackShifterAlg(
    std::string const& logCategory,
    SetupArg&& setupArg, MoreSetupArgs&&... moreSetupArgs
    );
  
  
  // --- BEGIN -- Setup --------------------------------------------------------
  /// @name Algorithm setup
  /// @{
  /**
   * @brief Sets up the algorithm.
   * @param geom geometry service provider
   * @param detProp detector properties information
   * 
   * The `setup()` call supersedes all the previous ones.
   * If the object is not set up, it will refuse to run the algorithm.
   */
  void setup(
    geo::GeometryCore const& geom,
    detinfo::DetectorPropertiesData detProp
    );
  
  /// Undoes the setup. Another `setup()` will be needed to use the algorithm.
  void unsetup();
  
  /// Returns whether the algorithm has been set up (at least once).
  bool isSetup() const;
  
  /// @}
  // --- END ---- Setup --------------------------------------------------------
  
  
  //@{
  /// Returns a new track shifted by the specified amount of time.
  recob::Track shiftTrack(recob::Track const& track, nanoseconds shift) const;
  recob::Track operator()(recob::Track const& track, nanoseconds shift) const
    { return shiftTrack(track, shift); }
  //@}
  
  
    private:
  
  /// Returns the most representative cryostat for the specified `track`.
  geo::CryostatGeo const* trackCryostat
    (recob::TrackTrajectory const& track) const;

  /// Returns a pointer to the cryostat closest to the specified `point`.
  geo::CryostatGeo const* findNearestCryostat(geo::Point_t const& point) const;
  
  /// Returns a pointer to the TPC in `cryo` closest to the specified `point`.
  geo::TPCGeo const* findNearestTPC
    (geo::Point_t const& point, geo::CryostatGeo const& cryo) const;
  
  /// Returns a pointer to the box closest to `point`, `nullptr` if no boxes.
  /// 
  /// The boxes must support a call to `Center()` returning a `geo::Point_t`.
  template <typename Boxes>
  static typename Boxes::value_type const* findNearestBox
    (geo::Point_t const& point, Boxes const& boxes);
  
}; // icarus::TrackShifterAlg



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename SetupArg, typename... MoreSetupArgs>
icarus::TrackShifterAlg::TrackShifterAlg(
  std::string const& logCategory,
  SetupArg&& setupArg, MoreSetupArgs&&... moreSetupArgs
  )
  : TrackShifterAlg{ logCategory }
{
  // this is just not to chase in the constructor all the changes in the setup:
  setup(
    std::forward<SetupArg>(setupArg),
    std::forward<MoreSetupArgs>(moreSetupArgs)...
    );
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_TRACKSHIFTERALG_H
