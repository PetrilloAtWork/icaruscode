/**
 * @file   icaruscode/Timing/IPMTTimingCorrectionService.h
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Andrea Scarpelli (ascarpell@bnl.gov), Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_TIMING_IPMTTIMINGCORRECTIONSERVICE_H
#define ICARUSCODE_TIMING_IPMTTIMINGCORRECTIONSERVICE_H

// ICARUS libraries
#include "icaruscode/Timing/PMTTimingCorrections.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceProviderWrappers.h"


// -----------------------------------------------------------------------------
namespace icarusDB {
  /**
   * @brief Provides access to PMT timing corrections.
   * 
   * The only thing this service does is to return its service provider of type
   * `icarusDB::PMTTimingCorrections`, which can then be used to query for
   * PMT timing corrections.
   * 
   * Use it to access the service provider like in:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * icarusDB::PMTTimingCorrections const* PMTtimingCorr
   *   = lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  using IPMTTimingCorrectionService
    = lar::ServiceProviderInterfaceWrapper<PMTTimingCorrections>;
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE(icarusDB::IPMTTimingCorrectionService, SHARED)


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_TIMING_IPMTTIMINGCORRECTIONSERVICE_H
