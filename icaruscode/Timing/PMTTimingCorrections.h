///////////////////////////////////////////////////////////////////////
///
/// \file   icaruscode/Timing/PMTTimingCorrections.h
///
/// \brief  Interface class between the calibration database and the PMT time corrections
///
/// \author A. Scarpelli
///
/// \mailto ascarpell@bnl.gov
///
////////////////////////////////////////////////////////////////////////

#ifndef ICARUSCODE_TIMING_PMTTIMINGCORRECTIONS_H
#define ICARUSCODE_TIMING_PMTTIMINGCORRECTIONS_H


#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"

namespace icarusDB {

  class PMTTimingCorrections: lar::UncopiableClass
  {
    public:
      
      virtual ~PMTTimingCorrections() noexcept = default;
      
      /**
       * @brief Returns the correction on global trigger signal digitization.
       * @returns the correction on global trigger signal digitization [us]
       * 
       * Corrections are saved with the sign corresponding to their time
       * direction: they are expected to be _added_ to a time to correct it.
       */
      virtual double getTriggerCableDelay( unsigned int channelID ) const = 0;
      
      /**
       * @brief Returns the correction on PMT readout board time counter.
       * @return correction on readout board time counter [us]
       * 
       * Corrections are saved with the sign corresponding to their time
       * direction: they are expected to be _added_ to a time to correct it.
       */
      virtual double getResetCableDelay( unsigned int channelID ) const = 0;

      virtual double getLaserCorrections( unsigned int channelID ) const = 0;

      virtual double getCosmicsCorrections( unsigned int channelID ) const = 0;

  }; // end class

}// end of namespace

DECLARE_ART_SERVICE_INTERFACE(icarusDB::PMTTimingCorrections, SHARED)

#endif
