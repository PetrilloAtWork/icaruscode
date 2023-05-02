/**
 * @file   icaruscode/Timing/PMTTimingCorrectionsProvider.h
 * @brief  Service for the PMT timing corrections.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 */

#ifndef ICARUSCODE_TIMING_PMTIMINGCORRECTIONSPROVIDER_H
#define ICARUSCODE_TIMING_PMTIMINGCORRECTIONSPROVIDER_H

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "icaruscode/Timing/PMTTimingCorrections.h"

// Database interface helpers
#include "wda.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

namespace icarusDB::details {
    
  struct PMTTimeCorrectionsDB {

    double triggerCableDelay=0; ///< Expect microseconds.

    double resetCableDelay=0; ///< Expect microseconds.

    double laserCableDelay=0; ///< Expect microseconds.

    double cosmicsCorrections=0; ///< Expect microseconds.
    
  };
  
} // icarusDB::details

namespace icarusDB{ class PMTTimingCorrectionsProvider; }
/**
 * @brief 
 * 
 * This service provider reads the PMT timing correction databases on a
 * run by run basis, and it returns corrections by the PMT channel.
 * 
 * This is an implementation of the `icarusDB::PMTTimingCorrections` interface.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `DatabaseURL` (text; default:
 *   `"https://dbdata0vm.fnal.gov:9443/icarus_con_prod/app/data?"`)
 *   URL of the database with trailing `?` for query.
 * * `Timeout` (seconds; default: `1000`): timeout on the response to database
 *   queries.
 * * `Verbose` (flag; default: `false`): if set, prints all the data loaded
 *   from the database.
 * * `LogCategory` (text, defauly: `"PMTTimingCorrection"`): name of the message
 *   facility stream to use for messages to console.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * This service provider is compatible with _art_ multithreading which requires
 * thread-safety only within the run boundaries: information is cached run by
 * run, and the information from one run replaces the one from the previous one.
 * 
 * The provider guarantees that calling exclusively the `const` interface is
 * thread-safe, while no guarantee is provided while the non-`const` interface
 * is called.
 */
class icarusDB::PMTTimingCorrectionsProvider : public PMTTimingCorrections {

    public: 

        PMTTimingCorrectionsProvider(const fhicl::ParameterSet& pset);

        void readTimeCorrectionDatabase(const art::Run& run);

        /**
         * @brief Returns the correction on global trigger signal digitization.
         * @return correction for global trigger signal digitization delay [us]
         * 
         * This correction removes the delay from global trigger signal
         * generation to its digitization.
         * 
         * This is the delay due to the cables conveying the global trigger
         * signal from the FPGA of the trigger crate that generates it
         * ("global", "third" FPGA) to the dedicated channel of the first
         * digitizer in each VME crate.
         * 
         * It includes an additional correction, the "phase correction", a fudge
         * factor holding local variation due to constant clock offsets
         * (it can be up to 8 ns).
         * It can be absorbed within other corrections if necessary.
         * 
         * Correction are saved with the sign corresponding to their time
         * direction: they are expected to be _added_ to a time to correct it.
         */
        double getTriggerCableDelay( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).triggerCableDelay;
        };

        /**
         * @brief Returns the correction on readout board time counter.
         * @return correction on readout board time counter [us]
         * 
         * This is the delay along the distribution line of the signal to reset
         * the TTT counters of the V1730 readout boards.
         * 
         * It includes an additional correction, the "phase correction", a fudge
         * factor holding local variation due to constant clock offsets
         * (it can be up to 8 ns).
         * It can be absorbed within other corrections if necessary.
         * 
         * Correction are saved with the sign corresponding to their time
         * direction: they are expected to be _added_ to a time to correct it.
         */
        
        double getResetCableDelay( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).resetCableDelay;
        };

        double getLaserCorrections( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).laserCableDelay;
        };

        double getCosmicsCorrections( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).cosmicsCorrections;
        };

    private:
        
        using PMTTimeCorrectionsDB = details::PMTTimeCorrectionsDB;

        std::string fUrl;

        unsigned int fTimeout;

        bool fCorrectCablesDelay;

        bool fVerbose = false; ///< Whether to print the configuration we read.
  
        std::string fLogCategory; ///< Category tag for messages.

        /// Interface to LArSoft configuration for detector timing.
        detinfo::DetectorClocksData const fClocksData;

        static constexpr PMTTimeCorrectionsDB CorrectionDefaults {};
        

        std::map<unsigned int, PMTTimeCorrectionsDB> fDatabaseTimingCorrections;
        
        
        /// Internal access to the channel correction record; returns defaults if not present.
        PMTTimeCorrectionsDB const& getChannelCorrOrDefault
            (unsigned int channelID) const
            {
                auto const it = fDatabaseTimingCorrections.find(channelID);
                return (it == fDatabaseTimingCorrections.end())? CorrectionDefaults: it->second;
            }

        int ConnectToDataset(const std::string& name, 
            uint32_t run, Dataset& dataset ) const;

        void ReadPMTCablesCorrections(uint32_t run);

        void ReadLaserCorrections(uint32_t run);

        void ReadCosmicsCorrections(uint32_t run);

}; // services class

#endif 
