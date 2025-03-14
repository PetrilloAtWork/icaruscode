////////////////////////////////////////////////////////////////////////
// Class:       SimTestPulse
// Module Type: producer
// File:        SimTestPulse_module.cc
//
// Imported to ICARUS on October 21, 2018
//
// Generated at Sat Feb 11 04:35:44 2017 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SummaryData/RunData.h"
#include <memory>

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "ParamHolder.h"

#include <TTree.h>
#include <TFile.h>

class SimTestPulse;

class SimTestPulse : public art::EDProducer
{
public:
    explicit SimTestPulse(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    SimTestPulse(SimTestPulse const &) = delete;
    SimTestPulse(SimTestPulse &&) = delete;
    
    SimTestPulse & operator = (SimTestPulse const &) = delete;
    SimTestPulse & operator = (SimTestPulse &&) = delete;
    
    // Required functions.
    void produce(art::Event & e) override;
    void beginRun(art::Run& run) override;
    
    void beginJob() override;
    void endJob()   override;

private:

    bool fVerbose; ///< verbosity
    
    //
    // Parameters to be configured
    //
    double fTriggerTime;                   ///< Trigger timing in electronics clock [us]
    std::vector<double> fSimTime_v;        ///< Charge timing in electronics clock [us]
    std::vector<double> fY_v;              ///< Charge Y position in detector coordinate [cm]
    std::vector<double> fZ_v;              ///< Charge Z position in detector coordinate [cm]
    std::vector<double> fNumElectrons_v;   ///< Charge amount in electron count
    
    //
    // Parameters to be calculated
    //
    std::vector<int>    fTick_v;           ///< Corresponding tick
    std::vector<int>    fPlane0Channel_v;  ///< Resulting plane 0 channel number where charge is deposited
    std::vector<int>    fPlane1Channel_v;  ///< Resulting plane 1 channel number where charge is deposited
    std::vector<int>    fPlane2Channel_v;  ///< Resulting plane 2 channel number where charge is deposited
    
    std::vector<int>    fPlane0Wire_v;     ///< Resulting plane 0 wire number where charge is deposited
    std::vector<int>    fPlane1Wire_v;     ///< Resulting plane 1 wire number where charge is deposited
    std::vector<int>    fPlane2Wire_v;     ///< Resulting plane 2 wire number where charge is deposited
    
    int _run, _subrun, _event;
    
    TTree*              fTupleTree;        ///< output analysis tree
};

SimTestPulse::SimTestPulse(fhicl::ParameterSet const & p)
  : EDProducer{p}
// Initialize member data here.
{
    produces< std::vector<sim::SimChannel> >();
    produces< std::vector<sim::SimEnergyDeposit> >();
    produces< std::vector<raw::Trigger> >();
    produces< sumdata::RunData, art::InRun >();
    
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTriggerTime = clockData.TriggerTime();

//    fTriggerTime    = p.get< double>              ("TriggerTime_us");
    fSimTime_v      = p.get< std::vector<double> >("SimTimeArray_us");
    fY_v            = p.get< std::vector<double> >("YArray_cm");
    fZ_v            = p.get< std::vector<double> >("ZArray_cm");
    fNumElectrons_v = p.get< std::vector<double> >("NumElectronsArray");
    fVerbose        = p.get< bool>                ("Verbose",false);
    
    assert( fSimTime_v.size() == fY_v.size() &&
            fSimTime_v.size() == fZ_v.size() &&
            fSimTime_v.size() == fNumElectrons_v.size() );
}

void SimTestPulse::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTupleTree = tfs->make<TTree>("simTestPulse", "Tree by SimTestPulse_module");
    fTupleTree->Branch("run",&_run,"run/I");
    fTupleTree->Branch("subrun",&_subrun,"subrun/I");
    fTupleTree->Branch("event",&_event,"event/I");
    fTupleTree->Branch("trigger_time",&fTriggerTime,"trigger_time/D");
    fTupleTree->Branch("charge_time_v","std::vector<double>",&fSimTime_v);
    fTupleTree->Branch("tick_v","std::vector<int>",&fTick_v);
    fTupleTree->Branch("y_v","std::vector<double>",&fY_v);
    fTupleTree->Branch("z_v","std::vector<double>",&fZ_v);
    fTupleTree->Branch("e_v","std::vector<double>",&fNumElectrons_v);
    
    fTupleTree->Branch("ch_plane0","std::vector<int>",&fPlane0Channel_v);
    fTupleTree->Branch("ch_plane1","std::vector<int>",&fPlane1Channel_v);
    fTupleTree->Branch("ch_plane2","std::vector<int>",&fPlane2Channel_v);
    
    fTupleTree->Branch("wire_plane0","std::vector<int>",&fPlane0Wire_v);
    fTupleTree->Branch("wire_plane1","std::vector<int>",&fPlane1Wire_v);
    fTupleTree->Branch("wire_plane2","std::vector<int>",&fPlane2Wire_v);
}

void SimTestPulse::endJob()
{
    alternative::ParamHolder::destroy();
}

void SimTestPulse::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using                                                                   
    art::ServiceHandle<geo::Geometry> geo;
    
    std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));
    
    run.put(std::move(runData), art::fullRun());
    
    return;
}

void SimTestPulse::produce(art::Event & e)
{
    std::unique_ptr<std::vector<raw::Trigger> > trigger_v(new std::vector<raw::Trigger> );
    trigger_v->push_back(raw::Trigger(0,fTriggerTime,fTriggerTime,1));
    
    std::unique_ptr<std::vector<sim::SimChannel> > simch_v(new std::vector<sim::SimChannel> );

    std::unique_ptr<std::vector<sim::SimEnergyDeposit>> simDep_v(new std::vector<sim::SimEnergyDeposit>);
    
    geo::Point_t chargeDepCoords = {0., 0., 0.};

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    
    fPlane0Channel_v.clear();
    fPlane1Channel_v.clear();
    fPlane2Channel_v.clear();
    fPlane0Wire_v.clear();
    fPlane1Wire_v.clear();
    fPlane2Wire_v.clear();
    fTick_v.clear();
    
    auto& pholder = alternative::ParamHolder::get();
    pholder.Clear();

    geo::WireReadoutGeom const& wireReadoutAlg = art::ServiceHandle<geo::WireReadout const>()->Get();
    const geo::PlaneGeo& planeGeo = wireReadoutAlg.Plane(geo::PlaneID(0,0,0)); // Get the coordinates of the first wire plane

    auto planeCoords = planeGeo.GetCenter();
    auto planeNormal = planeGeo.GetNormalDirection();

    for(size_t index=0; index < fSimTime_v.size(); ++index) {
        
        int tdc = fSimTime_v[index] / clockData.TPCClock().TickPeriod();
        fTick_v.push_back(clockData.TPCTDC2Tick(tdc));
        
        if(fVerbose)
            std::cout << "[BUFFOON!] Charge injection id " << index << " @ TDC=" << tdc
            << " @ Y=" << fY_v[index] << " @ Z=" << fZ_v[index]
            << " with " << fNumElectrons_v[index] << " electrons " << std::endl;
        
        if(tdc<0) {
            std::cerr << "\033[93m[WARNING BUFFOON!]\033[00m ignoring fSimTime " << fSimTime_v[index]
            << " as it results in negative TDC (invalid)" << std::endl;
            continue;
        }

        // Assume C=0, T=1
        chargeDepCoords = geo::Point_t(planeCoords.X() + planeNormal.X(),fY_v[index],fZ_v[index]);
        
        alternative::TruthHit pulse_record;
        pulse_record.tdc = tdc;
        pulse_record.num_electrons = fNumElectrons_v[index];
        pulse_record.tick = clockData.TPCTDC2Tick(tdc);

        double nElecADC = detProp.ElectronsToADC() * fNumElectrons_v[index];
        double eDeposit = 23.6 * fNumElectrons_v[index] / 1000.;    // This is a guesstimate for now but hopefully not far off, note we assume wirecell will handle recombination

        std::cout << "==> x position of plane 0: " << planeCoords.X() << ", normal: " << planeNormal.X() << ", nElec: " << fNumElectrons_v[index] << ", nElecADC: " << nElecADC << ", edep: " << eDeposit << std::endl;
        std::cout << "    x position of charge deposit: " << chargeDepCoords.X() << std::endl;

        simDep_v->emplace_back(0,
                               fNumElectrons_v[index],
                               0.,
                               eDeposit,
                               geo::Point_t(chargeDepCoords.X(),chargeDepCoords.Y()+0.001,chargeDepCoords.Z()+0.001),
                               geo::Point_t(chargeDepCoords.X(),chargeDepCoords.Y()-0.001,chargeDepCoords.Z()-0.001));

        geo::PlaneID collectionPlaneID(0,0,2);

        for(size_t plane=0; plane<3; ++plane) {
            geo::PlaneID planeID(0,0,plane);
            std::cout << "++ Searching for nearest wire id for planeID: " << planeID << ", xyz: "  << chargeDepCoords.X() << "," << chargeDepCoords.Y() << "," << chargeDepCoords.Z() << std::endl;
            geo::WireID nearestID = wireReadoutAlg.NearestWireID(chargeDepCoords,planeID);
            std::cout << "   WireID: " << nearestID << std::endl;
            std::cout << "** Searching for nearest channel with plane: " << plane << ", xyz: " << chargeDepCoords.X() << "," << chargeDepCoords.Y() << "," << chargeDepCoords.Z() << std::endl;
            auto channel = wireReadoutAlg.NearestChannel(chargeDepCoords,planeID);
            std::cout << "   Returned with channel: " << channel << std::endl;
            auto wire = wireReadoutAlg.ChannelToWire(channel).front().Wire;
            std::cout << "   which gives wire: " << wire << std::endl;
            if(fVerbose) std::cout << "[BUFFOON!]    plane " << plane << " channel "
                << channel << " ... wire " << wire << std::endl;
            pulse_record.channel_list[plane] = channel;
            switch(plane) {
                case 0: fPlane0Channel_v.push_back(channel); fPlane0Wire_v.push_back(wire); break;
                case 1: fPlane1Channel_v.push_back(channel); fPlane1Wire_v.push_back(wire); break;
                case 2: fPlane2Channel_v.push_back(channel); fPlane2Wire_v.push_back(wire); break;
                default:
                    std::cerr << "[BUFFOON!] unexpected plane! " << plane << std::endl;
                    throw std::exception();
            }
            sim::SimChannel sch(channel);
            unsigned planeTDC = clockData.TPCTick2TDC(pulse_record.tick + detProp.GetXTicksOffset(planeID) - detProp.GetXTicksOffset(collectionPlaneID));
            double const xyz[3] = {chargeDepCoords.X(), chargeDepCoords.Y(), chargeDepCoords.Z()};
            sch.AddIonizationElectrons(1,   /// track id, keep 0 = invalid
                                       (unsigned int)planeTDC,
                                       fNumElectrons_v[index],
                                       xyz,
                                       100.); /// energy, keep -1 = invalid
            simch_v->emplace_back(std::move(sch));
        }
        pholder.Register(std::move(pulse_record));
    }

    _run = e.id().run();
    _subrun = e.id().subRun();
    _event = e.id().event();
    
    fTupleTree->Fill();
    e.put(std::move(simch_v));
    e.put(std::move(simDep_v));
    e.put(std::move(trigger_v));
    
    return;
}

DEFINE_ART_MODULE(SimTestPulse)
