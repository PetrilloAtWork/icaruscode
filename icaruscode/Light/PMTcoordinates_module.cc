////////////////////////////////////////////////////////////////////////
// Class:       PMTcoordinates
// Plugin Type: analyzer (art v2_07_03)
// File:        PMTcoordinates_module.cc
//
// Generated at Mon Sep 25 15:10:31 2017 by Andrea Falcone using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larana/OpticalDetector/SimPhotonCounter.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "lardataobj/Simulation/SimChannel.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int nPMTs = 360;
const int PMTs_per_TPC = 90;
const int MaxPhotons = 10000;

namespace icarus {
  class PMTcoordinates;
}


class icarus::PMTcoordinates : public art::EDAnalyzer 
{
public:
explicit PMTcoordinates(fhicl::ParameterSet const & p);

// Plugins should not be copied or assigned.
PMTcoordinates(PMTcoordinates const &) = delete;
PMTcoordinates(PMTcoordinates &&) = delete;
PMTcoordinates & operator = (PMTcoordinates const &) = delete;
PMTcoordinates & operator = (PMTcoordinates &&) = delete;

// Required functions.
void analyze(art::Event const & e) override;

// Selected optional functions.
void beginJob();
//void reconfigure(fhicl::ParameterSet const & p);

private:

TTree* fTree;

int event;

int noPMT[nPMTs];

int turned_PMT;

double PMTx[nPMTs];
double PMTy[nPMTs];
double PMTz[nPMTs];

int Cryostat[nPMTs];	
int TPC[nPMTs];

float photons_collected[nPMTs];
float photon_time[nPMTs][MaxPhotons];
float photon_velocity[nPMTs][MaxPhotons];

float ics[nPMTs][MaxPhotons];
float ipsilon[nPMTs][MaxPhotons];
float zeta[nPMTs][MaxPhotons];

float dis[nPMTs][MaxPhotons];

float firstphoton_time[nPMTs];
float firstphoton_velocity[nPMTs];

float true_barycentre_y;
float true_barycentre_z;
float reco_barycentre_y;
float reco_barycentre_z;

float total_quenched_energy;
  
art::InputTag photonLabel;
art::InputTag chargeLabel;
};


icarus::PMTcoordinates::PMTcoordinates(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  photonLabel(p.get<art::InputTag>("fottoni", "largeant")),
  chargeLabel(p.get<art::InputTag>("carconi", "largeant"))
 // More initializers here.
{

}

void icarus::PMTcoordinates::analyze(art::Event const & evt)
{

art::ServiceHandle<geo::Geometry> geom;
 
std::vector<sim::SimPhotons> const& optical = *(evt.getValidHandle<std::vector<sim::SimPhotons>>(photonLabel));
std::vector<sim::SimChannel> const& charge  = *(evt.getValidHandle<std::vector<sim::SimChannel>>(chargeLabel));

event = evt.id().event();

for (int g=0; g<10000; g++)
{
	for (int u=0; u<360; u++)
	{
		photon_time[u][g]=0;
		photon_velocity[u][g]=0;

		firstphoton_time[u]=100000000;
		firstphoton_velocity[u]=0;
	}
} 

true_barycentre_y =0;
true_barycentre_z =0;

total_quenched_energy =0;

turned_PMT=0;

reco_barycentre_y=0;
reco_barycentre_z=0;

for (std::size_t chargechannel = 0;  chargechannel<charge.size(); ++chargechannel) //loop sui canali del SimChannel
{ 	
	auto const& channeltdcide = charge.at(chargechannel).TDCIDEMap();
	
	for (std::size_t TDCnu = 0;  TDCnu<channeltdcide.size(); ++TDCnu) 	//loop sui TDC dei canali
	{

		sim::TDCIDE const& tdcide = channeltdcide.at(TDCnu);

		for (std::size_t IDEnu = 0;  IDEnu<tdcide.second.size(); ++IDEnu) 	//loop sulle IDE dei TDC
		{
			sim::IDE const& ida = tdcide.second.at(IDEnu);

			std::cout << "IDA     " << ida.x << '\t' << ida.y << '\t' << ida.z << std::endl;

			true_barycentre_y = true_barycentre_y + ida.y*ida.energy;
			true_barycentre_z = true_barycentre_z + ida.z*ida.energy;
			total_quenched_energy      = total_quenched_energy + ida.energy;

		}
		
	}

}


true_barycentre_y = true_barycentre_y/total_quenched_energy;
true_barycentre_z = true_barycentre_z/total_quenched_energy;

total_quenched_energy = total_quenched_energy/3; 



for (std::size_t channel = 0; channel < optical.size(); ++channel) {

	sim::SimPhotons const& photon_vec = optical[channel];

	noPMT[channel] = channel;	

	photons_collected[channel]=photon_vec.size();

	if (photons_collected[channel]>0){turned_PMT++;}

	double xyz[3];

	geom->OpDetGeoFromOpChannel(channel).GetCenter(xyz);

	PMTx[channel] = xyz[0];
	PMTy[channel] = xyz[1];
	PMTz[channel] = xyz[2];

//	mf::LogVerbatim("PMTcoordinates") << "Channel #" << channel << ": " << photon_vec.size() << " photons";

	firstphoton_time[channel] = 100000000;
	firstphoton_velocity[channel] = 0;

	if (photons_collected[channel]>0)
	{
		for (size_t i = 0; i<photon_vec.size() && int(i)< MaxPhotons; ++i)
		{
			photon_time[channel][i]= photon_vec.at(i).Time;

			ics[channel][i] = photon_vec.at(i).InitialPosition.X()/10;
			ipsilon[channel][i] = photon_vec.at(i).InitialPosition.Y()/10;
			zeta[channel][i] = photon_vec.at(i).InitialPosition.Z()/10;

			dis[channel][i] = sqrt(((ics[channel][i]-PMTx[channel])*(ics[channel][i]-PMTx[channel]))+((ipsilon[channel][i]-PMTy[channel])*(ipsilon[channel][i]-PMTy[channel]))+((zeta[channel][i]-PMTz[channel])*(zeta[channel][i]-PMTz[channel])));

			photon_velocity[channel][i]= dis[channel][i]/photon_time[channel][i];

			if (photon_time[channel][i]<firstphoton_time[channel])
			{				
				firstphoton_time[channel]=photon_time[channel][i];
				firstphoton_velocity[channel]=dis[channel][i]/firstphoton_time[channel];
			}

		}

	}


//	std::cout << PMTx[channel] << '\t' << PMTy[channel] << '\t' << PMTz[channel] << std::endl;

	if (PMTx[channel]<0){Cryostat[channel]=0;}
	if (PMTx[channel]>0){Cryostat[channel]=1;}

	if (PMTx[channel]<-200){TPC[channel]=0;}
	if (PMTx[channel]>-200 && PMTx[channel]<0){TPC[channel]=1;}
	if (PMTx[channel]<200 && PMTx[channel]>0){TPC[channel]=2;}
	if (PMTx[channel]>200){TPC[channel]=3;}
    	
}
//	std::map <int,int> content = optical->DetectedPhotons;
	
//	std::vector<int> channels;
//	std::vector<int> collected_photons;

//	for(std::map<int,int>::iterator it = content.begin(); it != content.end(); ++it) {
 // 	channels.push_back(it->first);
 //	collected_photons.push_back(it->second);
//  	std::cout << it->first << '\t' << it->second << std::endl;}

fTree->Fill();
}

void icarus::PMTcoordinates::beginJob()
{

art::ServiceHandle<art::TFileService> tfs;
fTree = tfs->make<TTree>("lighttree","tree for the light response");

fTree->Branch("event",&event,"event/I");
fTree->Branch("turned_PMT",&turned_PMT,"turned_PMT/I");
fTree->Branch("noPMT",&noPMT,("noPMT[" + std::to_string(nPMTs) + "]/I").c_str());
fTree->Branch("Cryostat",&Cryostat,("Cryostat[" + std::to_string(nPMTs) + "]/I").c_str());
fTree->Branch("TPC",&TPC,("TPC[" + std::to_string(nPMTs) + "]/I").c_str());
fTree->Branch("PMTx",&PMTx,("PMTx[" + std::to_string(nPMTs) + "]/D").c_str());
fTree->Branch("PMTy",&PMTy,("PMTy[" + std::to_string(nPMTs) + "]/D").c_str());
fTree->Branch("PMTz",&PMTz,("PMTz[" + std::to_string(nPMTs) + "]/D").c_str());
fTree->Branch("photons_colleted",&photons_collected,("photons_collected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_time",&firstphoton_time,("firstphoton_time[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("photon_time",&photon_time,"photon_time[360][10000]/F");
fTree->Branch("dis",&dis,"dis[360][10000]/F");
fTree->Branch("firstphoton_velocity",&firstphoton_velocity,("firstphoton_velocity[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("photon_velocity",&photon_velocity,"photon_velocity[360][10000]/F");
fTree->Branch("total_quenced_energy",&total_quenced_energy,"total_quenced_energy");
fTree->Branch("true_barycentre_y",&true_barycentre_y,"true_barycentre_y/F");
fTree->Branch("true_barycentre_z",&true_barycentre_z,"true_barycentre_z/F");
//fTree->Branch("reco_barycentre_y",reco_barycentre_y,"reco_barycentre_y/F");
//fTree->Branch("reco_barycentre_z",reco_barycentre_z,"reco_barycentre_z/F");

}

DEFINE_ART_MODULE(icarus::PMTcoordinates)
