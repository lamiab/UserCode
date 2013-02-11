// -*- C++ -*-
//
// Package:    HiZAnalyzer
// Class:      HiZAnalyzer
// 
/**\class HiZAnalyzer HiZAnalyzer.cc UserCode/tdahms/HiAnalysis/HiOnia/plugins/HiZAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Torsten Dahms,40 4-A32,+41227671635,
//         Created:  Mon Nov 29 03:13:35 CET 2010
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "HiAnalysis/HiOnia/interface/MyCommonHistoManager.h"

// adding Event Plane by dmoon 
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

//
// class declaration
//

class HiZAnalyzer : public edm::EDAnalyzer {
public:
  explicit HiZAnalyzer(const edm::ParameterSet&);
  ~HiZAnalyzer();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  
  void InitEvent();
  void InitTree();

  void makeCuts(int sign) ;
  bool checkCuts(const pat::CompositeCandidate* cand, const pat::Muon* muon1,  const pat::Muon* muon2, bool(HiZAnalyzer::* callFunc1)(const pat::Muon*), bool(HiZAnalyzer::* callFunc2)(const pat::Muon*)); 

  void fillGenInfo();
  int getParentId(const reco::GenParticle* p);

  void fillRecoMuons();
  bool isMuonInAccept(const pat::Muon* aMuon);

  bool selGlobalMuon(const pat::Muon* aMuon);
  bool selTrackerMuon(const pat::Muon* aMuon);
  bool selStaMuon(const pat::Muon* aMuon);

  void fillRecoHistos(int lastSign);

  void fillTreeMuon(const pat::Muon* muon, int iType, int trigBits);
  void fillTreeJpsi(int iSign, int count);

  void checkTriggers(const pat::CompositeCandidate* aJpsiCand);

  TLorentzVector lorentzMomentum(const reco::Candidate::LorentzVector& p);
  // ----------member data ---------------------------
  enum StatBins {
    BIN_nEvents = 0,
    BIN_HLT_PAL1DoubleMuOpen = 1,
    BIN_HLT_PAL1DoubleMu0_HighQ = 2,
    BIN_HLT_PAL2DoubleMu3 = 3,
    BIN_HLT_PAMu3 = 4,
    BIN_HLT_PAMu7 = 5,
    BIN_HLT_PAMu12 = 6
  };

  enum dimuonCategories {
    GlbGlb = 0,
    GlbTrk = 1,
    TrkTrk = 2,
    GlbCal = 3,
    TrkCal = 4,
    CalCal = 5,
    StaSta = 6,
    GlbSta = 7
  };

  std::vector<std::string> theCentralities;
  std::vector<std::string> theTriggerNames;
  std::vector<std::string> theSign;

  float etaMax;

  // TFile
  TFile* fOut;

  // TTree
  TTree* myTree;

  TClonesArray* Reco_mu_4mom;
  TClonesArray* Reco_mu_3vec;
  TClonesArray* Reco_QQ_4mom;
  TClonesArray* Reco_QQ_mupl_4mom;
  TClonesArray* Reco_QQ_mumi_4mom;

  TClonesArray* Gen_mu_4mom;
  TClonesArray* Gen_mu_3vec;
  TClonesArray* Gen_QQ_4mom;
  TClonesArray* Gen_QQ_mupl_4mom;
  TClonesArray* Gen_QQ_mumi_4mom;

  static const int Max_QQ_size = 100;
  static const int Max_mu_size = 100;

  int Gen_QQ_size; // number of generated Onia
  int Gen_QQ_mu_size;
  int Gen_QQ_type[100]; // Onia type: prompt, non-prompt, unmatched
  
  int Gen_mu_size; // number of generated muons
  int Gen_mu_charge[100]; // muon charge
  int Gen_mu_type[100]; // muon type: prompt, non-prompt, unmatched

  int Reco_QQ_size;       // Number of reconstructed Onia 
  int Reco_QQ_type[100];   // Onia category: GG, GT, TT
  int Reco_QQ_sign[100];   /* Mu Mu combinations sign:
			     0 = +/- (signal)
			     1 = +/+
			     2 = -/- 
			  */
  int Reco_QQ_trig[100];      // Vector of trigger bits matched to the Onia
  float Reco_QQ_VtxProb[100]; // chi2 probability of vertex fitting 

  int Reco_QQ_mupl_IsTrackerMu[100]; //Cut variables for mupl: isTrackerMuon()
  int Reco_QQ_mupl_MuValidHits[100]; //Number of MuonValidHits
  int Reco_QQ_mupl_StationsMatched[100]; //Muon segments in n muon stations. This implies that the muon is also an arbitrated tracker muon.
  int Reco_QQ_mupl_TrackerHits[100]; //Number of Inner Tracker hits
  int Reco_QQ_mupl_TrackerLayers[100]; //Number of Tracker layers with measurement
  int Reco_QQ_mupl_PixelHits[100]; //Number of valid pixel hits
  int Reco_QQ_mupl_PixelLayers[100]; //Number of Pixel layers with measurement
  float Reco_QQ_mupl_GlobalChi2[100]; //Global Track fit chi2
  float Reco_QQ_mupl_InnerTrackChi2[100]; //Inner Track fit chi2
  float Reco_QQ_mupl_Dz[100]; //Dz respect to primary vertex
  float Reco_QQ_mupl_Dxy[100]; //Dxy respect to primary vertex
  float Reco_QQ_mupl_PtError[100]; //PtError of muon (not divided by pt)

  int Reco_QQ_mumi_IsTrackerMu[100];
  int Reco_QQ_mumi_MuValidHits[100]; //Cut variables for mumi
  int Reco_QQ_mumi_StationsMatched[100];
  int Reco_QQ_mumi_TrackerHits[100];
  int Reco_QQ_mumi_TrackerLayers[100];
  int Reco_QQ_mumi_PixelHits[100];
  int Reco_QQ_mumi_PixelLayers[100];
  float Reco_QQ_mumi_GlobalChi2[100];
  float Reco_QQ_mumi_InnerTrackChi2[100];
  float Reco_QQ_mumi_Dz[100];
  float Reco_QQ_mumi_Dxy[100];
  float Reco_QQ_mumi_PtError[100];

  int Reco_mu_size;           // Number of reconstructed muons
  int Reco_mu_trig[100];      // Vector of trigger bits matched to the muons
  int Reco_mu_charge[100];  // Vector of charge of muons
  int Reco_mu_type[100];  // Vector of type of muon (global=0, tracker=1, calo=2)
  int Reco_mu_IsTrackerMu[100];
  int Reco_mu_MuValidHits[100]; //Cut variables for all reco muons
  int Reco_mu_StationsMatched[100];
  int Reco_mu_TrackerHits[100];
  int Reco_mu_TrackerLayers[100];
  int Reco_mu_PixelHits[100];
  int Reco_mu_PixelLayers[100];
  float Reco_mu_GlobalChi2[100];
  float Reco_mu_InnerTrackChi2[100];
  float Reco_mu_Dz[100];
  float Reco_mu_Dxy[100];
  float Reco_mu_PtError[100];

  // histos
  TH1F* hGoodMuonsNoTrig;
  TH1F* hGoodMuons;
  TH1F* hL1DoubleMuOpen;
  TH1F* hL2DoubleMu3;
  TH1F* hL2Mu20;

  // event counters
  TH1F* hStats;

  // centrality
  TH1F *hCent;

  // number of primary vertices
  TH1F* hPileUp;

  // z vertex distribution
  TH1F* hZVtx;

  // centrality
  CentralityProvider* centralityProvider_;
  int centBin;
  int theCentralityBin;
  int hiNpix, hiNpixelTracks, hiNtracks;
  float hiHF, hiHFplus, hiHFminus, hiHFplusEta4, hiHFminusEta4, hiEB, hiET, hiEE, hiEEplus, hiEEminus, hiZDC, hiZDCplus, hiZDCminus;

  // handles
  edm::Handle<pat::CompositeCandidateCollection> collJpsi;
  edm::Handle<pat::MuonCollection> collMuon;
  edm::Handle<pat::MuonCollection> collMuonNoTrig;

  edm::Handle<reco::GenParticleCollection> collGenParticles;

  // data members
  edm::InputTag       _patMuon;
  edm::InputTag       _patMuonNoTrig;
  edm::InputTag       _patJpsi;
  edm::InputTag       _genParticle;
  edm::InputTag       _thePVs;
  std::string         _histfilename;
  edm::InputTag       _centralitytag;

  std::vector<double> _centralityranges;
  bool           _applycuts;
  bool           _useBS;
  bool           _removeSignal;
  bool           _removeMuons;
  bool           _storeSs;
  bool           _fillTree;
  bool           _fillSingleMuons;
  bool           _isHI;
  bool           _isMC;
  bool           _isPromptMC;

  int _oniaPDG;

  std::vector<unsigned int>                     _thePassedCats[3];
  std::vector<const pat::CompositeCandidate*>   _thePassedCands[3];

  // number of events
  unsigned int nEvents;

  unsigned int runNb;
  unsigned int eventNb;
  unsigned int lumiSection;

  math::XYZPoint RefVtx;
  float zVtx;
  float nPV;

 // Triger stuff
  // PUT HERE THE *LAST FILTERS* OF THE BITS YOU LIKE
  static const unsigned int sNTRIGGERS = 10;
  unsigned int NTRIGGERS;
  // MC 8E29
  bool isTriggerMatched[sNTRIGGERS];
  std::string HLTLastFilters[sNTRIGGERS];
  bool alreadyFilled[sNTRIGGERS];
  int HLTriggers;

  const edm::ParameterSet _iConfig;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HiZAnalyzer::HiZAnalyzer(const edm::ParameterSet& iConfig):
  _patMuon(iConfig.getParameter<edm::InputTag>("srcMuon")),
  _patMuonNoTrig(iConfig.getParameter<edm::InputTag>("srcMuonNoTrig")),
  _patJpsi(iConfig.getParameter<edm::InputTag>("src")),
  _genParticle(iConfig.getParameter<edm::InputTag>("genParticles")),
  _thePVs(iConfig.getParameter<edm::InputTag>("primaryVertexTag")),
  _histfilename(iConfig.getParameter<std::string>("histFileName")),
  _centralitytag(iConfig.getParameter<edm::InputTag>("srcCentrality")),
  _centralityranges(iConfig.getParameter< std::vector<double> >("centralityRanges")),			
  _applycuts(iConfig.getParameter<bool>("applyCuts")),	
  _useBS(iConfig.getParameter<bool>("useBeamSpot")),
  _removeSignal(iConfig.getUntrackedParameter<bool>("removeSignalEvents",false)),
  _removeMuons(iConfig.getUntrackedParameter<bool>("removeTrueMuons",false)),
  _storeSs(iConfig.getUntrackedParameter<bool>("storeSameSign",false)), 
  _fillTree(iConfig.getParameter<bool>("fillTree")),   
  _fillSingleMuons(iConfig.getParameter<bool>("fillSingleMuons")),
  _isHI(iConfig.getUntrackedParameter<bool>("isHI",true) ),
  _isMC(iConfig.getUntrackedParameter<bool>("isMC",false) ),
  _isPromptMC(iConfig.getUntrackedParameter<bool>("isPromptMC",true) ),
  _oniaPDG(iConfig.getParameter<int>("oniaPDG")),
  NTRIGGERS(iConfig.getParameter<uint32_t>("NumberOfTriggers")),
  _iConfig(iConfig)
{
   //now do what ever initialization is needed
  nEvents = 0;
  if (_isHI) {

   centralityProvider_ = 0;

   std::stringstream centLabel;
   for (unsigned int iCent=0; iCent<_centralityranges.size(); ++iCent) {
    if (iCent==0)
      centLabel << "00" << _centralityranges.at(iCent);
    else
      centLabel << _centralityranges.at(iCent-1) << _centralityranges.at(iCent);

    theCentralities.push_back(centLabel.str());
    centLabel.str("");
   }
   theCentralities.push_back("MinBias");
  }

  theSign.push_back("pm");
  if (_storeSs) {
    theSign.push_back("pp");
    theSign.push_back("mm");
  }

  isTriggerMatched[0]=true; // first entry 'hardcoded' true to accept "all" events
  for (unsigned int iTr = 1; iTr<NTRIGGERS; ++iTr) {
    isTriggerMatched[iTr] = false;
  }

  HLTLastFilters[0] = "";
  HLTLastFilters[1] = "hltL1fL1sPAL1DoubleMuOpenL1Filtered0";		// BIT HLT_PAL1DoubleMuOpen
  HLTLastFilters[2] = "hltL1fL1sPAL1DoubleMu0HighQL1FilteredHighQ";	// BIT HLT_PAL1DoubleMu0_HighQ
  HLTLastFilters[3] = "hltL2fL1sPAL2DoubleMu3L2Filtered3";		// BIT HLT_PAL2DoubleMu3
  HLTLastFilters[4] = "hltL3fL2sMu3L3Filtered3";		// BIT HLT_PAMu3
  HLTLastFilters[5] = "hltL3fL2sMu7L3Filtered7";		// BIT HLT_PAMu7
  HLTLastFilters[6] = "hltL3fL2sMu12L3Filtered12";		// BIT HLT_PAMu12

  theTriggerNames.push_back("NoTrigger");
  theTriggerNames.push_back("HLT_PAL1DoubleMuOpen");
  theTriggerNames.push_back("HLT_PAL1DoubleMu0_HighQ");
  theTriggerNames.push_back("HLT_PAL2DoubleMu3");
  theTriggerNames.push_back("HLT_PAMu3");
  theTriggerNames.push_back("HLT_PAMu7");
  theTriggerNames.push_back("HLT_PAMu12");

  etaMax = 2.4;  
}


HiZAnalyzer::~HiZAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HiZAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   using namespace edm;
  InitEvent();

  nEvents++;
  hStats->Fill(BIN_nEvents);

  //std::cout<< "-------------------- Next event: " << nEvents << " --------------------"<<std::endl;
   
  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lumiSection = iEvent.luminosityBlock();
  
  edm::Handle<reco::VertexCollection> privtxs;
  iEvent.getByLabel(_thePVs, privtxs);
  reco::VertexCollection::const_iterator privtx;

  nPV = privtxs->size();
   
  if ( privtxs->begin() != privtxs->end() ) {
    privtx=privtxs->begin();
    RefVtx = privtx->position();
  } else {
    RefVtx.SetXYZ(0.,0.,0.);
  }

  zVtx = RefVtx.Z();

  hZVtx->Fill(zVtx);
  if (fabs(zVtx) > _iConfig.getParameter< double > ("maxAbsZ")) return;
  hPileUp->Fill(nPV);

  if(_isHI) {
    edm::Handle<reco::Centrality> centrality_;
    iEvent.getByLabel(_centralitytag,centrality_);

    if(!centralityProvider_) centralityProvider_ = new CentralityProvider(iSetup);
    centralityProvider_->newEvent(iEvent,iSetup); // make sure you do this first in every event
    centBin = centralityProvider_->getBin();
    hCent->Fill(centBin);

    hiNpix = centrality_->multiplicityPixel();
    hiNpixelTracks = centrality_->NpixelTracks();
    hiNtracks = centrality_->Ntracks();
    hiHF = centrality_->EtHFtowerSum();
    hiHFplus = centrality_->EtHFtowerSumPlus();
    hiHFminus = centrality_->EtHFtowerSumMinus();
    hiHFplusEta4 = centrality_->EtHFtruncatedPlus();
    hiHFminusEta4 = centrality_->EtHFtruncatedMinus();
    hiZDC = centrality_->zdcSum();
    hiZDCplus = centrality_->zdcSumPlus();
    hiZDCminus = centrality_->zdcSumMinus();
    hiEEplus = centrality_->EtEESumPlus();
    hiEEminus = centrality_->EtEESumMinus();
    hiEE = centrality_->EtEESum();
    hiEB = centrality_->EtEBSum();
    hiET = centrality_->EtMidRapiditySum();

    for (unsigned int iCent=0; iCent<_centralityranges.size(); ++iCent) {
     if (centBin<_centralityranges.at(iCent)) {
       theCentralityBin=iCent;
       break;
     }
    }
  }

  iEvent.getByLabel(_patJpsi,collJpsi); 
  iEvent.getByLabel(_patMuon,collMuon);
  iEvent.getByLabel(_patMuonNoTrig,collMuonNoTrig);

  if (_isMC) {
    iEvent.getByLabel(_genParticle,collGenParticles);
    fillGenInfo();
  }
  
  // APPLY CUTS
  int lastSign = 0;
  makeCuts(0);
  if (_storeSs) {
    makeCuts(1);
    makeCuts(2);
    lastSign = 2;
  }

  if (_fillSingleMuons)
    fillRecoMuons();

  fillRecoHistos(lastSign);

  if (_fillTree)
    myTree->Fill();

  return;
}

void
HiZAnalyzer::fillRecoHistos(int lastSign) {
   
    for (int iSign = 0; iSign <= lastSign; ++iSign) {
      for( unsigned int count = 0; count < _thePassedCands[iSign].size(); count++) { 
	const pat::CompositeCandidate* aJpsiCand = _thePassedCands[iSign].at(count); 

	checkTriggers(aJpsiCand);
	if (_fillTree)
	  fillTreeJpsi(iSign, count);
      }
    }

  return;
}

void
HiZAnalyzer::fillTreeMuon(const pat::Muon* muon, int iType, int trigBits) {
  if (Reco_mu_size >= Max_mu_size) {
    std::cout << "Too many muons: " << Reco_mu_size << std::endl;
    std::cout << "Maximum allowed: " << Max_mu_size << std::endl;
    return;
  }

  Reco_mu_charge[Reco_mu_size] = muon->charge();
  Reco_mu_type[Reco_mu_size] = iType;
  
  TLorentzVector vMuon = lorentzMomentum(muon->p4());
  new((*Reco_mu_4mom)[Reco_mu_size])TLorentzVector(vMuon);

  Reco_mu_trig[Reco_mu_size] = trigBits;

  Reco_mu_IsTrackerMu[Reco_mu_size] = muon->isTrackerMuon();
  Reco_mu_MuValidHits[Reco_mu_size] = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
  Reco_mu_StationsMatched[Reco_mu_size] = muon->numberOfMatchedStations();
  Reco_mu_TrackerHits[Reco_mu_size] = muon->innerTrack()->found();
  Reco_mu_TrackerLayers[Reco_mu_size] = muon->globalTrack()->hitPattern().trackerLayersWithMeasurement();
  Reco_mu_PixelHits[Reco_mu_size] = muon->innerTrack()->hitPattern().numberOfValidPixelHits();
  Reco_mu_PixelLayers[Reco_mu_size] = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
  Reco_mu_GlobalChi2[Reco_mu_size] = muon->globalTrack()->chi2()/muon->globalTrack()->ndof();
  Reco_mu_InnerTrackChi2[Reco_mu_size] = muon->innerTrack()->chi2()/muon->innerTrack()->ndof();
  Reco_mu_Dz[Reco_mu_size] = muon->innerTrack()->dz(RefVtx);
  Reco_mu_Dxy[Reco_mu_size] = muon->innerTrack()->dxy(RefVtx);
  Reco_mu_PtError[Reco_mu_size] = muon->innerTrack()->ptError();

  Reco_mu_size++;
  return;
}

void
HiZAnalyzer::fillTreeJpsi(int iSign, int count) {
  if (Reco_QQ_size >= Max_QQ_size) {
    std::cout << "Too many dimuons: " << Reco_QQ_size << std::endl;
    std::cout << "Maximum allowed: " << Max_QQ_size << std::endl;
    return;
  }

  const pat::CompositeCandidate* aJpsiCand = _thePassedCands[iSign].at(count);
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aJpsiCand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aJpsiCand->daughter("muon2"));


  int trigBits=0;
  for (unsigned int iTr=1; iTr<NTRIGGERS; ++iTr) {
    if (isTriggerMatched[iTr])
      trigBits += pow(2,iTr-1);
  }

  Reco_QQ_sign[Reco_QQ_size] = iSign;
  Reco_QQ_type[Reco_QQ_size] = _thePassedCats[iSign].at(count);

  Reco_QQ_trig[Reco_QQ_size] = trigBits;

  TLorentzVector vMuon1 = lorentzMomentum(muon1->p4());
  TLorentzVector vMuon2 = lorentzMomentum(muon2->p4());

  int istrackermu1 = 0;
  int istrackermu2 = 0;
  int muvalidhits1 = -1;
  int muvalidhits2 = -1;
  int stationsmatched1 = -1;
  int stationsmatched2 = -1;
  int trackerhits1 = -1;
  int trackerhits2 = -1;
  int trackerlayers1 = -1;
  int trackerlayers2 = -1;
  int pixelhits1 = -1;
  int pixelhits2 = -1;
  int pixellayers1 = -1;
  int pixellayers2 = -1;
  float glbchi1 = -1;
  float glbchi2 = -1;
  float innchi1 = -1;
  float innchi2 = -1;
  float dz1 = 99;
  float dz2 = 99;
  float dxy1 = 99;
  float dxy2 = 99;
  float pterr1 = -1;
  float pterr2 = -1;

  istrackermu1 = muon1->isTrackerMuon();
  istrackermu2 = muon2->isTrackerMuon();
  muvalidhits1 = muon1->globalTrack()->hitPattern().numberOfValidMuonHits();
  muvalidhits2 = muon2->globalTrack()->hitPattern().numberOfValidMuonHits();
  stationsmatched1 = muon1->numberOfMatchedStations();
  stationsmatched2 = muon2->numberOfMatchedStations();
  trackerhits1 = muon1->innerTrack()->found();
  trackerhits2 = muon2->innerTrack()->found();
  trackerlayers1 = muon1->globalTrack()->hitPattern().trackerLayersWithMeasurement();
  trackerlayers2 = muon2->globalTrack()->hitPattern().trackerLayersWithMeasurement();
  pixelhits1 = muon1->innerTrack()->hitPattern().numberOfValidPixelHits();
  pixelhits2 = muon2->innerTrack()->hitPattern().numberOfValidPixelHits();
  pixellayers1 = muon1->innerTrack()->hitPattern().pixelLayersWithMeasurement();
  pixellayers2 = muon2->innerTrack()->hitPattern().pixelLayersWithMeasurement();
  glbchi1 = muon1->globalTrack()->chi2()/muon1->globalTrack()->ndof();
  glbchi2 = muon2->globalTrack()->chi2()/muon2->globalTrack()->ndof();
  innchi1 = muon1->innerTrack()->chi2()/muon1->innerTrack()->ndof();
  innchi2 = muon2->innerTrack()->chi2()/muon2->innerTrack()->ndof();
  dz1 = muon1->innerTrack()->dz(RefVtx);
  dz2 = muon2->innerTrack()->dz(RefVtx);
  dxy1 = muon1->innerTrack()->dxy(RefVtx);
  dxy2 = muon2->innerTrack()->dxy(RefVtx);
  pterr1 = muon1->innerTrack()->ptError();
  pterr2 = muon2->innerTrack()->ptError();

  if (muon1->charge() > muon2->charge()) {
    new((*Reco_QQ_mupl_4mom)[Reco_QQ_size])TLorentzVector(vMuon1);
    new((*Reco_QQ_mumi_4mom)[Reco_QQ_size])TLorentzVector(vMuon2);
    Reco_QQ_mupl_MuValidHits[Reco_QQ_size]=muvalidhits1;
    Reco_QQ_mumi_MuValidHits[Reco_QQ_size]=muvalidhits2;
    Reco_QQ_mupl_StationsMatched[Reco_QQ_size]=stationsmatched1;
    Reco_QQ_mumi_StationsMatched[Reco_QQ_size]=stationsmatched2;
    Reco_QQ_mupl_TrackerHits[Reco_QQ_size]=trackerhits1;
    Reco_QQ_mumi_TrackerHits[Reco_QQ_size]=trackerhits2;
    Reco_QQ_mupl_TrackerLayers[Reco_QQ_size]=trackerlayers1;
    Reco_QQ_mumi_TrackerLayers[Reco_QQ_size]=trackerlayers2;
    Reco_QQ_mupl_PixelHits[Reco_QQ_size]=pixelhits1;
    Reco_QQ_mumi_PixelHits[Reco_QQ_size]=pixelhits2;
    Reco_QQ_mupl_PixelLayers[Reco_QQ_size]=pixellayers1;
    Reco_QQ_mumi_PixelLayers[Reco_QQ_size]=pixellayers2;
    Reco_QQ_mupl_GlobalChi2[Reco_QQ_size]=glbchi1;
    Reco_QQ_mumi_GlobalChi2[Reco_QQ_size]=glbchi2;
    Reco_QQ_mupl_InnerTrackChi2[Reco_QQ_size]=innchi1;
    Reco_QQ_mumi_InnerTrackChi2[Reco_QQ_size]=innchi2;
    Reco_QQ_mupl_Dz[Reco_QQ_size]=dz1;
    Reco_QQ_mumi_Dz[Reco_QQ_size]=dz2;
    Reco_QQ_mupl_Dxy[Reco_QQ_size]=dxy1;
    Reco_QQ_mumi_Dxy[Reco_QQ_size]=dxy2;
    Reco_QQ_mupl_PtError[Reco_QQ_size]=pterr1;
    Reco_QQ_mumi_PtError[Reco_QQ_size]=pterr2;
    Reco_QQ_mupl_IsTrackerMu[Reco_QQ_size]=istrackermu1;
    Reco_QQ_mumi_IsTrackerMu[Reco_QQ_size]=istrackermu2;
  }
  else {
    new((*Reco_QQ_mupl_4mom)[Reco_QQ_size])TLorentzVector(vMuon2);
    new((*Reco_QQ_mumi_4mom)[Reco_QQ_size])TLorentzVector(vMuon1);
    Reco_QQ_mupl_MuValidHits[Reco_QQ_size]=muvalidhits2;
    Reco_QQ_mumi_MuValidHits[Reco_QQ_size]=muvalidhits1;
    Reco_QQ_mupl_StationsMatched[Reco_QQ_size]=stationsmatched2;
    Reco_QQ_mumi_StationsMatched[Reco_QQ_size]=stationsmatched1;
    Reco_QQ_mupl_TrackerHits[Reco_QQ_size]=trackerhits2;
    Reco_QQ_mumi_TrackerHits[Reco_QQ_size]=trackerhits1;
    Reco_QQ_mupl_TrackerLayers[Reco_QQ_size]=trackerlayers2;
    Reco_QQ_mumi_TrackerLayers[Reco_QQ_size]=trackerlayers1;
    Reco_QQ_mupl_PixelHits[Reco_QQ_size]=pixelhits2;
    Reco_QQ_mumi_PixelHits[Reco_QQ_size]=pixelhits1;
    Reco_QQ_mupl_PixelLayers[Reco_QQ_size]=pixellayers2;
    Reco_QQ_mumi_PixelLayers[Reco_QQ_size]=pixellayers1;
    Reco_QQ_mupl_GlobalChi2[Reco_QQ_size]=glbchi2;
    Reco_QQ_mumi_GlobalChi2[Reco_QQ_size]=glbchi1;
    Reco_QQ_mupl_InnerTrackChi2[Reco_QQ_size]=innchi2;
    Reco_QQ_mumi_InnerTrackChi2[Reco_QQ_size]=innchi1;
    Reco_QQ_mupl_Dz[Reco_QQ_size]=dz2;
    Reco_QQ_mumi_Dz[Reco_QQ_size]=dz1;
    Reco_QQ_mupl_Dxy[Reco_QQ_size]=dxy2;
    Reco_QQ_mumi_Dxy[Reco_QQ_size]=dxy1;
    Reco_QQ_mupl_PtError[Reco_QQ_size]=pterr2;
    Reco_QQ_mumi_PtError[Reco_QQ_size]=pterr1;
    Reco_QQ_mupl_IsTrackerMu[Reco_QQ_size]=istrackermu2;
    Reco_QQ_mumi_IsTrackerMu[Reco_QQ_size]=istrackermu1;
  }
  
  TLorentzVector vJpsi = lorentzMomentum(aJpsiCand->p4());
  new((*Reco_QQ_4mom)[Reco_QQ_size])TLorentzVector(vJpsi);

  Reco_QQ_VtxProb[Reco_QQ_size] = aJpsiCand->userFloat("vProb");

  Reco_QQ_size++;
  return;
}

void
HiZAnalyzer::checkTriggers(const pat::CompositeCandidate* aJpsiCand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(aJpsiCand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(aJpsiCand->daughter("muon2"));

  // Trigger passed
  for (unsigned int iTr = 1; iTr<NTRIGGERS; ++iTr) {
    const pat::TriggerObjectStandAloneCollection mu1HLTMatchesFilter = muon1->triggerObjectMatchesByFilter( HLTLastFilters[iTr] );
    const pat::TriggerObjectStandAloneCollection mu2HLTMatchesFilter = muon2->triggerObjectMatchesByFilter( HLTLastFilters[iTr] );
    
    const pat::TriggerObjectStandAloneCollection mu1HLTMatchesPath = muon1->triggerObjectMatchesByPath( theTriggerNames.at(iTr) );
    const pat::TriggerObjectStandAloneCollection mu2HLTMatchesPath = muon2->triggerObjectMatchesByPath( theTriggerNames.at(iTr) );
    
    bool pass1 = false;
    bool pass2 = false;
    //    if (iTr<7) { // apparently matching by path gives false positives so we use matching by filter for all triggers for which we know the filter name
    pass1 = mu1HLTMatchesFilter.size() > 0;
    pass2 = mu2HLTMatchesFilter.size() > 0;
    //    }
    //    else {
    //   pass1 = mu1HLTMatchesPath.size() > 0;
    //   pass2 = mu2HLTMatchesPath.size() > 0;
    // }

    if (iTr > 3) {  // single triggers here
      isTriggerMatched[iTr] = pass1 || pass2;
    } else {        // double triggers here
      isTriggerMatched[iTr] = pass1 && pass2;
    }
  }

  for (unsigned int iTr=1;iTr<NTRIGGERS;++iTr) {
    if (isTriggerMatched[iTr]) {
      // fill event counting histogram only once per event, also if several muons fired trigger
      if (alreadyFilled[iTr]) continue;
      hStats->Fill(iTr);
      HLTriggers += pow(2,iTr-1);
      alreadyFilled[iTr]=true;
    }
  }

  return;
}

void
HiZAnalyzer::makeCuts(int sign) {

  if (collJpsi.isValid()) {
    for(std::vector<pat::CompositeCandidate>::const_iterator it=collJpsi->begin();
	it!=collJpsi->end(); ++it) {
      
      const pat::CompositeCandidate* cand = &(*it);	
      //if (fabs(cand->rapidity()) >= etaMax) continue;

      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(cand->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(cand->daughter("muon2"));

      //if (fabs(muon1->rapidity()) >= etaMax ||
	//  fabs(muon2->rapidity()) >= etaMax) continue;
      
      bool thisSign = ( (sign == 0 && muon1->charge() + muon2->charge() == 0) || 
			(sign == 1 && muon1->charge() + muon2->charge() == 2) || 
			(sign == 2 && muon1->charge() + muon2->charge() == -2) );

      if (thisSign) {
	
	// global + global?
	if (checkCuts(cand,muon1,muon2,&HiZAnalyzer::selGlobalMuon,&HiZAnalyzer::selGlobalMuon)){
	  _thePassedCats[sign].push_back(GlbGlb);  _thePassedCands[sign].push_back(cand);
	  continue;
	}

	/*
	// global + tracker? (x2)    
	if (checkCuts(cand,muon1,muon2,&HiZAnalyzer::selGlobalMuon,&HiZAnalyzer::selTrackerMuon)){
	  _thePassedCats[sign].push_back(GlbTrk);  _thePassedCands[sign].push_back(cand);
	  continue;
	}

	if (checkCuts(cand,muon2,muon1,&HiZAnalyzer::selGlobalMuon,&HiZAnalyzer::selTrackerMuon)){
	  _thePassedCats[sign].push_back(GlbTrk);  _thePassedCands[sign].push_back(cand);
	  continue;
	}*/

	// tracker + tracker?  
	/*if (checkCuts(cand,muon1,muon2,&HiZAnalyzer::selTrackerMuon,&HiZAnalyzer::selTrackerMuon)){
	  _thePassedCats[sign].push_back(TrkTrk);  _thePassedCands[sign].push_back(cand);
	  continue;
	}*/

	// sta + sta? meaning both are only sta (not global, not tracker)
	/*if (checkCuts(cand,muon1,muon2,&HiZAnalyzer::selStaMuon,&HiZAnalyzer::selStaMuon)){
	  _thePassedCats[sign].push_back(StaSta);  _thePassedCands[sign].push_back(cand);
	  //continue;
	}

	// global + sta (2x)
	if (checkCuts(cand,muon1,muon2,&HiZAnalyzer::selGlobalMuon,&HiZAnalyzer::selStaMuon)){
	  _thePassedCats[sign].push_back(GlbSta);  _thePassedCands[sign].push_back(cand);
	  //continue;
	}

	if (checkCuts(cand,muon2,muon1,&HiZAnalyzer::selGlobalMuon,&HiZAnalyzer::selStaMuon)){
	  _thePassedCats[sign].push_back(GlbSta);  _thePassedCands[sign].push_back(cand);
	  //continue;
	}*/

      }
    }
  }
  
  return;
}


bool
HiZAnalyzer::checkCuts(const pat::CompositeCandidate* cand, const pat::Muon* muon1,  const pat::Muon* muon2, bool(HiZAnalyzer::* callFunc1)(const pat::Muon*), bool(HiZAnalyzer::* callFunc2)(const pat::Muon*)) {
  //if ( (  (this->*callFunc1)(muon1) &&  (this->*callFunc2)(muon2) ) && (!_applycuts || cand->userFloat("vProb") > 0.01) )
  if ( (this->*callFunc1)(muon1) &&  (this->*callFunc2)(muon2) )
    return true;
  else
    return false;
}

bool
HiZAnalyzer::isMuonInAccept(const pat::Muon* aMuon) {
  return (fabs(aMuon->eta()) < 2.4 &&
	  aMuon->pt() >= 10);
}

bool
HiZAnalyzer::selGlobalMuon(const pat::Muon* aMuon) {
  
  if(!aMuon->isGlobalMuon())
    return false;
  
  if(!_applycuts)
    return true;

  reco::TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();

  reco::TrackRef gTrack = aMuon->globalTrack();
  const reco::HitPattern& q = gTrack->hitPattern();

  // Tight muon cuts January 2013
  return (isMuonInAccept(aMuon) &&
	  gTrack->normalizedChi2() < 10. &&
	  q.numberOfValidMuonHits() > 0 &&
	  aMuon->numberOfMatchedStations() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 0.2 &&
	  fabs(iTrack->dz(RefVtx)) < 0.5 &&
 	  p.numberOfValidPixelHits() > 0 &&
	  q.trackerLayersWithMeasurement() > 5 );
}

bool 
HiZAnalyzer::selTrackerMuon(const pat::Muon* aMuon) { //we don't use this now
  
  if(!aMuon->isTrackerMuon())
    return false;

  if(!_applycuts)
    return true;

  reco::TrackRef iTrack = aMuon->innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();

  return (isMuonInAccept(aMuon) &&
	  iTrack->found() > 10 &&
	  iTrack->chi2()/iTrack->ndof() < 4.0 &&
	  aMuon->muonID("TrackerMuonArbitrated") &&
	  aMuon->muonID("TMLastStationAngTight") &&
	  p.pixelLayersWithMeasurement() > 0 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
	  fabs(iTrack->dz(RefVtx)) < 15.0 );
}

bool 
HiZAnalyzer::selStaMuon(const pat::Muon* aMuon) {
  
  if(aMuon->isTrackerMuon() || aMuon->isGlobalMuon() || !aMuon->isStandAloneMuon())
    return false;

  if(!_applycuts)
    return true;

  return (isMuonInAccept(aMuon) &&
	  aMuon->outerTrack()->numberOfValidHits()>6);

}


void
HiZAnalyzer::InitEvent()
{
  for (unsigned int iTr=1;iTr<NTRIGGERS;++iTr) {
    alreadyFilled[iTr]=false;
  }
  HLTriggers = 0;

  _thePassedCats[0].clear();      _thePassedCands[0].clear();
  _thePassedCats[1].clear();      _thePassedCands[1].clear();
  _thePassedCats[2].clear();      _thePassedCands[2].clear();

  Reco_QQ_size = 0;
  Reco_mu_size = 0;

  Gen_QQ_size = 0;
  Gen_QQ_mu_size = 0;
  Gen_mu_size = 0;

  Reco_QQ_4mom->Clear();
  Reco_QQ_mupl_4mom->Clear();
  Reco_QQ_mumi_4mom->Clear();
  Reco_mu_4mom->Clear();
  Reco_mu_3vec->Clear();

  if (_isMC) {
    Gen_QQ_4mom->Clear();
    Gen_QQ_mupl_4mom->Clear();
    Gen_QQ_mumi_4mom->Clear();
    Gen_mu_4mom->Clear();
    Gen_mu_3vec->Clear();
  }

  return;
}

void
HiZAnalyzer::fillGenInfo()
{
  if (Gen_QQ_size >= Max_QQ_size || Gen_QQ_mu_size >= Max_QQ_size) {
    std::cout << "Too many dimuons: " << Gen_QQ_size << " or " << Gen_QQ_mu_size << std::endl;
    std::cout << "Maximum allowed: " << Max_QQ_size << std::endl;
    return;
  }

  if (Gen_mu_size >= Max_mu_size) {
    std::cout << "Too many muons: " << Gen_mu_size << std::endl;
    std::cout << "Maximum allowed: " << Max_mu_size << std::endl;
    return;
  }
  
  bool filledplus = false;
  bool filledminus = false;

  if (collGenParticles.isValid()) {
    for(std::vector<reco::GenParticle>::const_iterator it=collGenParticles->begin();
	it!=collGenParticles->end();++it) {
      const reco::GenParticle* gen = &(*it);

      if (abs(gen->pdgId()) == 13 && gen->status() == 1) {
	Gen_mu_type[Gen_mu_size] = _isPromptMC ? 0 : 1; // prompt: 0, non-prompt: 1
	Gen_mu_charge[Gen_mu_size] = gen->charge();

	TLorentzVector vMuon = lorentzMomentum(gen->p4());
	new((*Gen_mu_4mom)[Gen_mu_size])TLorentzVector(vMuon);

	Gen_mu_size++;

	std::cout << "Muon from: " << getParentId(gen) << endl;

	if (getParentId(gen) == _oniaPDG || fabs(getParentId(gen)) == 15) { //muons from tau are also accepted
	    TLorentzVector vMuonZ = lorentzMomentum(gen->p4());    
	    if (gen->charge() > 0) {
	      new((*Gen_QQ_mupl_4mom)[Gen_QQ_mu_size])TLorentzVector(vMuonZ);
	      //std::cout<<"Positive muon with pt: " << vMuonZ.Pt() << endl;
	      filledplus = true;
	    }
	    else if (gen->charge() < 0) {
	      new((*Gen_QQ_mumi_4mom)[Gen_QQ_mu_size])TLorentzVector(vMuonZ);
	      //std::cout<<"Negative muon with pt: " << vMuonZ.Pt() << endl;
	      filledminus = true;
	    }
	    if (filledplus && filledminus) Gen_QQ_mu_size++;
	}
      }

      if (abs(gen->pdgId()) == _oniaPDG && gen->status() == 3) {
	Gen_QQ_type[Gen_QQ_size] = _isPromptMC ? 0 : 1; // prompt: 0, non-prompt: 1
	
	TLorentzVector vJpsi = lorentzMomentum(gen->p4());
	new((*Gen_QQ_4mom)[Gen_QQ_size])TLorentzVector(vJpsi);

	Gen_QQ_size++;
      }

    }
  }

  if(!(filledplus && filledminus)) std::cout << "Error: one or both muons from Z not filled" << std::endl;
  if(Gen_QQ_size != Gen_QQ_mu_size) std::cout << "Error: gen sizes don't match" << std::endl;

  return;
}

int HiZAnalyzer::getParentId(const reco::GenParticle* p)
{
 int parentId = -99;

 if ( p->numberOfMothers()>0 )
   {
     parentId = p->mother()->pdgId();
     if (parentId == p->pdgId() && p->mother()->numberOfMothers()> 0 )
       {
         parentId = p->mother()->mother()->pdgId();
         if (parentId == p->pdgId() && p->mother()->mother()->numberOfMothers()> 0 )
           {
             parentId = p->mother()->mother()->mother()->pdgId();
           }

       }
   }
 return parentId;
}

void
HiZAnalyzer::fillRecoMuons()
{
  int nL1DoubleMuOpenMuons=0;
  int nL2DoubleMu3Muons=0;
  int nL2Mu20Muons=0;
  int nGoodMuons=0;
  int nGoodMuonsNoTrig=0;

  if (collMuonNoTrig.isValid()) {
    for(std::vector<pat::Muon>::const_iterator it=collMuonNoTrig->begin();
	it!=collMuonNoTrig->end();++it) {
      const pat::Muon* muon = &(*it);

      if (muon->isGlobalMuon() &&
	  selGlobalMuon(muon))
	nGoodMuonsNoTrig++;
    }
  }

  if (collMuon.isValid()) {
    for(vector<pat::Muon>::const_iterator it=collMuon->begin();
	it!=collMuon->end();++it) {
      const pat::Muon* muon = &(*it);

      if (muon->isGlobalMuon() &&
	  selGlobalMuon(muon)) {
	
	nGoodMuons++;

	int trigBits=0;
	for (unsigned int iTr=1; iTr<NTRIGGERS; ++iTr) {
	  const pat::TriggerObjectStandAloneCollection muHLTMatchesFilter = muon->triggerObjectMatchesByFilter(  HLTLastFilters[iTr] );
	  const pat::TriggerObjectStandAloneCollection muHLTMatchesPath = muon->triggerObjectMatchesByPath( theTriggerNames.at(iTr) );

	  // apparently matching by path gives false positives so we use matching by filter for all triggers for which we know the filter name
	  if ( muHLTMatchesFilter.size() > 0 ) {

	    trigBits += pow(2,iTr-1);

	    if (iTr==1) nL1DoubleMuOpenMuons++;
	    if (iTr==3) nL2DoubleMu3Muons++;
	    if (iTr==6) nL2Mu20Muons++;
	  }
	}
	if (_fillTree)
	  fillTreeMuon(muon, 0, trigBits);
      }
    }
  }
  
  hGoodMuonsNoTrig->Fill(nGoodMuonsNoTrig);
  hGoodMuons->Fill(nGoodMuons);
  hL1DoubleMuOpen->Fill(nL1DoubleMuOpenMuons);
  hL2DoubleMu3->Fill(nL2DoubleMu3Muons);
  hL2Mu20->Fill(nL2Mu20Muons);

  return;
}

void
HiZAnalyzer::InitTree()
{
  Reco_mu_4mom = new TClonesArray("TLorentzVector", 100);
  Reco_mu_3vec = new TClonesArray("TVector3", 100);
  Reco_QQ_4mom = new TClonesArray("TLorentzVector",10);
  Reco_QQ_mupl_4mom = new TClonesArray("TLorentzVector",10);
  Reco_QQ_mumi_4mom = new TClonesArray("TLorentzVector",10);

  if (_isMC) {
    Gen_mu_4mom = new TClonesArray("TLorentzVector", 10);
    Gen_mu_3vec = new TClonesArray("TVector3", 2);
    Gen_QQ_4mom = new TClonesArray("TLorentzVector", 1);
    Gen_QQ_mupl_4mom = new TClonesArray("TLorentzVector", 1);
    Gen_QQ_mumi_4mom = new TClonesArray("TLorentzVector", 1);
  }

  myTree = new TTree("myTree","My TTree of dimuons");
  
  myTree->Branch("eventNb", &eventNb,   "eventNb/i");
  myTree->Branch("runNb",   &runNb,     "runNb/i");
  myTree->Branch("LS",      &lumiSection, "LS/i"); 
  myTree->Branch("zVtx",    &zVtx,        "zVtx/F"); 
  myTree->Branch("HLTriggers", &HLTriggers, "HLTriggers/I");
  myTree->Branch("Centrality", &centBin, "Centrality/I");

  myTree->Branch("hiHF",&hiHF,"hiHF/F");
  myTree->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
  myTree->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
  myTree->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
  myTree->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");
  myTree->Branch("hiZDC",&hiZDC,"hiZDC/F");
  myTree->Branch("hiZDCplus",&hiZDCplus,"hiZDCplus/F");
  myTree->Branch("hiZDCminus",&hiZDCminus,"hiZDCminus/F");
  myTree->Branch("hiET",&hiET,"hiET/F");
  myTree->Branch("hiEE",&hiEE,"hiEE/F");
  myTree->Branch("hiEB",&hiEB,"hiEB/F");
  myTree->Branch("hiEEplus",&hiEEplus,"hiEEplus/F");
  myTree->Branch("hiEEminus",&hiEEminus,"hiEEminus/F");
  myTree->Branch("hiNpix",&hiNpix,"hiNpix/I");
  myTree->Branch("hiNpixelTracks",&hiNpixelTracks,"hiNpixelTracks/I");
  myTree->Branch("hiNtracks",&hiNtracks,"hiNtracks/I");

  myTree->Branch("Reco_QQ_size", &Reco_QQ_size,  "Reco_QQ_size/I");
  myTree->Branch("Reco_QQ_type", Reco_QQ_type,   "Reco_QQ_type[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_sign", Reco_QQ_sign,   "Reco_QQ_sign[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_4mom", "TClonesArray", &Reco_QQ_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_mupl_4mom", "TClonesArray", &Reco_QQ_mupl_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_mumi_4mom", "TClonesArray", &Reco_QQ_mumi_4mom, 32000, 0);
  myTree->Branch("Reco_QQ_trig", Reco_QQ_trig,   "Reco_QQ_trig[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_VtxProb", Reco_QQ_VtxProb,   "Reco_QQ_VtxProb[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mupl_IsTrackerMu", Reco_QQ_mupl_IsTrackerMu, "Reco_QQ_mupl_IsTrackerMu[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_IsTrackerMu", Reco_QQ_mumi_IsTrackerMu, "Reco_QQ_mumi_IsTrackerMu[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_MuValidHits", Reco_QQ_mupl_MuValidHits, "Reco_QQ_mupl_MuValidHits[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_MuValidHits", Reco_QQ_mumi_MuValidHits, "Reco_QQ_mumi_MuValidHits[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_StationsMatched",Reco_QQ_mupl_StationsMatched, "Reco_QQ_mupl_StationsMatched[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_StationsMatched",Reco_QQ_mumi_StationsMatched, "Reco_QQ_mumi_StationsMatched[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_TrackerHits", Reco_QQ_mupl_TrackerHits, "Reco_QQ_mupl_TrackerHits[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_TrackerHits", Reco_QQ_mumi_TrackerHits, "Reco_QQ_mumi_TrackerHits[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_TrackerLayers", Reco_QQ_mupl_TrackerLayers, "Reco_QQ_mupl_TrackerLayers[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_TrackerLayers", Reco_QQ_mumi_TrackerLayers, "Reco_QQ_mumi_TrackerLayers[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_PixelHits", Reco_QQ_mupl_PixelHits, "Reco_QQ_mupl_PixelHits[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_PixelHits", Reco_QQ_mumi_PixelHits, "Reco_QQ_mumi_PixelHits[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_PixelLayers", Reco_QQ_mupl_PixelLayers, "Reco_QQ_mupl_PixelLayers[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mumi_PixelLayers", Reco_QQ_mumi_PixelLayers, "Reco_QQ_mumi_PixelLayers[Reco_QQ_size]/I");
  myTree->Branch("Reco_QQ_mupl_GlobalChi2", Reco_QQ_mupl_GlobalChi2, "Reco_QQ_mupl_GlobalChi2[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mumi_GlobalChi2", Reco_QQ_mumi_GlobalChi2, "Reco_QQ_mumi_GlobalChi2[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mupl_InnerTrackChi2", Reco_QQ_mupl_InnerTrackChi2, "Reco_QQ_mupl_InnerTrackChi2[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mumi_InnerTrackChi2", Reco_QQ_mumi_InnerTrackChi2, "Reco_QQ_mumi_InnerTrackChi2[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mupl_Dz", Reco_QQ_mupl_Dz, "Reco_QQ_mupl_Dz[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mumi_Dz", Reco_QQ_mumi_Dz, "Reco_QQ_mumi_Dz[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mupl_Dxy", Reco_QQ_mupl_Dxy, "Reco_QQ_mupl_Dxy[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mumi_Dxy", Reco_QQ_mumi_Dxy, "Reco_QQ_mumi_Dxy[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mupl_PtError", Reco_QQ_mupl_PtError, "Reco_QQ_mupl_PtError[Reco_QQ_size]/F");
  myTree->Branch("Reco_QQ_mumi_PtError", Reco_QQ_mumi_PtError, "Reco_QQ_mumi_PtError[Reco_QQ_size]/F");

  myTree->Branch("Reco_mu_size", &Reco_mu_size,  "Reco_mu_size/I");
  myTree->Branch("Reco_mu_type", Reco_mu_type,   "Reco_mu_type[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_charge", Reco_mu_charge,   "Reco_mu_charge[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_4mom", "TClonesArray", &Reco_mu_4mom, 32000, 0);
  myTree->Branch("Reco_mu_trig", Reco_mu_trig,   "Reco_mu_trig[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_IsTrackerMu", Reco_mu_IsTrackerMu,   "Reco_mu_IsTrackerMu[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_MuValidHits", Reco_mu_MuValidHits, "Reco_mu_MuValidHits[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_StationsMatched",Reco_mu_StationsMatched, "Reco_mu_StationsMatched[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_TrackerHits", Reco_mu_TrackerHits, "Reco_mu_TrackerHits[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_TrackerLayers", Reco_mu_TrackerLayers, "Reco_mu_TrackerLayers[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_PixelHits", Reco_mu_PixelHits, "Reco_mu_PixelHits[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_PixelLayers", Reco_mu_PixelLayers, "Reco_mu_PixelLayers[Reco_mu_size]/I");
  myTree->Branch("Reco_mu_GlobalChi2", Reco_mu_GlobalChi2, "Reco_mu_GlobalChi2[Reco_mu_size]/F");
  myTree->Branch("Reco_mu_InnerTrackChi2", Reco_mu_InnerTrackChi2, "Reco_mu_InnerTrackChi2[Reco_mu_size]/F");
  myTree->Branch("Reco_mu_Dz", Reco_mu_Dz, "Reco_mu_Dz[Reco_mu_size]/F");
  myTree->Branch("Reco_mu_Dxy", Reco_mu_Dxy, "Reco_mu_Dxy[Reco_mu_size]/F");
  myTree->Branch("Reco_mu_PtError", Reco_mu_PtError, "Reco_mu_PtError[Reco_mu_size]/F");

  if (_isMC) {
    myTree->Branch("Gen_QQ_size",      &Gen_QQ_size,    "Gen_QQ_size/I");
    myTree->Branch("Gen_QQ_mu_size",      &Gen_QQ_mu_size,    "Gen_QQ_mu_size/I");
    myTree->Branch("Gen_QQ_type",      Gen_QQ_type,    "Gen_QQ_type[Gen_QQ_size]/I");
    myTree->Branch("Gen_QQ_4mom",      "TClonesArray", &Gen_QQ_4mom, 32000, 0);
    myTree->Branch("Gen_QQ_mupl_4mom", "TClonesArray", &Gen_QQ_mupl_4mom, 32000, 0);
    myTree->Branch("Gen_QQ_mumi_4mom", "TClonesArray", &Gen_QQ_mumi_4mom, 32000, 0);

    myTree->Branch("Gen_mu_size",   &Gen_mu_size,  "Gen_mu_size/I");
    myTree->Branch("Gen_mu_type",   Gen_mu_type,   "Gen_mu_type[Gen_mu_size]/I");
    myTree->Branch("Gen_mu_charge", Gen_mu_charge, "Gen_mu_charge[Gen_mu_size]/I");
    myTree->Branch("Gen_mu_4mom",   "TClonesArray", &Gen_mu_4mom, 32000, 0);
  }

}

// ------------ method called once each job just before starting event loop  ------------
void 
HiZAnalyzer::beginJob()
{
  fOut = new TFile(_histfilename.c_str(), "RECREATE");
  InitTree();

  // book histos
  hGoodMuonsNoTrig = new TH1F("hGoodMuonsNoTrig","hGoodMuonsNoTrig",10,0,10);
  hGoodMuons = new TH1F("hGoodMuons","hGoodMuons",10,0,10);
  hL1DoubleMuOpen = new TH1F("hL1DoubleMuOpen","hL1DoubleMuOpen",10,0,10);
  hL2DoubleMu3    = new TH1F("hL1DoubleMu3","hL2DoubleMu3",10,0,10);
  hL2Mu20         = new TH1F("hL2Mu20","hL1Mu20",10,0,10);
  
  hGoodMuonsNoTrig->Sumw2();
  hGoodMuons->Sumw2();
  hL1DoubleMuOpen->Sumw2();
  hL2DoubleMu3->Sumw2();
  hL2Mu20->Sumw2();

  hStats = new TH1F("hStats","hStats;;Number of Events",20,0,20);
  hStats->GetXaxis()->SetBinLabel(1,"All");
  for (int i=2; i< (int) theTriggerNames.size()+1; ++i) {
    hStats->GetXaxis()->SetBinLabel(i,theTriggerNames.at(i-1).c_str());
  }
  hStats->Sumw2();

  if (_isHI) {
   hCent = new TH1F("hCent","hCent;centrality bin;Number of Events",100,0,100);
   hCent->Sumw2();
  }

  hPileUp = new TH1F("hPileUp","Number of Primary Vertices;n_{PV};counts", 50, 0, 50);
  hPileUp->Sumw2();

  hZVtx = new TH1F("hZVtx","Primary z-vertex distribution;z_{vtx} [cm];counts", 120, -30, 30);
  hZVtx->Sumw2();

  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiZAnalyzer::endJob() {
  std::cout << "Total number of events = " << nEvents << std::endl;

  fOut->cd();
  hStats->Write();
  if (_isHI) hCent->Write();
  hPileUp->Write();
  hZVtx->Write();

  if (_fillTree)
    myTree->Write();

  hGoodMuonsNoTrig->Write();
  hGoodMuons->Write();
  hL1DoubleMuOpen->Write();
  hL2DoubleMu3->Write();
  hL2Mu20->Write();

  return;
}

TLorentzVector
HiZAnalyzer::lorentzMomentum(const reco::Candidate::LorentzVector& p) {
  TLorentzVector res;
  res.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());

  return res;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiZAnalyzer);
