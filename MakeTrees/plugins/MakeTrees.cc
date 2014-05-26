// -*- C++ -*-
//
// Package:    L1Tutorial/MakeTrees
// Class:      MakeTrees
// 
/**\class MakeTrees MakeTrees.cc L1Analyzer/MakeTrees/plugins/MakeTrees.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Pablo Martinez (ETHZ) [pablom]
//         Created:  Tue, 13 May 2014 07:56:56 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#define MAX_OBJECT 10


//
// class declaration
//

//========================Main class===========================//
class MakeTrees : public edm::EDAnalyzer {
   public:
      explicit MakeTrees(const edm::ParameterSet&);
      ~MakeTrees();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      std::vector<const reco::Candidate*> getCollections(const edm::Event&, edm::InputTag);
      int MatchingJets(std::vector<const reco::Candidate*>, std::vector<const reco::Candidate*>, float);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void setBranches();

        
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      TTree *outputTree;
      bool doReco;
};
//====================End of Main class========================//


//=====================Store L1 object=========================//
class L1Object {

  public:

    L1Object(){};
    ~L1Object(){};
    void reset();
    
    int nMuon;
    int nCentralJet;
    int nForwardJet;
    int nTau;
    int nIsoEG;
    int nNonIsoEG;
    int nMET;
    int nMHT;
    int nPFJet;

    int passes_Jet_Met_trigger;
    int match_L1Jet_Reco; 
    float matched_pt; 

    float ptMuon[MAX_OBJECT];
    float etMuon[MAX_OBJECT];
    float phiMuon[MAX_OBJECT];
    float etaMuon[MAX_OBJECT];
  
    float ptCentralJet[MAX_OBJECT];
    float etCentralJet[MAX_OBJECT];
    float phiCentralJet[MAX_OBJECT];
    float etaCentralJet[MAX_OBJECT];
  
    float ptForwardJet[MAX_OBJECT];
    float etForwardJet[MAX_OBJECT];
    float phiForwardJet[MAX_OBJECT];
    float etaForwardJet[MAX_OBJECT];
    
    float ptPFJet[MAX_OBJECT];
    float etPFJet[MAX_OBJECT];
    float phiPFJet[MAX_OBJECT];
    float etaPFJet[MAX_OBJECT];
  
    float ptTau[MAX_OBJECT];
    float etTau[MAX_OBJECT];
    float phiTau[MAX_OBJECT];
    float etaTau[MAX_OBJECT];
  
    float ptIsoEG[MAX_OBJECT];
    float etIsoEG[MAX_OBJECT];
    float phiIsoEG[MAX_OBJECT];
    float etaIsoEG[MAX_OBJECT];

    float ptNonIsoEG[MAX_OBJECT];
    float etNonIsoEG[MAX_OBJECT];
    float phiNonIsoEG[MAX_OBJECT];
    float etaNonIsoEG[MAX_OBJECT];
  
    float ptMET[MAX_OBJECT];
    float etMET[MAX_OBJECT];
    float phiMET[MAX_OBJECT];
    float etaMET[MAX_OBJECT];
  
    float ptMHT[MAX_OBJECT];
    float etMHT[MAX_OBJECT];
    float phiMHT[MAX_OBJECT];
    float etaMHT[MAX_OBJECT];

};
//====================End of L1 object=========================//


//
// constants, enums and typedefs
//

//
// static data member definitions
//
class L1Object L1Event;
edm::Service<TFileService> fs;
//
// constructors and destructor
//
MakeTrees::MakeTrees(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  doReco = iConfig.getParameter<bool>("doReco");

}


MakeTrees::~MakeTrees()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //delete outputTree;
  
}


//
// member functions
//



// Turn a set of InputTags into a colleciton of candidate pointers.
std::vector<const reco::Candidate*> MakeTrees::getCollections(const edm::Event& evt, edm::InputTag collection) {
  
  std::vector<const reco::Candidate*> output;

  edm::Handle<edm::View<reco::Candidate> > handle;
  evt.getByLabel(collection, handle);

  // Loop over objects in current collection
  for (size_t j = 0; j < handle->size(); ++j) {
    const reco::Candidate& object = handle->at(j);
    output.push_back(&object);
  }

  return output;

}





int MakeTrees::MatchingJets(std::vector<const reco::Candidate*> l1, std::vector<const reco::Candidate*> reco, float DRMax) {

  //Assuming jets are ordered by pt
  for(size_t j = 0; j < l1.size(); ++j) {
    for(size_t i = 0; i < reco.size(); ++i) {
      float DeltaPhi = l1[j]->phi()-reco[j]->phi();
      float DeltaEta = l1[j]->eta()-reco[j]->eta();
      if(sqrt(DeltaPhi*DeltaPhi + DeltaEta*DeltaEta) < DRMax) return i;
    }
  }
  return -1;

}


bool ptSort(const reco::Candidate * i, const reco::Candidate * j) { return i->pt() > j->pt(); }


// ------------ method called for each event  ------------
void
MakeTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;


  edm::InputTag tagJetCentral("l1extraParticlesUCT", "Central");
  edm::InputTag tagJetForward("l1extraParticlesUCT", "Forward");
  edm::InputTag tagTau("l1extraParticlesUCT", "Tau");
  edm::InputTag tagMuon("l1extraParticlesUCT", "");
  edm::InputTag tagIsoEG("l1extraParticlesUCT", "Isolated");
  edm::InputTag tagNoIsoEG("l1extraParticlesUCT", "NonIsolated");
  edm::InputTag tagMET("l1extraParticlesUCT", "MET");
  edm::InputTag tagMHT("l1extraParticlesUCT", "MHT");

  std::vector<const reco::Candidate*> JetCentral = getCollections(iEvent, tagJetCentral);
  std::vector<const reco::Candidate*> JetForward = getCollections(iEvent, tagJetForward);
  std::vector<const reco::Candidate*> Tau = getCollections(iEvent, tagTau);
  std::vector<const reco::Candidate*> Muon = getCollections(iEvent, tagMuon);
  std::vector<const reco::Candidate*> IsoEG = getCollections(iEvent, tagIsoEG);
  std::vector<const reco::Candidate*> NonIsoEG = getCollections(iEvent, tagNoIsoEG);
  std::vector<const reco::Candidate*> MET = getCollections(iEvent, tagMET);
  std::vector<const reco::Candidate*> MHT = getCollections(iEvent, tagMHT);

  L1Event.reset();


  //JetCentral loop  
  std::cout << "Jet Central: " << JetCentral.size() << std::endl; 
  L1Event.nCentralJet = JetCentral.size();
  for(size_t j = 0; j < JetCentral.size(); ++j) {
    L1Event.ptCentralJet[j] = JetCentral[j]->pt();
    L1Event.etCentralJet[j] = JetCentral[j]->et();
    L1Event.phiCentralJet[j] = JetCentral[j]->phi();
    L1Event.etaCentralJet[j] = JetCentral[j]->eta();
  }

  //JetForward loop  
  std::cout << "Jet Forward: " << JetForward.size() << std::endl; 
  L1Event.nForwardJet = JetForward.size();
  for(size_t j = 0; j < JetForward.size(); ++j) {
    L1Event.ptForwardJet[j] = JetForward[j]->pt();
    L1Event.etForwardJet[j] = JetForward[j]->et();
    L1Event.phiForwardJet[j] = JetForward[j]->phi();
    L1Event.etaForwardJet[j] = JetForward[j]->eta();
  }

  //Muon loop  
  std::cout << "Muon " << Muon.size() << std::endl; 
  L1Event.nMuon = Muon.size();
  for(size_t j = 0; j < Muon.size(); ++j) {
    L1Event.ptMuon[j] = Muon[j]->pt();
    L1Event.etMuon[j] = Muon[j]->et();
    L1Event.phiMuon[j] = Muon[j]->phi();
    L1Event.etaMuon[j] = Muon[j]->eta();
  }

  //Tau loop  
  std::cout << "Tau " << Tau.size() << std::endl; 
  L1Event.nTau = Tau.size();
  for(size_t j = 0; j < Tau.size(); ++j) {
    L1Event.ptTau[j] = Tau[j]->pt();
    L1Event.etTau[j] = Tau[j]->et();
    L1Event.phiTau[j] = Tau[j]->phi();
    L1Event.etaTau[j] = Tau[j]->eta();
  }

  //IsoEG loop  
  std::cout << "IsoEG " << IsoEG.size() << std::endl; 
  L1Event.nIsoEG = IsoEG.size();
  for(size_t j = 0; j < IsoEG.size(); ++j) {
    L1Event.ptIsoEG[j] = IsoEG[j]->pt();
    L1Event.etIsoEG[j] = IsoEG[j]->et();
    L1Event.phiIsoEG[j] = IsoEG[j]->phi();
    L1Event.etaIsoEG[j] = IsoEG[j]->eta();
  }
 
  //NonIsoEG loop  
  std::cout << "IsoNonEG " << NonIsoEG.size() << std::endl; 
  L1Event.nNonIsoEG = NonIsoEG.size();
  for(size_t j = 0; j < NonIsoEG.size(); ++j) {
    L1Event.ptNonIsoEG[j] = NonIsoEG[j]->pt();
    L1Event.etNonIsoEG[j] = NonIsoEG[j]->et();
    L1Event.phiNonIsoEG[j] = NonIsoEG[j]->phi();
    L1Event.etaNonIsoEG[j] = NonIsoEG[j]->eta();
  }

  //MET loop  
  std::cout << "MET " << MET.size() << std::endl; 
  L1Event.nMET = MET.size();
  for(size_t j = 0; j < MET.size(); ++j) {
    L1Event.ptMET[j] = MET[j]->pt();
    L1Event.etMET[j] = MET[j]->et();
    L1Event.phiMET[j] = MET[j]->phi();
    L1Event.etaMET[j] = MET[j]->eta();
  }

  //MHT loop  
  std::cout << "MHT " << MHT.size() << std::endl; 
  L1Event.nMHT = MHT.size();
  for(size_t j = 0; j < MHT.size(); ++j) {
    L1Event.ptMHT[j] = MHT[j]->pt();
    L1Event.etMHT[j] = MHT[j]->et();
    L1Event.phiMHT[j] = MHT[j]->phi();
    L1Event.etaMHT[j] = MHT[j]->eta();
  }




  //An easy toy-jet-met trigger
  float jet_threshold = 50.0;
  float met_threshold = 100.0;

  bool isValidJet = false;
  bool isValidMET = false;

  for(size_t j = 0; j < JetCentral.size(); ++j) {
    if(JetCentral[j]->pt() > jet_threshold) {
      isValidJet = true;
      break;
    }
  }
  for(size_t j = 0; j < MET.size(); ++j) {
    if(MET[j]->pt() > met_threshold) {
      isValidMET = true;
      break;
    }
  }
  
  if(isValidJet && isValidMET) L1Event.passes_Jet_Met_trigger = 1;



  //Now we get PFRecoJets
  if(doReco) {
   
    edm::InputTag tagPFJets("ak5PFJets", "");
    std::vector<const reco::Candidate*> PFJets = getCollections(iEvent, tagPFJets);
    std::vector<const reco::Candidate*> L1Jets;

    std::sort(PFJets.begin(), PFJets.end(), ptSort);

    //PfJet loop  
    L1Event.nPFJet = PFJets.size();
    for(size_t j = 0; j < PFJets.size(); ++j) {
      L1Event.ptPFJet[j] = PFJets[j]->pt();
      L1Event.etPFJet[j] = PFJets[j]->et();
      L1Event.phiPFJet[j] = PFJets[j]->phi();
      L1Event.etaPFJet[j] = PFJets[j]->eta();
    }

    float threshold = 50.0;
    for(size_t j = 0; j < JetCentral.size(); ++j) {
      if(JetCentral[j]->pt() > threshold) L1Jets.push_back(JetCentral[j]);
    }
    for(size_t j = 0; j < JetForward.size(); ++j) {
      if(JetForward[j]->pt() > threshold) L1Jets.push_back(JetForward[j]);
    }

    int jetIndex = MatchingJets(L1Jets, PFJets, 0.5);
    if(jetIndex == -1) {
      L1Event.match_L1Jet_Reco = 0;
      if(PFJets.size() > 0 ) L1Event.matched_pt = PFJets[0]->pt();
    } else {
      L1Event.match_L1Jet_Reco = 1;
      L1Event.matched_pt = PFJets[jetIndex]->pt();
    }
  
  }

  outputTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTrees::beginJob()
{
  outputTree = fs->make<TTree>( "L1event"  , "L1event");
  setBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeTrees::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
MakeTrees::beginRun(edm::Run const &r, edm::EventSetup const &es)
{
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
MakeTrees::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MakeTrees::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MakeTrees::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeTrees::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void MakeTrees::setBranches() {
 
  outputTree->Branch("nMuon",&L1Event.nMuon,"nMuon/I");
  outputTree->Branch("nCentralJet",&L1Event.nCentralJet,"nCentralJet/I");
  outputTree->Branch("nForwardJet",&L1Event.nForwardJet,"nForwardJet/I");
  outputTree->Branch("nTau",&L1Event.nTau,"nTau/I");
  outputTree->Branch("nIsoEG",&L1Event.nIsoEG,"nIsoEG/I");
  outputTree->Branch("nNonIsoEG",&L1Event.nNonIsoEG,"nNonIsoEG/I");
  outputTree->Branch("nMET",&L1Event.nMET,"nMET/I");
  outputTree->Branch("nMHT",&L1Event.nMHT,"nMHT/I");
  outputTree->Branch("nPFJet",&L1Event.nPFJet,"nPFJet/I");

  outputTree->Branch("passes_Jet_Met_trigger", &L1Event.passes_Jet_Met_trigger, "passes_Jet_Met_trigger/I");
  outputTree->Branch("match_L1Jet_Reco", &L1Event.match_L1Jet_Reco, "match_L1Jet_Reco/I");
  outputTree->Branch("matched_pt", &L1Event.matched_pt, "matched_pt/F");

  outputTree->Branch("ptMuon", L1Event.ptMuon, "ptMuon[nMuon]/F");
  outputTree->Branch("etMuon", L1Event.etMuon, "etMuon[nMuon]/F");
  outputTree->Branch("phiMuon", L1Event.phiMuon, "phiMuon[nMuon]/F");
  outputTree->Branch("etaMuon", L1Event.etaMuon, "etaMuon[nMuon]/F");
  
  outputTree->Branch("ptCentralJet", L1Event.ptCentralJet, "ptCentralJet[nCentralJet]/F");
  outputTree->Branch("etCentralJet", L1Event.etCentralJet, "etCentralJet[nCentralJet]/F");
  outputTree->Branch("phiCentralJet", L1Event.phiCentralJet, "phiCentralJet[nCentralJet]/F");
  outputTree->Branch("etaCentralJet", L1Event.etaCentralJet, "etaCentralJet[nCentralJet]/F");

  outputTree->Branch("ptForwardJet", L1Event.ptForwardJet, "ptForwardJet[nForwardJet]/F");
  outputTree->Branch("etForwardJet", L1Event.etForwardJet, "etForwardJet[nForwardJet]/F");
  outputTree->Branch("phiForwardJet", L1Event.phiForwardJet, "phiForwardJet[nForwardJet]/F");
  outputTree->Branch("etaForwardJet", L1Event.etaForwardJet, "etaForwardJet[nForwardJet]/F");
  
  outputTree->Branch("ptPFJet", L1Event.ptPFJet, "ptPFJet[nPFJet]/F");
  outputTree->Branch("etPFJet", L1Event.etPFJet, "etPFJet[nPFJet]/F");
  outputTree->Branch("phiPFJet", L1Event.phiPFJet, "phiPFJet[nPFJet]/F");
  outputTree->Branch("etaPFJet", L1Event.etaPFJet, "etaPFJet[nPFJet]/F");

  outputTree->Branch("ptTau", L1Event.ptTau, "ptTau[nTau]/F");
  outputTree->Branch("etTau", L1Event.etTau, "etTau[nTau]/F");
  outputTree->Branch("phiTau", L1Event.phiTau, "phiTau[nTau]/F");
  outputTree->Branch("etaTau", L1Event.etaTau, "etaTau[nTau]/F");

  outputTree->Branch("ptIsoEG", L1Event.ptIsoEG, "ptIsoEG[nIsoEG]/F");
  outputTree->Branch("etIsoEG", L1Event.etIsoEG, "etIsoEG[nIsoEG]/F");
  outputTree->Branch("phiIsoEG", L1Event.phiIsoEG, "phiIsoEG[nIsoEG]/F");
  outputTree->Branch("etaIsoEG", L1Event.etaIsoEG, "etaIsoEG[nIsoEG]/F");
  
  outputTree->Branch("ptNonIsoEG", L1Event.ptNonIsoEG, "ptNonIsoEG[nNonIsoEG]/F");
  outputTree->Branch("etNonIsoEG", L1Event.etNonIsoEG, "etNonIsoEG[nNonIsoEG]/F");
  outputTree->Branch("phiNonIsoEG", L1Event.phiNonIsoEG, "phiNonIsoEG[nNonIsoEG]/F");
  outputTree->Branch("etaNonIsoEG", L1Event.etaNonIsoEG, "etaNonIsoEG[nNonIsoEG]/F");
  
  outputTree->Branch("ptMET", L1Event.ptMET, "ptMET[nMET]/F");
  outputTree->Branch("etMET", L1Event.etMET, "etMET[nMET]/F");
  outputTree->Branch("phiMET", L1Event.phiMET, "phiMET[nMET]/F");
  outputTree->Branch("etaMET", L1Event.etaMET, "etaMET[nMET]/F");
  
  outputTree->Branch("ptMHT", L1Event.ptMHT, "ptMHT[nMHT]/F");
  outputTree->Branch("etMHT", L1Event.etMHT, "etMHT[nMHT]/F");
  outputTree->Branch("phiMHT", L1Event.phiMHT, "phiMHT[nMHT]/F");
  outputTree->Branch("etaMHT", L1Event.etaMHT, "etaMHT[nMHT]/F");

}


void L1Object::reset() {

  nMuon = 0;
  nCentralJet = 0;
  nForwardJet = 0;
  nTau = 0;
  nIsoEG = 0;
  nNonIsoEG = 0;
  nMET = 0;
  nMHT = 0;
  passes_Jet_Met_trigger = 0;
  match_L1Jet_Reco = 0;
  matched_pt = 0;
  
  for(size_t j = 0; j < MAX_OBJECT; ++j) {

    ptMuon[j] = 0;
    etMuon[j] = 0;
    phiMuon[j] = 0;
    etaMuon[j] = 0;

    ptCentralJet[j] = 0;
    etCentralJet[j] = 0;
    phiCentralJet[j] = 0;
    etaCentralJet[j] = 0;

    ptForwardJet[j] = 0;
    etForwardJet[j] = 0;
    phiForwardJet[j] = 0;
    etaForwardJet[j] = 0;
    
    ptPFJet[j] = 0;
    etPFJet[j] = 0;
    phiPFJet[j] = 0;
    etaPFJet[j] = 0;

    ptTau[j] = 0;
    etTau[j] = 0;
    phiTau[j] = 0;
    etaTau[j] = 0;

    ptIsoEG[j] = 0;
    etIsoEG[j] = 0;
    phiIsoEG[j] = 0;
    etaIsoEG[j] = 0;

    ptNonIsoEG[j] = 0;
    etNonIsoEG[j] = 0;
    phiNonIsoEG[j] = 0;
    etaNonIsoEG[j] = 0;

    ptMET[j] = 0;
    etMET[j] = 0;
    phiMET[j] = 0;
    etaMET[j] = 0;

    ptMHT[j] = 0;
    etMHT[j] = 0;
    phiMHT[j] = 0;
    etaMHT[j] = 0;

  }

}




//define this as a plug-in
DEFINE_FWK_MODULE(MakeTrees);
