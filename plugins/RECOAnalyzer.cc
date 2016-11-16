// -*- C++ -*-
//
// Package:    TriggerEff/RECOAnalyzer
// Class:      RECOAnalyzer
// 
/**\class RECOAnalyzer RECOAnalyzer.cc TriggerEff/RECOAnalyzer/plugins/RECOAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dylan Sheldon Rankin
//         Created:  Tue, 12 Jul 2016 15:16:15 GMT
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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "DataFormats/METReco/interface/PFMET.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

typedef math::XYZPoint Point;

//
// class declaration
//

class RECOAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RECOAnalyzer(const edm::ParameterSet&);
      ~RECOAnalyzer();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;
      edm::EDGetTokenT<std::vector<reco::PFJet>> jetToken_;
      edm::EDGetTokenT<std::vector<reco::PFMET>> mToken_;
      edm::EDGetTokenT<std::vector<reco::GsfElectron>> elToken_;
      edm::EDGetTokenT<reco::BeamSpot> bsToken_;
      edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
/*
 * These are the collections that we plan to import and use in the rest of the analyzer.
 */

      TH1F *nvtx;
      TH1F *elnum, *elden;
      TH1F *elnum3, *elden3;

      TH2F *elcutflow;
      TH1F *eletanum, *eletaden;
      TH1F *eletanum3, *eletaden3;
      TH2F *eletacutflow;
      TH2F *eldist;
      TH2F *jtdist, *mtdist;
      TH1F *jtpt, *jteta, *mtpt;

/*
 * These are the histograms that we will fill in the analyzer.
 */

      bool debug;

  struct JetRefCompare :
    public std::binary_function<edm::RefToBase<reco::Jet>, edm::RefToBase<reco::Jet>, bool> {
    inline bool operator () (const edm::RefToBase<reco::Jet> &j1,
			     const edm::RefToBase<reco::Jet> &j2) const
    { return j1.id() < j2.id() || (j1.id() == j2.id() && j1.key() < j2.key()); }
  };

  typedef std::map<edm::RefToBase<reco::Jet>, unsigned int, JetRefCompare> FlavourMap;

};

RECOAnalyzer::RECOAnalyzer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    genToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genparticles"))),
    jetToken_(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    mToken_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("met"))),
    elToken_(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
    bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
    convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
    debug(iConfig.getParameter<bool>("debug"))
/*
 * This takes the information given in the cmsRun config file (runEff.py) and stores it in the tokens we defined earlier.
 */

{
    edm::Service<TFileService> fs;
/*
 * This defines the fileservice, which allows us to store the histograms we are going to define in a file.
 */

    nvtx      = fs->make<TH1F>( "nvtx", "nvtx" , 250, 0, 250);
/*
 * This is an example of how to define the histograms we are interested in and specify the name/title/binning
 * The syntax should be self evident (the" fs->make<TH1F>" is required in order for the fileservice file to recognize which histogram you want to save.
 * A similar syntax can be used to save TTrees and other TObjects
 */
    elnum     = fs->make<TH1F>( "elnum", "elnum" , 100, 0, 1000);
    elden     = fs->make<TH1F>( "elden", "elden" , 100, 0, 1000);
    elnum3     = fs->make<TH1F>( "elnum3", "elnum3" , 100, 0, 1000);
    elden3     = fs->make<TH1F>( "elden3", "elden3" , 100, 0, 1000);
/*
 * The "3" suffix here is meant to indicate the genParticle status code used to define the efficiencies
 */
    eletanum     = fs->make<TH1F>( "eletanum", "eletanum" , 60, -3, 3);
    eletaden     = fs->make<TH1F>( "eletaden", "eletaden" , 60, -3, 3);
    eletanum3     = fs->make<TH1F>( "eletanum3", "eletanum3" , 60, -3, 3);
    eletaden3     = fs->make<TH1F>( "eletaden3", "eletaden3" , 60, -3, 3);
    elcutflow = fs->make<TH2F>( "elcutflow", "elcutflow", 50, 0., 500., 10, 0.5, 10.5);
    eletacutflow = fs->make<TH2F>( "eletacutflow", "eletacutflow", 60, -3., 3., 10, 0.5, 10.5);
    eldist = fs->make<TH2F>( "eldist", "eldist", 50, 0., 500., 60, -3., 3.);
    jtdist = fs->make<TH2F>( "jtdist", "jtdist", 50, 0., 2000., 60, -3., 3.);
    jtpt = fs->make<TH1F>(  "jtpt", "jtpt", 50, 0., 2000.);
    jteta = fs->make<TH1F>(  "jteta", "jteta", 60, 5., -5.);
    mtdist = fs->make<TH2F>( "mtdist", "mtdist", 50, 0., 1000., 60, -3., 3.);
    mtpt = fs->make<TH1F>(  "mtpt", "mtpt", 50, 0., 1000.);
}

RECOAnalyzer::~RECOAnalyzer()
{
}

void
RECOAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    
    std::vector<TLorentzVector> rmv;
    std::vector<TLorentzVector> rmtv;
    std::vector<TLorentzVector> rev;
    std::vector<TLorentzVector> retv;
    TLorentzVector tmpvec;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
/*
 * This is the standard way to get the collection you are interested in from the input file.
 */

    if (vertices->empty()) return; // skip the event if no PV found
    nvtx->Fill(vertices->size());
    Point PVtx = vertices->at(0).position();

    edm::Handle<std::vector<reco::PFJet>> jets;
    iEvent.getByToken(jetToken_, jets);
    for (const reco::PFJet &j : *jets) {
        if (j.pt() < 10) continue;
        jtpt->Fill(j.pt());
        jteta->Fill(j.eta());
        jtdist->Fill(j.pt(),j.eta());
    }

    edm::Handle<std::vector<reco::PFMET>> hMet;
    iEvent.getByToken(mToken_, hMet);
    mtpt->Fill(hMet->at(0).pt());
    mtdist->Fill(hMet->at(0).pt(),hMet->at(0).phi());

    //std::cout<<"MET = "<<hMet->at(0).pt()<<std::endl;

    edm::Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(convToken_, conversions);
    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByToken(bsToken_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();

    edm::Handle<std::vector<reco::GsfElectron>> electrons;
    iEvent.getByToken(elToken_, electrons);
    //int ec = 0;
    for (const reco::GsfElectron &e : *electrons) {
        //ec++;
        //std::cout<<"\nEl"<<ec<<" ";
        if (e.pt() < 10) continue;
        //std::cout<<"pt ";
        //std::cout<<"conv ";
        bool fail = false;
        elcutflow->Fill(e.pt(),1);
        eletacutflow->Fill(e.eta(),1);
        double Ooemoop = 999.;
        if (e.ecalEnergy()==0) Ooemoop = 999.;
        else if (!std::isfinite(e.ecalEnergy())) Ooemoop = 998.;
        else Ooemoop = (1.0/e.ecalEnergy() - e.eSuperClusterOverP()/e.ecalEnergy());
        if (fabs(e.superCluster()->eta())<1.4442) {
            if (e.full5x5_sigmaIetaIeta() >= 0.0103) {
                fail = true;
                elcutflow->Fill(e.pt(),2);
                eletacutflow->Fill(e.eta(),2);
            }
            //std::cout<<"sihih ";
            if (fabs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.0105) {
                fail = true;
                elcutflow->Fill(e.pt(),3);
                eletacutflow->Fill(e.eta(),3);
            }
            //std::cout<<"deta ";
            if (fabs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.115) {
                fail = true;
                elcutflow->Fill(e.pt(),4);
                eletacutflow->Fill(e.eta(),4);
            }
            //std::cout<<"dphi ";
            if (e.hcalOverEcal() >= 0.104) {
                fail = true;
                elcutflow->Fill(e.pt(),5);
                eletacutflow->Fill(e.eta(),5);
            }
            //std::cout<<"hoe ";
            if (Ooemoop >= 0.102) {
                fail = true;
                elcutflow->Fill(e.pt(),6);
                eletacutflow->Fill(e.eta(),6);
            }
            //std::cout<<"eop ";
            if (fabs(e.gsfTrack()->dxy(PVtx)) >= 0.0261) {
                fail = true;
                elcutflow->Fill(e.pt(),7);
                eletacutflow->Fill(e.eta(),7);
            }
            //std::cout<<"dxy ";
            if (fabs(e.gsfTrack()->dz(PVtx)) >= 0.41) {
                fail = true;
                elcutflow->Fill(e.pt(),8);
                eletacutflow->Fill(e.eta(),8);
            }
            //std::cout<<"dz ";
            if (e.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2) {
                fail = true;
                elcutflow->Fill(e.pt(),9);
                eletacutflow->Fill(e.eta(),9);
            }
            //std::cout<<"nhits ";
        }
        else if (fabs(e.superCluster()->eta())<1.556) continue;
        else if (fabs(e.superCluster()->eta())<2.5) {
            if (e.full5x5_sigmaIetaIeta() >= 0.0301) {
                fail = true;
                elcutflow->Fill(e.pt(),2);
                eletacutflow->Fill(e.eta(),2);
            }
            //std::cout<<"sihih ";
            if (fabs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00814) {
                fail = true;
                elcutflow->Fill(e.pt(),3);
                eletacutflow->Fill(e.eta(),3);
            }
            //std::cout<<"deta ";
            if (fabs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.182) {
                fail = true;
                elcutflow->Fill(e.pt(),4);
                eletacutflow->Fill(e.eta(),4);
            }
            //std::cout<<"dphi ";
            if (e.hcalOverEcal() >= 0.0897) {
                fail = true;
                elcutflow->Fill(e.pt(),5);
                eletacutflow->Fill(e.eta(),5);
            }
            //std::cout<<"hoe ";
            if (Ooemoop >= 0.126) {
                fail = true;
                elcutflow->Fill(e.pt(),6);
                eletacutflow->Fill(e.eta(),6);
            }
            //std::cout<<"eop ";
            if (fabs(e.gsfTrack()->dxy(PVtx)) >= 0.118) {
                fail = true;
                elcutflow->Fill(e.pt(),7);
                eletacutflow->Fill(e.eta(),7);
            }
            //std::cout<<"dxy ";
            if (fabs(e.gsfTrack()->dz(PVtx)) >= 0.822) {
                fail = true;
                elcutflow->Fill(e.pt(),8);
                eletacutflow->Fill(e.eta(),8);
            }
            //std::cout<<"dz ";
            if (e.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) {
                fail = true;
                elcutflow->Fill(e.pt(),9);
                eletacutflow->Fill(e.eta(),9);
            }
            //std::cout<<"nhits ";
        }
        else continue;
        if (ConversionTools::hasMatchedConversion(e, conversions, beamspot.position())) {
            fail = true;
            elcutflow->Fill(e.pt(),10);
            eletacutflow->Fill(e.eta(),10);
            }
        //if (fabs(e.convDist()) < 0.02 && fabs(e.convDcot()) < 0.02) fail = true;
        if (fail) continue;
        tmpvec.SetPtEtaPhiE(e.pt(),e.eta(),e.phi(),e.energy());
        rev.push_back(tmpvec);
        if (fabs(e.superCluster()->eta())<1.4442) {
            if (e.full5x5_sigmaIetaIeta() >= 0.0101) continue;
            if (fabs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00926) continue;
            if (fabs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.0336) continue;
            if (e.hcalOverEcal() >= 0.0597) continue;
            if (Ooemoop >= 0.012) continue;
            if (fabs(e.gsfTrack()->dxy(PVtx)) >= 0.0111) continue;
            if (fabs(e.gsfTrack()->dz(PVtx)) >= 0.0466) continue;
            if (e.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2) continue;
        }
        else {
            if (e.full5x5_sigmaIetaIeta() >= 0.0279) continue;
            if (fabs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00724) continue;
            if (fabs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.0918) continue;
            if (e.hcalOverEcal() >= 0.0615) continue;
            if (Ooemoop >= 0.00999) continue;
            if (fabs(e.gsfTrack()->dxy(PVtx)) >= 0.0351) continue;
            if (fabs(e.gsfTrack()->dz(PVtx)) >= 0.417) continue;
            if (e.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) continue;
        }
        retv.push_back(tmpvec);
    }

    edm::Handle<std::vector<reco::GenParticle>> genparts;
    iEvent.getByToken(genToken_, genparts);

    for (const reco::GenParticle &g : *genparts) {
        //std::cout<<"\nGenPart ID"<<g.pdgId()<<" Stat"<<g.status()<<" Pt"<<g.pt();
        if (g.pt() < 10 or g.status()!=1) continue;
        if (abs(g.pdgId())!=11 and abs(g.pdgId())!=13) continue;
        tmpvec.SetPtEtaPhiE(g.pt(),g.eta(),g.phi(),g.energy());
        int pid = abs(g.pdgId());
        bool match = false;
        if (pid==11) {
            elden->Fill(g.pt());
            eletaden->Fill(g.eta());
            eldist->Fill(g.pt(),g.eta());
            for (unsigned int i=0; i<rev.size(); i++) {
                if (tmpvec.DeltaR(rev[i])<0.3 && fabs((g.pt()-rev[i].Pt())/g.pt()) < 0.5 ) {//this is the matching criteria (dR<0.3 and dpT/pT<0.5)
                    match=true;
                    break;
                }
            }
            if (match) {
                elnum->Fill(g.pt());
                eletanum->Fill(g.eta());
            }
        }
    }

    for (const reco::GenParticle &g : *genparts) {
        if (debug) std::cout<<"\nGenPart ID"<<g.pdgId()<<" Stat"<<g.status()<<" Pt"<<g.pt();
        if (g.pt() < 10 or (g.status()!=23 and g.status()!=3)) continue;
        if (abs(g.pdgId())!=11 and abs(g.pdgId())!=13) continue;
        tmpvec.SetPtEtaPhiE(g.pt(),g.eta(),g.phi(),g.energy());
        int pid = abs(g.pdgId());
        bool match = false;
        if (pid==11) {
            elden3->Fill(g.pt());
            eletaden3->Fill(g.eta());
            for (unsigned int i=0; i<rev.size(); i++) {
                if (tmpvec.DeltaR(rev[i])<0.3 && fabs((g.pt()-rev[i].Pt())/g.pt()) < 0.5 ) {
                    match=true;
                    break;
                }
            }
            if (match) {
                elnum3->Fill(g.pt());
                eletanum3->Fill(g.eta());
            }
        }
    }

}

//define this as a plug-in
DEFINE_FWK_MODULE(RECOAnalyzer);
