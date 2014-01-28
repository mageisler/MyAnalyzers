// -*- C++ -*-
//
// Package:    ZDaughtersFinder
// Class:      ZDaughtersFinder
// 
/**\class DaughtersFinder ZDaughtersFinder.cc MGeisler/DaughtersFinder/plugins/ZDaughtersFinder.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Wed Mar 13 10:54:35 CET 2013
// $Id$
//
//

// system include files
#include <memory>
#include <string>

// user include files
#include "MGeisler/DaughtersFinder/interface/ZDaughtersFinder.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidateT.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// constants, enums and typedefs
//

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

typedef math::PtEtaPhiMLorentzVectorD TLV;

//
// static data member definitions
//

//
// constructors and destructor
//

ZDaughtersFinder::ZDaughtersFinder(const edm::ParameterSet& iConfig)
{
   
  //parameters for vs_eta plots
  minEta  = -2.5;  
  maxEta  = 2.5;
  nintEta = 50;

  //parameters for vs_pt plots
  minPt  = 0.1;
  maxPt  = 100.;
  nintPt = 40;
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;

  //configure TP selectors

  using namespace reco::modules;

  ParameterSet generalTpSignalSelectorPSet = iConfig.getParameter<ParameterSet>("generalTpSelector");

  generalTpSignalSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpSignalSelectorPSet));
  
  //now do what ever initialization is needed

  input_recoTrackColls_ = iConfig.getParameter< vector<InputTag> >("recoTrackColls");

  input_recoJetColls_ = iConfig.getParameter< vector<InputTag> >("recoJetColls");

  input_pdgIds_ = iConfig.getParameter<vector<int> >("pdgIds");

  input_genParticles_ = iConfig.getParameter<InputTag>("GenParticles");

  input_genJets_ = iConfig.getParameter<InputTag>("GenJets");

  input_trackingParticles_ = iConfig.getParameter<InputTag>("TrackingParticles");

  input_puSummary_ = iConfig.getParameter<InputTag>("PUInfo");

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  //tracks
  for(unsigned tcl=0; tcl<input_recoTrackColls_.size(); tcl++){

    InputTag input_recoTrackColl_ = input_recoTrackColls_[tcl];
    string dirName = "";
    dirName += input_recoTrackColl_.label();

    subDir->push_back(tfs->mkdir(dirName));

    //book profiles
    p_trackEfficiency.push_back(subDir->at(tcl).make<TProfile3D>("p_trackEfficiency", "efficiency vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
    
    p_trackEfficiency[tcl]->Sumw2();

  }

  //jets
  for(unsigned jcl=0; jcl<input_recoJetColls_.size(); jcl++){
  
    int index = input_recoTrackColls_.size() + jcl;

    InputTag input_recoJetColl_ = input_recoJetColls_[jcl];
    string dirName = "";
    dirName += input_recoJetColl_.label();

    subDir->push_back(tfs->mkdir(dirName));

    //book profiles
    p_jetEfficiency.push_back(subDir->at(index).make<TProfile3D>("p_jetEfficiency", "efficiency vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
    p_jetPtResponse.push_back(subDir->at(index).make<TProfile3D>("p_jetPtResponse", "p_{t}-response vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));

  }

}


ZDaughtersFinder::~ZDaughtersFinder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool
ZDaughtersFinder::isInVector(int ref_pdg, std::vector<int> input_pdgIds_)
{

  for ( unsigned pdg_idx=0; pdg_idx!=input_pdgIds_.size(); pdg_idx++ ) {
  
    if ( ref_pdg==input_pdgIds_[pdg_idx] ) return true;
  
  }
  
  return false;

}

void
ZDaughtersFinder::addDaughters(GenParticle inGenParticle, GenParticleCollection* gen_Zdaughters, std::vector<int> input_pdgIds_)
{

  if ( inGenParticle.numberOfDaughters()>0 ) {
    
    const RefVector<vector<GenParticle> > igpds = inGenParticle.daughterRefVector();  
    RefVector<vector<GenParticle> >::const_iterator d_ite = igpds.begin();
      
    for ( ; d_ite!=igpds.end(); d_ite++ ) {
    
      GenParticle gp_tmp = **d_ite;
    
      if ( isInVector( gp_tmp.pdgId(), input_pdgIds_ ) ) return;
          
    }
      
    for ( d_ite = igpds.begin(); d_ite!=igpds.end(); d_ite++ ) {
    
      GenParticle gp_tmp = **d_ite;
        
      addDaughters(gp_tmp,gen_Zdaughters, input_pdgIds_); 
          
    }
      
  } else {

    if ( !isInVector( inGenParticle.pdgId(), input_pdgIds_ ) ) gen_Zdaughters->push_back( inGenParticle );
      
  }
    
  return;

}

void
ZDaughtersFinder::isGenJet(HepMC::GenParticle tpg, GenJetCollection* inGenJetColl, vector<GenJetRef>* outGenJetColl)
{

  for ( unsigned jet_ite=0; jet_ite<inGenJetColl->size(); jet_ite++ ) {
 
    GenJetRef jetRef(inGenJetColl,jet_ite);

    vector<const GenParticle*> constituents = jetRef->getGenConstituents();

    for ( unsigned const_ite=0; const_ite<constituents.size(); const_ite++ ) {

      if( isGenPart(tpg,*constituents[const_ite]) ){
      
        outGenJetColl->push_back( jetRef );
//         inGenJetColl->erase( inGenJetColl->begin() + jet_ite );
        break;
        
      }

    }
	
  }

  return;

}

bool
ZDaughtersFinder::isGenPart(HepMC::GenParticle tpg, reco::GenParticle gpd)
{

  return( ( fabs(tpg.pdg_id()) == fabs(gpd.pdgId()) ) && 
          ( fabs(tpg.momentum().theta() - gpd.theta())<=0.1 ) && 
          ( fabs(tpg.momentum().phi() - gpd.phi())<=0.1 ) && 
          ( fabs(tpg.momentum().perp() - gpd.pt())<=0.1 ) );

}

JetRef 
ZDaughtersFinder::get_best_matching_recoJet(GenJetRef refJet, Handle<JetCollection> jetsH)
{

  TLV refJetLV( refJet->pt(),
		refJet->eta(),
		refJet->phi(),
		0.0 );

  double minDR = 0.3;
  JetRef bestMatch = *(new JetRef());

  for (unsigned rj_ite=0; rj_ite<jetsH->size(); rj_ite++){

    JetRef rjet(jetsH,rj_ite);

    TLV recoJetLV( rjet->pt(),
	  	   rjet->eta(),
		   rjet->phi(),
		   0.0 );

    double dR = deltaR(refJetLV,recoJetLV);

    if (dR<minDR){

      minDR = dR;
      bestMatch = rjet;

    }
	
  }
	 
  return bestMatch;

}

// ------------ method called for each event  ------------
void
ZDaughtersFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //get the pileup information  
  Handle< vector<PileupSummaryInfo> > puinfoH;
  iEvent.getByLabel(input_puSummary_, puinfoH);
  PileupSummaryInfo puinfo;      
  
  for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      puinfo=(*puinfoH)[puinfo_ite];
      break;
    }
  }

  int npu = puinfo.getPU_NumInteractions();

  //get the gen particles   
  Handle<GenParticleCollection> GPCollectionH;
  iEvent.getByLabel(input_genParticles_, GPCollectionH);
  
  GenParticleCollection gen_Zdaughters;
  
  GenParticleCollection::const_iterator gp_ite = GPCollectionH->begin();
  
  for ( ; gp_ite!=GPCollectionH->end(); gp_ite++ ) {
  
    for ( unsigned pdg_idx=0; pdg_idx!=input_pdgIds_.size(); pdg_idx++ ) {
    
      int input_pdgId_ = input_pdgIds_.at(pdg_idx);
  
      if ( fabs(gp_ite->pdgId())==input_pdgId_ ) {
    
        addDaughters(*gp_ite, &gen_Zdaughters, input_pdgIds_); 	  
  	  
      }
      
    }
	  
  } 

  //get the tracking particles   
  Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByLabel(input_trackingParticles_, TPCollectionH);  
  
  vector<TrackingParticleRef> tp_Zdaughters;
  
  for ( unsigned tp_idx=0; tp_idx<TPCollectionH->size(); tp_idx++ ) {
    
    TrackingParticleRef tpr(TPCollectionH,tp_idx);
  
    TrackingParticle::GenParticleRefVector tp_gens = tpr->genParticle();
    
    TrackingParticle::genp_iterator tpg_ite = tp_gens.begin();
  
    for ( ; tpg_ite!=tp_gens.end(); tpg_ite++ ) {
    
      HepMC::GenParticle tpg = **tpg_ite;
  
      GenParticleCollection::const_iterator gpd_ite = gen_Zdaughters.begin();
  
      for ( ; gpd_ite!=gen_Zdaughters.end(); gpd_ite++ ) {
    
        reco::GenParticle gpd = *gpd_ite;
  
        if ( isGenPart(tpg, gpd) ) {
        
          tp_Zdaughters.push_back( tpr );
        
        }
        
      }
    
    }
  
  }

  //get the needed jets 
  Handle<GenJetCollection> genJetCollH;
  iEvent.getByLabel(input_genJets_, genJetCollH);
  GenJetCollection genJetColl = *(genJetCollH.product());
  
  vector<GenJetRef> ZJets;

  for ( unsigned tp_idx=0; tp_idx<tp_Zdaughters.size(); tp_idx++){
    
    TrackingParticleRef tpr = tp_Zdaughters.at(tp_idx);
  
    TrackingParticle::GenParticleRefVector tp_gens = tpr->genParticle();
    
    TrackingParticle::genp_iterator tpg_ite = tp_gens.begin();
  
    for ( ; tpg_ite!=tp_gens.end(); tpg_ite++ ) {
    
      HepMC::GenParticle tpg = **tpg_ite;
      isGenJet(tpg, &genJetColl, &ZJets);
      
    }
     
  }
    
  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits", theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product(); 
  
  for ( unsigned rtc_idx=0; rtc_idx!=input_recoTrackColls_.size(); rtc_idx++ ) { 

    //get track collection from the event
    Handle<View<Track> >  trackCollectionH;
    iEvent.getByLabel(input_recoTrackColls_[rtc_idx], trackCollectionH);

    SimToRecoCollection simRecColl;
    simRecColl=theTrackAssociator_->associateSimToReco(trackCollectionH,TPCollectionH,&iEvent,&iSetup); 

    for ( unsigned tp_idx=0; tp_idx<tp_Zdaughters.size(); tp_idx++){
    
      TrackingParticleRef tpr = tp_Zdaughters.at(tp_idx);
      
      if((*generalTpSignalSelector)(*tpr)){
      
        double tp_eta = tpr->eta();
        double tp_pt = tpr->pt();
      
        vector<pair<RefToBase<Track>, double> > rt;
      
        int effic = 0;

        if ( simRecColl.find(tpr) != simRecColl.end() ) {
      
 	      rt = (vector<pair<RefToBase<Track>, double> >) simRecColl[tpr];
	      if ( rt.size()!=0 ) {
	        effic = 1;
          }
          
        }  
      
        p_trackEfficiency[rtc_idx]->Fill(tp_eta, tp_pt, npu, effic);
        
      }
             
    }
  
  }
  
  for ( unsigned rjc_idx=0; rjc_idx!=input_recoJetColls_.size(); rjc_idx++ ) {

    //get reco jet collection from the event
    Handle<JetCollection> recoJetCollH;
    iEvent.getByLabel(input_recoJetColls_[rjc_idx],recoJetCollH);
  
    for (unsigned gj_ite=0; gj_ite<ZJets.size(); gj_ite++){

      GenJetRef gjr = ZJets.at(gj_ite);
      
      double gj_eta = gjr->eta();
      double gj_pt = gjr->pt();

      JetRef matchedRecoJet = get_best_matching_recoJet(gjr,recoJetCollH); 
           
      if ( matchedRecoJet.isNull() ) {
      
        p_jetEfficiency[rjc_idx]->Fill(gj_eta, gj_pt, npu, 0.);
        
      } else {
      
        double pt_response = gj_pt*1/matchedRecoJet->pt();
      
        p_jetEfficiency[rjc_idx]->Fill(gj_eta, gj_pt, npu, 1.);
        p_jetPtResponse[rjc_idx]->Fill(gj_eta, gj_pt, npu, pt_response);  
           
      }

    }        
    
  } 
	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZDaughtersFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZDaughtersFinder);
