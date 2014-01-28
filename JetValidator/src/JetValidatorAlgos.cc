#include "MGeisler/JetValidator/interface/JetValidatorAlgos.h"

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
   
using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

//
// constants, enums and typedefs
//

typedef math::PtEtaPhiMLorentzVectorD TLV; 

JetValidatorAlgos::JetValidatorAlgos(const edm::ParameterSet& iConfig)
{
  //parameters for eta plots
  minEta  = -5.;  
  maxEta  = 5.;
  nintEta = 50;

  //parameters for pt plots
  minPt  = 0.1;
  maxPt  = 1000.;
  nintPt = 40;
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;
  
  //parameters for Jet plots
  minJetcount  = -0.5;
  maxJetcount  = 99.5;
  nintJetcount = 100;

  // fix for the LogScale by Ryan
  useLogpt_ = iConfig.getParameter<bool>("UseLogPt");

  if(useLogpt_){
    maxPt=log10(maxPt);
    if(minPt > 0){
      minPt=log10(minPt);
    }else{
      minPt=log10(0.1);
    }
  }

  // fix for the LogScale by Ryan
  MaxNumberOfjetsPerEvent_ = iConfig.getParameter<unsigned>("MaxNumberOfjetsPerEvent");

}


void 
JetValidatorAlgos::CreateIntervalVectors()
{

  // eta vectors

  double eta_step=(maxEta-minEta)/nintEta;
  etaIntervals.push_back(minEta);

  for (int k=1;k<nintEta+1;k++) {

    double d=minEta+k*eta_step;
    etaIntervals.push_back(d);

  }

  // pt vectors

  double pt_step=(maxPt-minPt)/nintPt;
  ptIntervals.push_back(minPt);

  for (int k=1;k<nintPt+1;k++) {

    double d;
    if(useLogpt_){
      d=pow(10,minPt+k*pt_step);
    }else{
      d=minPt+k*pt_step;
    }
    ptIntervals.push_back(d);

  }

  // npu vectors

  int stepVertcount=(maxVertcount-minVertcount)/nintVertcount;
  vertcountIntervals.push_back(minVertcount);

  for (int k=1;k<nintVertcount+1;k++) {

    int d=minVertcount+k*stepVertcount;
    vertcountIntervals.push_back(d);

  }   

}

void 
JetValidatorAlgos::setUpVectors()
{

  // eta vectors
  vector<int> etaIntervalsh;

  for (int k=1;k<nintEta+1;k++) {
    etaIntervalsh.push_back(0);
  }

  reco_jets_eta.push_back(etaIntervalsh);
  matchedReco_jets_eta.push_back(etaIntervalsh);
  gen_jets_eta.push_back(etaIntervalsh);
  matchedGen_jets_eta.push_back(etaIntervalsh);


  // pt vectors
  vector<int> ptIntervalsh;

  for (int k=1;k<nintPt+1;k++) {
    ptIntervalsh.push_back(0);
  }

  reco_jets_pt.push_back(ptIntervalsh);
  matchedReco_jets_pt.push_back(ptIntervalsh);
  gen_jets_pt.push_back(ptIntervalsh);
  matchedGen_jets_pt.push_back(ptIntervalsh);


  // npu vectors
  vector<int> vertcountIntervalsh;

  for (int k=1;k<nintVertcount+1;k++) {
    vertcountIntervalsh.push_back(0);
  }   

  reco_jets_npu.push_back(vertcountIntervalsh);
  matchedReco_jets_npu.push_back(vertcountIntervalsh);
  gen_jets_npu.push_back(vertcountIntervalsh);
  matchedGen_jets_npu.push_back(vertcountIntervalsh);

}

void 
JetValidatorAlgos::BookHistos(TFileDirectory subDir) 
{

  //Book generation related histograms

  h_gen_vertex_num.push_back(subDir.make<TH1F>("gen_vertex_num", "Number of genulated vertices", nintVertcount, minVertcount, maxVertcount));

  h_gen_jets_num.push_back(subDir.make<TH1F>("gen_jets_num", "Number of genulated jets", nintJetcount, minJetcount, maxJetcount));

  h_gen_jets_eta.push_back(subDir.make<TH1F>("gen_jets_eta", "Number of genulated jets vs eta", nintEta, minEta, maxEta));
  h_gen_jets_pt.push_back(subDir.make<TH1F>("gen_jets_pt", "Number of genulated jets vs pt", nintPt, minPt, maxPt));
  h_gen_jets_npu.push_back(subDir.make<TH1F>("gen_jets_npu", "Number of genulated jets vs npu", nintVertcount, minVertcount, maxVertcount));


  //Book reconstruction related histograms

  h_reco_jets_num.push_back(subDir.make<TH1F>("reco_jets_num", "Number of reconstructed jets", nintJetcount, minJetcount, maxJetcount));

  h_reco_jets_eta.push_back(subDir.make<TH1F>("reco_jets_eta", "Number of reconstructed jets vs eta", nintEta, minEta, maxEta));
  h_reco_jets_pt.push_back(subDir.make<TH1F>("reco_jets_pt", "Number of reconstructed jets vs pt", nintPt, minPt, maxPt));
  h_reco_jets_npu.push_back(subDir.make<TH1F>("reco_jets_npu", "Number of reconstructed jets vs npu", nintVertcount, minVertcount, maxVertcount));


  //Book association related histograms

  h_matchedGen_jets_eta.push_back(subDir.make<TH1F>("num_assoc(genToReco)_eta", "Number of associated genulated jets vs eta", nintEta, minEta, maxEta));
  h_matchedGen_jets_pt.push_back(subDir.make<TH1F>("num_assoc(genToReco)_pt", "Number of associated genulated jets vs pt", nintPt, minPt, maxPt));
  h_matchedGen_jets_npu.push_back(subDir.make<TH1F>("num_assoc(genToReco)_npu", "Number of associated genulated jets vs npu", nintVertcount, minVertcount, maxVertcount));

  h_matchedReco_jets_eta.push_back(subDir.make<TH1F>("num_assoc(recoTogen)_eta", "Number of associated reconstructed jets vs eta", nintEta, minEta, maxEta));
  h_matchedReco_jets_pt.push_back(subDir.make<TH1F>("num_assoc(recoTogen)_pt", "Number of associated reconstructed jets vs pt", nintPt, minPt, maxPt));
  h_matchedReco_jets_npu.push_back(subDir.make<TH1F>("num_assoc(recoTogen)_npu", "Number of associated reconstructed jets vs npu", nintVertcount, minVertcount, maxVertcount));


  //Book efficiency and fakerate profiles

  p_efficiency.push_back(subDir.make<TProfile3D>("p_efficiency", "efficiency vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_fakerate.push_back(subDir.make<TProfile3D>("p_fakerate", "fake rate vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));


  //Book resolution, response and deltaR profiles

  p_ptResolution.push_back(subDir.make<TProfile3D>("p_ptResolution", "pt resolution vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_ptResponse.push_back(subDir.make<TProfile3D>("p_ptResponse", "pt response vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));

  p_constResolution.push_back(subDir.make<TProfile3D>("p_constResolution", "number of constituents resolution vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_constResponse.push_back(subDir.make<TProfile3D>("p_constResponse", "number of constituents response vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));

  p_deltaR.push_back(subDir.make<TProfile3D>("p_deltaR", "deltaR vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));


  //Book resolution and response profiles for charged only

  p_ptResolutionCharged.push_back(subDir.make<TProfile3D>("p_ptResolutionCharged", "charged pt resolution vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_ptResponseCharged.push_back(subDir.make<TProfile3D>("p_ptResponseCharged", "charged pt response vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));

  p_constResolutionCharged.push_back(subDir.make<TProfile3D>("p_constResolutionCharged", "number of charged constituents resolution vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_constResponseCharged.push_back(subDir.make<TProfile3D>("p_constResponseCharged", "number of charged constituents response vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));

}

void 
JetValidatorAlgos::fill_independent_histos(int counter, int npu, int rj, int sj)
{

  h_gen_vertex_num.at(counter)->Fill(npu);
  h_gen_jets_num.at(counter)->Fill(sj);

  h_reco_jets_num.at(counter)->Fill(rj);

}

JetRef 
JetValidatorAlgos::get_best_matching_recoJet(GenJetRef refJet, Handle<JetCollection> jetsH)
{

  TLV refJetLV( refJet->pt(),
		refJet->eta(),
		refJet->phi(),
		0.0 );

  double minDR = 0.3;
  JetRef bestMatch;

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

GenJetRef 
JetValidatorAlgos::get_best_matching_genJet(JetRef refJet, Handle<GenJetCollection> jetsH)
{

  TLV refJetLV( refJet->pt(),
		refJet->eta(),
		refJet->phi(),
		0.0 );

  double minDR = 0.3;
  GenJetRef bestMatch;

  for (unsigned gj_ite=0; gj_ite<jetsH->size(); gj_ite++){

    GenJetRef gjet(jetsH,gj_ite);

    TLV genJetLV( gjet->pt(),
  	          gjet->eta(),
		  gjet->phi(),
		  0.0 );

    double dR = deltaR(refJetLV,genJetLV);

    if (dR<minDR){

      minDR = dR;
      bestMatch = gjet;

    }
	
  }
	 
  return bestMatch;

}

GenJet
JetValidatorAlgos::get_charged_gen_jet(GenJetRef genJetRef)
{

  GenJet::LorentzVector chargedLV;
  vector< Ptr<Candidate> > chargedConstituents;
  const GenJet::Specific jetSpecific;

  for (unsigned c_ite=0; c_ite<genJetRef->numberOfDaughters(); c_ite++){

    const GenParticle* consti = genJetRef->getGenConstituent(c_ite);

    if ( consti->charge()!=0 ){

      chargedLV+= consti->p4();
      chargedConstituents.push_back((Ptr<Candidate>) genJetRef->daughterPtr(c_ite));

    }
    
  }

  return GenJet(chargedLV, genJetRef->vertex(), jetSpecific, chargedConstituents);

}

pat::Jet
JetValidatorAlgos::get_charged_reco_jet(JetRef jetRef)
{

  pat::Jet::LorentzVector chargedLV;
  vector< Ptr<Candidate> > chargedConstituents;

  for (unsigned c_ite=0; c_ite<jetRef->numberOfDaughters(); c_ite++){

    PFCandidatePtr consti = jetRef->getPFConstituent(c_ite);

    if ( consti->charge()!=0 ){

      chargedLV+= consti->p4();
      chargedConstituents.push_back((Ptr<Candidate>) jetRef->daughterPtr(c_ite));

    }
    
  }



  return PATObject< reco::Jet >::Jet(chargedLV, jetRef->vertex(), chargedConstituents);

}

void 
JetValidatorAlgos::fill_recoAssociated_genJet_histos(int counter, reco::GenJetRef genJet, JetRef matchedRecoJet, int npu, unsigned jet_number)
{

  bool isMatched = !matchedRecoJet.isNull();

  double genPt = genJet->pt();
  double genEta = genJet->eta();

  // vs eta
  for(unsigned int f=0; f<etaIntervals.size()-1; f++){
    if( genEta>=etaIntervals[f] &&
        genEta<etaIntervals[f+1] ){
      gen_jets_eta[counter][f]++;
      if(isMatched){
	matchedGen_jets_eta[counter][f]++; 
      }
      break;
    }
  } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

  // vs pt
  for(unsigned int f=0; f<ptIntervals.size()-1; f++){
    if( genPt>=ptIntervals[f] &&
        genPt<ptIntervals[f+1] ){
      gen_jets_pt[counter][f]++;
      if(isMatched){
	matchedGen_jets_pt[counter][f]++; 
      }
      break;
    }
  } // END for (unsigned int f=0; f<ptIntervals.size()-1; f++){

  // vs num pileup vertices
  for(unsigned int f=0; f<vertcountIntervals.size()-1; f++){
    if(npu == vertcountIntervals[f]){
      gen_jets_npu[counter][f]++;
      if(isMatched){
        matchedGen_jets_npu[counter][f]++;
      }
      break;
    }    
  }// End for (unsigned int f=0; f<vertcountIntervals.size()-1; f++){

  //efficiency vs eta vs pt vs npu
  double eff = 0.;
  if (isMatched) eff = 1.;

  p_efficiency[counter]->Fill(genEta,genPt,npu,eff);

  if (jet_number>=MaxNumberOfjetsPerEvent_) return;

  //resolution, response and deltaR vs eta vs pt vs npu
  if (isMatched){

    double recoPt = matchedRecoJet->pt();

    double pt_resolution = (recoPt - genPt) / genPt;
    double pt_response = recoPt/genPt;

    p_ptResolution[counter]->Fill(genEta,genPt,npu,pt_resolution);
    p_ptResponse[counter]->Fill(genEta,genPt,npu,pt_response);

    int genConst = genJet->nConstituents();
    int recoConst = matchedRecoJet->nConstituents();

    double const_resolution = (recoConst - genConst) *1./ genConst;
    double const_response = recoConst *1./genConst;

    p_constResolution[counter]->Fill(genEta,genPt,npu,const_resolution);
    p_constResponse[counter]->Fill(genEta,genPt,npu,const_response);

    double delR = deltaR(*genJet,*matchedRecoJet);

    p_deltaR[counter]->Fill(genEta,genPt,npu,delR);

    //charged only

    GenJet genJetCharged = get_charged_gen_jet(genJet);
    pat::Jet matchedRecoJetCharged = get_charged_reco_jet(matchedRecoJet);

    double genPtCharged = genJetCharged.pt();
    double recoPtCharged = matchedRecoJetCharged.pt();

    double pt_resolutionCharged = (recoPtCharged - genPtCharged) / genPtCharged;
    double pt_responseCharged = recoPtCharged/genPtCharged;

    p_ptResolutionCharged[counter]->Fill(genEta,genPtCharged,npu,pt_resolutionCharged);
    p_ptResponseCharged[counter]->Fill(genEta,genPtCharged,npu,pt_responseCharged);

    int genConstCharged = genJetCharged.nConstituents();
    int recoConstCharged = matchedRecoJetCharged.nConstituents();

    double const_resolutionCharged = (recoConstCharged - genConstCharged) *1./ genConstCharged;
    double const_responseCharged = recoConstCharged *1./genConstCharged;

    p_constResolutionCharged[counter]->Fill(genEta,genPtCharged,npu,const_resolutionCharged);
    p_constResponseCharged[counter]->Fill(genEta,genPtCharged,npu,const_responseCharged);

  }

}

void 
JetValidatorAlgos::fill_genAssociated_recoJet_histos(int counter, JetRef recoJet, reco::GenJetRef matchedGenJet, int npu)
{

  bool isMatched = !matchedGenJet.isNull();

  double recoPt = recoJet->pt();
  double recoEta = recoJet->eta();

  // vs eta
  for(unsigned int f=0; f<etaIntervals.size()-1; f++){
    if( recoEta>=etaIntervals[f] &&
        recoEta<etaIntervals[f+1] ){
      reco_jets_eta[counter][f]++;
      if(isMatched){
	matchedReco_jets_eta[counter][f]++; 
      }
      break;
    }
  } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

  // vs pt
  for(unsigned int f=0; f<ptIntervals.size()-1; f++){
    if( recoPt>=ptIntervals[f] &&
        recoPt<ptIntervals[f+1] ){
      reco_jets_pt[counter][f]++;
      if(isMatched){
	matchedReco_jets_pt[counter][f]++; 
      }
      break;
    }
  } // END for (unsigned int f=0; f<ptIntervals.size()-1; f++){

  // vs num pileup vertices
  for(unsigned int f=0; f<vertcountIntervals.size()-1; f++){
    if(npu == vertcountIntervals[f]){
      reco_jets_npu[counter][f]++;
      if(isMatched){
        matchedReco_jets_npu[counter][f]++;
      }
      break;
    }    
  }// End for (unsigned int f=0; f<vertcountIntervals.size()-1; f++){

  //fake rate vs eta vs pt vs npu
  double fr = 1.;
  if (isMatched) fr = 0.;

  p_fakerate[counter]->Fill(recoEta,recoPt,npu,fr);

}

void 
JetValidatorAlgos::fillHistosFromVectors(int counter)
{

  fillPlotFromVector(h_gen_jets_eta[counter],gen_jets_eta[counter]);
  fillPlotFromVector(h_gen_jets_pt[counter],gen_jets_pt[counter]);
  fillPlotFromVector(h_gen_jets_npu[counter],gen_jets_npu[counter]);
  fillPlotFromVector(h_matchedGen_jets_eta[counter],matchedGen_jets_eta[counter]);
  fillPlotFromVector(h_matchedGen_jets_pt[counter],matchedGen_jets_pt[counter]);
  fillPlotFromVector(h_matchedGen_jets_npu[counter],matchedGen_jets_npu[counter]);

  fillPlotFromVector(h_reco_jets_eta[counter],reco_jets_eta[counter]);
  fillPlotFromVector(h_reco_jets_pt[counter],reco_jets_pt[counter]);
  fillPlotFromVector(h_reco_jets_npu[counter],reco_jets_npu[counter]);
  fillPlotFromVector(h_matchedReco_jets_eta[counter],matchedReco_jets_eta[counter]);
  fillPlotFromVector(h_matchedReco_jets_pt[counter],matchedReco_jets_pt[counter]);
  fillPlotFromVector(h_matchedReco_jets_npu[counter],matchedReco_jets_npu[counter]);

}

void
JetValidatorAlgos::fillPlotFromVector(TH1F* histo,vector<int> vIn)
{

  for(unsigned ite=0; ite<vIn.size();ite++){
    histo->SetBinContent(ite+1,vIn.at(ite));
  } 

}