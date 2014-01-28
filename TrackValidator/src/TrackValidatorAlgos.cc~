#include "MGeisler/TrackValidator/interface/TrackValidatorAlgos.h"

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

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

// ROOT include files
#include <TH1F.h>
#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

TrackValidatorAlgos::TrackValidatorAlgos(const edm::ParameterSet& iConfig)
{
  //parameters for vs_eta plots
  minEta  = -2.5;  
  maxEta  = 2.5;
  nintEta = 50;
  
  //parameters for vs_etaErr plots
  minEtaErr  = -0.2;  
  maxEtaErr  = 0.2;
  nintEtaErr = 40;

  //parameters for vs_pt plots
  minPt  = 0.1;
  maxPt  = 100.;
  nintPt = 40;

  //parameters for vs_ptErr plots
  minPtErr  = 0.1;
  maxPtErr  = 100.;
  nintPtErr = 40;

  //parameters for vs_ptErrRel plots
  minPtErrRel  = 0.;
  maxPtErrRel  = 2.;
  nintPtErrRel = 40;
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;
  
  //parameters for track number plots
  minTrackcount  = -0.5;
  maxTrackcount  = 999.5;
  nintTrackcount = 1000;
  
  //parameters for track number plots
  minTW  = -0.025;
  maxTW = 1.025;
  nintTW = 21;
  
  //parameters for p contribution for jets plots
  minContribution  = 0.;
  maxContribution  = 1.;
  nintContribution = 200;

  //configure TP selectors

  using namespace reco::modules;

  ParameterSet generalTpSignalSelectorPSet = iConfig.getParameter<ParameterSet>("generalTpSelector");
  

  generalTpSignalSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpSignalSelectorPSet));

  ParameterSet generalTpPUSelectorPSet = generalTpSignalSelectorPSet;
  Entry sOname("signalOnly",false,true);
  generalTpPUSelectorPSet.insert(true,"signalOnly",sOname);

  generalTpPUSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpPUSelectorPSet));
  

  ParameterSet generalTpSignalSelectorPtPSet = generalTpSignalSelectorPSet;
  Entry ptOname("ptMin",0.1,true);
  generalTpSignalSelectorPtPSet.insert(true,"ptMin",ptOname);

  generalTpSignalSelectorPt = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpSignalSelectorPtPSet));

  // fix for the LogScale by Ryan
  useLogpt_ = iConfig.getParameter<bool>("UseLogPt");

  useJetWeighting_ = iConfig.getUntrackedParameter<bool>("useJetWeighting",false);
  genJetCollLabel_ = iConfig.getUntrackedParameter<string>("genJetCollLabel");
  jetCollLabel_ = iConfig.getUntrackedParameter<string>("jetCollLabel");

  if(useLogpt_){
    maxPt = log10(maxPt);
    maxPtErr = log10(maxPtErr);
    if ( minPt > 0 ) {
      minPt = log10(minPt);
    } else {
      minPt = log10(0.1);
    }
    if ( minPtErr > 0 ) {
      minPtErr = log10(minPtErr);
    } else {
      minPtErr = log10(0.1);
    }
  }

}


void 
TrackValidatorAlgos::CreateIntervalVectors()
{

  // eta vectors

  double eta_step=(maxEta-minEta)/nintEta;
  etaintervals.push_back(minEta);

  for (int k=1;k<nintEta+1;k++) {

    double d=minEta+k*eta_step;
    etaintervals.push_back(d);

  }

  // pt vectors

  double pt_step=(maxPt-minPt)/nintPt;
  ptintervals.push_back(minPt);

  for (int k=1;k<nintPt+1;k++) {

    double d;
    if(useLogpt_){
      d=pow(10,minPt+k*pt_step);
    }else{
      d=minPt+k*pt_step;
    }
    ptintervals.push_back(d);

  }

  // npu vectors

  int stepVertcount= 1;
  vertcountintervals.push_back(0);

  for (int k=1;k<nintVertcount+1;k++) {

    int d=minVertcount+k*stepVertcount;
    vertcountintervals.push_back(d);

  }   

}

void 
TrackValidatorAlgos::GetInputCollections(const Event& iEvent){

  //get jet collection from the event
  if(useJetWeighting_){
    iEvent.getByLabel(genJetCollLabel_,genJetCollH);
    iEvent.getByLabel(jetCollLabel_,jetCollH);
  }

}

void 
TrackValidatorAlgos::setUpVectors()
{

  // eta vectors
  vector<double> etaintervalsh;

  for (int k=1;k<nintEta+1;k++) {
    etaintervalsh.push_back(0.);
  }

  allSignalTP_eta.push_back(etaintervalsh);
  allRT_eta.push_back(etaintervalsh);
  assSignalTP_eta.push_back(etaintervalsh);
  assSignalRT_eta.push_back(etaintervalsh);
  allSigRT_eta.push_back(etaintervalsh);
  allPURT_eta.push_back(etaintervalsh);
  allAssPURT_eta.push_back(etaintervalsh);
  allRemovedRT_eta.push_back(etaintervalsh);
  removedSigRT_eta.push_back(etaintervalsh);
  removedPURT_eta.push_back(etaintervalsh);

  // pt vectors
  vector<double> ptintervalsh;

  for (int k=1;k<nintPt+1;k++) {
    ptintervalsh.push_back(0.);
  }

  allSignalTP_pt.push_back(ptintervalsh);
  allRT_pt.push_back(ptintervalsh);
  assSignalTP_pt.push_back(ptintervalsh);
  assSignalRT_pt.push_back(ptintervalsh);
  allSigRT_pt.push_back(ptintervalsh);
  allPURT_pt.push_back(ptintervalsh);
  allRemovedRT_pt.push_back(ptintervalsh);
  removedSigRT_pt.push_back(ptintervalsh);
  removedPURT_pt.push_back(ptintervalsh);

  // npu vectors
  vector<double> vertcountintervalsh;

  for (int k=1;k<nintVertcount+1;k++) {
    vertcountintervalsh.push_back(0.);
  }   

  allSignalTP_npu.push_back(vertcountintervalsh);
  allRT_npu.push_back(vertcountintervalsh);
  assSignalTP_npu.push_back(vertcountintervalsh);
  assSignalRT_npu.push_back(vertcountintervalsh);
  allPURT_npu.push_back(vertcountintervalsh);
  allSigRT_npu.push_back(vertcountintervalsh);
  allRemovedRT_npu.push_back(vertcountintervalsh);
  removedSigRT_npu.push_back(vertcountintervalsh);
  removedPURT_npu.push_back(vertcountintervalsh);


  //counting vectors

  sim_tracks.push_back(0);

}

void 
TrackValidatorAlgos::BookHistos(TFileDirectory subDir) 
{

  effic_npu_Contr.push_back(subDir.make<TH2F>("effic_npu_Contr", "efficiency vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  num_simul_tracks_npu_Contr.push_back(subDir.make<TH2F>("num_simul_tracks_npu_Contr", "Number of simulated tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));
  num_assoc_npu_Contr.push_back(subDir.make<TH2F>("num_assoc(simToReco)_npu_Contr", "Number of associated simulated tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  fakerate_npu_Contr.push_back(subDir.make<TH2F>("fakerate_npu_Contr", "fakerate vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));
  fakerate_npu_Contr_help.push_back(subDir.make<TH2F>("fakerate_npu_Contr_help", "HELP vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  num_reco_tracks_npu_Contr.push_back(subDir.make<TH2F>("num_reco_tracks_npu_Contr", "Number of reconstructed tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));
  num_assoc2_npu_Contr.push_back(subDir.make<TH2F>("num_assoc(recoToSim)_npu_Contr", "Number of associated  reco tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  weights.push_back(subDir.make<TH1F>("weights", "weights", 1000, 0., 2.));

  h_tweight.push_back(subDir.make<TH1F>("h_tweight", "Number of tracks vs. track weight; track weight; number of tracks", nintTW,minTW,maxTW));

  //Book PileUp related histograms

  PU_effic_eta.push_back(subDir.make<TH1F>("PU_effic_eta", "PU_effic vs eta", nintEta, minEta, maxEta));
  PU_effic_pt.push_back(subDir.make<TH1F>("PU_effic_pt", "PU_effic vs pt", nintPt, minPt, maxPt));
  PU_effic_npu.push_back(subDir.make<TH1F>("PU_effic_npu", "PU_effic vs npu", nintVertcount, minVertcount, maxVertcount));

  PU_fakerate_1_eta.push_back(subDir.make<TH1F>("PU_fakerate_1_eta", "PU_fakerate 1 vs eta", nintEta, minEta, maxEta));
  PU_fakerate_1_pt.push_back(subDir.make<TH1F>("PU_fakerate_1_pt", "PU_fakerate 1 vs pt", nintPt, minPt, maxPt));
  PU_fakerate_1_npu.push_back(subDir.make<TH1F>("PU_fakerate_1_npu", "PU_fakerate 1 vs npu", nintVertcount, minVertcount, maxVertcount));

  PU_fakerate_2_eta.push_back(subDir.make<TH1F>("PU_fakerate_2_eta", "PU_fakerate 2 vs eta", nintEta, minEta, maxEta));
  PU_fakerate_2_pt.push_back(subDir.make<TH1F>("PU_fakerate_2_pt", "PU_fakerate 2 vs pt", nintPt, minPt, maxPt));
  PU_fakerate_2_npu.push_back(subDir.make<TH1F>("PU_fakerate_2_npu", "PU_fakerate 2 vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book efficiency and fakerate histograms

  effic_eta.push_back(subDir.make<TH1F>("effic_eta", "effic vs eta", nintEta, minEta, maxEta));
  effic_pt.push_back(subDir.make<TH1F>("effic_pt", "effic vs pt", nintPt, minPt, maxPt));
  effic_npu.push_back(subDir.make<TH1F>("effic_npu", "effic vs npu", nintVertcount, minVertcount, maxVertcount));

  fakerate_eta.push_back(subDir.make<TH1F>("fakerate_eta", "fakerate vs eta", nintEta, minEta, maxEta));
  fakerate_pt.push_back(subDir.make<TH1F>("fakerate_pt", "fakerate vs pt", nintPt, minPt, maxPt));
  fakerate_npu.push_back(subDir.make<TH1F>("fakerate_npu", "fakerate vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book simulation related histograms

  num_simul_tracks.push_back(subDir.make<TH1F>("num_simul_tracks", "Number of simulated tracks", nintTrackcount, minTrackcount, maxTrackcount));

  num_track_simul_eta.push_back(subDir.make<TH1F>("num_track_simul_eta", "Number of simulated tracks vs eta", nintEta, minEta, maxEta));
  num_track_simul_pt.push_back(subDir.make<TH1F>("num_track_simul_pt", "Number of simulated tracks vs pt", nintPt, minPt, maxPt));
  num_track_simul_npu.push_back(subDir.make<TH1F>("num_track_simul_npu", "Number of simulated tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  num_simul_vertex.push_back(subDir.make<TH1F>("num_simul_vertex", "Number of simulated vertices", nintVertcount, minVertcount, maxVertcount));

  //Book reconstruction related histograms

  num_reco_tracks.push_back(subDir.make<TH1F>("num_reco_tracks", "Number of reconstructed tracks", nintTrackcount, minTrackcount, maxTrackcount));

  num_track_reco_eta.push_back(subDir.make<TH1F>("num_track_reco_eta", "Number of reconstructed tracks vs eta", nintEta, minEta, maxEta));
  num_track_reco_pt.push_back(subDir.make<TH1F>("num_track_reco_pt", "Number of reconstructed tracks vs pt", nintPt, minPt, maxPt));
  num_track_reco_npu.push_back(subDir.make<TH1F>("num_track_reco_npu", "Number of reconstructed tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  num_removed_reco_signal_eta.push_back(subDir.make<TH1F>("num_removed_reco_signal_eta", "Number of removed reconstructed signal tracks vs eta", nintEta, minEta, maxEta));
  num_removed_reco_eta.push_back(subDir.make<TH1F>("num_removed_reco_eta", "Number of removed reconstructed tracks vs eta", nintEta, minEta, maxEta));

  num_removed_reco_PU_eta.push_back(subDir.make<TH1F>("num_removed_reco_PU_eta", "Number of removed reconstructed pileup tracks vs eta", nintEta, minEta, maxEta));
  num_reco_PU_eta.push_back(subDir.make<TH1F>("num_reco_PU_eta", "Number of reconstructed pileup tracks vs eta", nintEta, minEta, maxEta));

  num_removed_reco_signal_pt.push_back(subDir.make<TH1F>("num_removed_reco_signal_pt", "Number of removed reconstructed signal tracks vs pt", nintPt, minPt, maxPt));
  num_removed_reco_pt.push_back(subDir.make<TH1F>("num_removed_reco_pt", "Number of removed reconstructed tracks vs pt", nintPt, minPt, maxPt));
  num_removed_reco_PU_pt.push_back(subDir.make<TH1F>("num_removed_reco_PU_pt", "Number of removed reconstructed pileup tracks vs pt", nintPt, minPt, maxPt));

  num_removed_reco_signal_npu.push_back(subDir.make<TH1F>("num_removed_reco_signal_npu", "Number of removed reconstructed signal tracks vs npu", nintVertcount, minVertcount, maxVertcount));
  num_removed_reco_npu.push_back(subDir.make<TH1F>("num_removed_reco_npu", "Number of removed reconstructed tracks vs npu", nintVertcount, minVertcount, maxVertcount));
  num_removed_reco_PU_npu.push_back(subDir.make<TH1F>("num_removed_reco_PU_npu", "Number of removed reconstructed pileup tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  num_track_reco_PU_eta.push_back(subDir.make<TH1F>("num_track_reco_PU_eta", "Number of reco tracks from PileUp vs eta", nintEta, minEta, maxEta));
  num_track_reco_signal_eta.push_back(subDir.make<TH1F>("num_track_reco_signal_eta", "Number of reco tracks from Signal vs eta", nintEta, minEta, maxEta));

  num_track_reco_PU_pt.push_back(subDir.make<TH1F>("num_track_reco_PU_pt", "Number of reco tracks from PileUp vs pt", nintPt, minPt, maxPt));
  num_track_reco_signal_pt.push_back(subDir.make<TH1F>("num_track_reco_signal_pt", "Number of reco tracks from Signal vs pt", nintPt, minPt, maxPt));

  num_track_reco_PU_npu.push_back(subDir.make<TH1F>("num_track_reco_PU_npu", "Number of reco tracks from PileUp vs npu", nintVertcount, minVertcount, maxVertcount));
  num_track_reco_signal_npu.push_back(subDir.make<TH1F>("num_track_reco_signal_npu", "Number of reco tracks from Signal vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book association related histograms

  num_assoc_eta.push_back(subDir.make<TH1F>("num_assoc(simToReco)_eta", "Number of associated simulated tracks vs eta", nintEta, minEta, maxEta));
  num_assoc2_eta.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_eta", "Number of associated reconstructed tracks vs eta", nintEta, minEta, maxEta));

  num_assoc_pt.push_back(subDir.make<TH1F>("num_assoc(simToReco)_pt", "Number of associated simulated tracks vs pt", nintPt, minPt, maxPt));
  num_assoc2_pt.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_pt", "Number of associated reconstructed tracks vs pt", nintPt, minPt, maxPt));

  num_assoc_npu.push_back(subDir.make<TH1F>("num_assoc(simToReco)_npu", "Number of associated simulated tracks vs npu", nintVertcount, minVertcount, maxVertcount));
  num_assoc2_npu.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_npu", "Number of associated reconstructed tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book profiles

  p_efficiency_val.push_back(subDir.make<TProfile3D>("p_efficiency_val", "efficiency vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_fakerate_val.push_back(subDir.make<TProfile3D>("p_fakerate_val", "fake rate vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));

  p_efficiencyPileUp_val.push_back(subDir.make<TProfile3D>("p_efficiencyPileUp_val", "PileUp efficiency vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  p_fakeratePileUp_val.push_back(subDir.make<TProfile3D>("p_fakeratePileUp_val", "PileUp fake rate vs eta vs pt vs npu", nintEta, minEta, maxEta, nintPt, minPt, maxPt, nintVertcount, minVertcount, maxVertcount));
  
  
  p_fakerate_err.push_back(subDir.make<TProfile3D>("p_fakerate_err", "fake rate vs #sigma_{eta} vs #sigma_{pt} vs npu", nintEtaErr, minEtaErr, maxEtaErr, nintPtErr, minPtErr, maxPtErr, nintVertcount, minVertcount, maxVertcount));

  p_efficiencyPileUp_err.push_back(subDir.make<TProfile3D>("p_efficiencyPileUp_err", "PileUp efficiency vs #sigma_{eta} vs #sigma_{pt} vs npu", nintEtaErr, minEtaErr, maxEtaErr, nintPtErr, minPtErr, maxPtErr, nintVertcount, minVertcount, maxVertcount));
  p_fakeratePileUp_err.push_back(subDir.make<TProfile3D>("p_fakeratePileUp_err", "PileUp fake rate vs #sigma_{eta} vs #sigma_{pt} vs npu", nintEtaErr, minEtaErr, maxEtaErr, nintPtErr, minPtErr, maxPtErr, nintVertcount, minVertcount, maxVertcount));
  

  p_fakerate_relerr.push_back(subDir.make<TProfile3D>("p_fakerate_relerr", "fake rate vs #sigma_{eta} vs rel. #sigma_{pt} vs npu", nintEtaErr, minEtaErr, maxEtaErr, nintPtErrRel, minPtErrRel, maxPtErrRel, nintVertcount, minVertcount, maxVertcount));

  p_efficiencyPileUp_relerr.push_back(subDir.make<TProfile3D>("p_efficiencyPileUp_relerr", "PileUp efficiency vs #sigma_{eta} vs rel. #sigma_{pt} vs npu", nintEtaErr, minEtaErr, maxEtaErr, nintPtErrRel, minPtErrRel, maxPtErrRel, nintVertcount, minVertcount, maxVertcount));
  p_fakeratePileUp_relerr.push_back(subDir.make<TProfile3D>("p_fakeratePileUp_relerr", "PileUp fake rate vs #sigma_{eta} vs rel. #sigma_{pt} vs npu", nintEtaErr, minEtaErr, maxEtaErr, nintPtErrRel, minPtErrRel, maxPtErrRel, nintVertcount, minVertcount, maxVertcount));

  p_purity_tw.push_back(subDir.make<TProfile>("p_purity_tw", "Purity vs track weight; track weight; purity", nintTW,minTW,maxTW));

}

void 
TrackValidatorAlgos::fill_independent_histos(int counter, int npu, int rt, unsigned st)
{

  num_simul_vertex.at(counter)->Fill(npu);
  num_simul_tracks.at(counter)->Fill(st);
  num_reco_tracks.at(counter)->Fill(rt);

}

void
TrackValidatorAlgos::fill_recoAssociated_simTrack_histos(int counter, TrackingParticle* tp, const Track* track, int npu, unsigned* st)
{

  bool isMatched = track;  
  double tp_eta = tp->momentum().eta();  
  double tp_pt = sqrt(tp->momentum().perp2());

  double weight = 1.;
  if(useJetWeighting_) weight = getTrackWeight(tp,&(*genJetCollH));

  if((*generalTpSignalSelector)(*tp)){

    (*st)++;

    num_simul_tracks_npu_Contr[counter]->Fill(npu,weight);
    if(isMatched) num_assoc_npu_Contr[counter]->Fill(npu,weight);

    //effic vs eta
    for(unsigned int f=0; f<etaintervals.size()-1; f++){
      if(tp_eta>etaintervals[f]&&
	 tp_eta<etaintervals[f+1]){
	allSignalTP_eta[counter][f]+=weight;
	if(isMatched){
	  assSignalTP_eta[counter][f]+=weight;
	}
        break;
      }
    } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

    //effic vs num pileup vertices
    for(unsigned int f=0; f<vertcountintervals.size()-1; f++){
      if(npu>=vertcountintervals[f]&&
         npu<vertcountintervals[f+1]){
        allSignalTP_npu[counter][f]+=weight;
        if(isMatched){
          assSignalTP_npu[counter][f]+=weight;
        }
        break;
      }    
    }// End for (unsigned int f=0; f<vertcountintervals.size()-1; f++){

    //efficiency vs eta vs pt vs npu
    int eff = 0;
    if (isMatched) eff = 1;
  
    p_efficiency_val[counter]->Fill(tp_eta, tp_pt, npu, eff, weight);

  }

  if((*generalTpSignalSelectorPt)(*tp)){

    //effic vs pt
    for(unsigned int f=0; f<ptintervals.size()-1; f++){
      if (sqrt(tp->momentum().perp2())>ptintervals[f]&&
	      sqrt(tp->momentum().perp2())<ptintervals[f+1]){
        allSignalTP_pt[counter][f]+=weight; 
        if ( isMatched ) {
	      assSignalTP_pt[counter][f]+=weight;
        }	
        break;      
      }
    } // End for (unsigned int f=0; f<ptintervals.size()-1; f++){
  
  }


}

void 
TrackValidatorAlgos::fill_simAssociated_recoTrack_histos(int counter, const Track& track, bool isMatched, bool isSigMatched, int npu, TpDoubV tp)
{

  double weight = 1.;      
  if(useJetWeighting_) weight = getTrackWeightReco(track,&(*jetCollH));
  
  double trk_eta = track.eta();
  double trk_etaerr = track.etaError();
  
  double trk_pt = sqrt(track.momentum().perp2());
  double trk_pterr = track.ptError();
  double trk_ptrelerr = trk_pterr/trk_pt;

  num_reco_tracks_npu_Contr[counter]->Fill(npu,weight);
  if(isSigMatched) num_assoc2_npu_Contr[counter]->Fill(npu,weight);

  //fake rate vs eta
  for (unsigned int f=0; f<etaintervals.size()-1; f++){
    if (track.eta()>etaintervals[f]&&
        track.eta()<etaintervals[f+1]) {
      allRT_eta[counter][f]+=weight;
      if (isSigMatched){
	assSignalRT_eta[counter][f]+=weight;
      }else{
        if(isMatched) allAssPURT_eta[counter][f]+=weight;
      }
      break;
    }
  } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

  //fake rate vs pt
  for (unsigned int f=0; f<ptintervals.size()-1; f++){
    if (sqrt(track.momentum().perp2())>ptintervals[f]&&
        sqrt(track.momentum().perp2())<ptintervals[f+1]) {
      allRT_pt[counter][f]+=weight;
      if (isSigMatched){
	assSignalRT_pt[counter][f]+=weight;
      }	  
      break;    
    }
  } // End for (unsigned int f=0; f<ptintervals.size()-1; f++){

  //fake rate vs num pileup vertices
  for (unsigned int f=0; f<vertcountintervals.size()-1; f++){
    if(npu>=vertcountintervals[f]&&
       npu<vertcountintervals[f+1]){
      allRT_npu[counter][f]+=weight;
      if (isSigMatched){
        assSignalRT_npu[counter][f]+=weight;
      }
      break;
    }    
  }// End for (unsigned int f=0; f<vertcountintervals.size()-1; f++){

  //fake rate vs eta vs pt vs npu
  double fr = 1.;
  if (isSigMatched) fr = 0.;

  p_fakerate_val[counter]->Fill(trk_eta, trk_pt, npu, fr, weight);
  p_fakerate_err[counter]->Fill(trk_etaerr, trk_pterr, npu, fr, weight);
  p_fakerate_relerr[counter]->Fill(trk_etaerr, trk_ptrelerr, npu, fr, weight);


}

void 
TrackValidatorAlgos::fill_removedRecoTrack_histos(int counter, const Track& refTrack, bool isSigMatched,  bool isRemoved, int npu, TpDoubV tp)
{

  double weight = 1.;
  if(useJetWeighting_){
    weight = getTrackWeightReco(refTrack,&(*jetCollH));
    weights[counter]->Fill(weight);
  }
  
  double trk_eta = refTrack.eta();
  double trk_etaerr = refTrack.etaError();
  
  double trk_pt = sqrt(refTrack.momentum().perp2());
  double trk_pterr = refTrack.ptError();
  double trk_ptrelerr = trk_pterr/trk_pt;

  // vs eta
  for (unsigned int f=0; f<etaintervals.size()-1; f++){
    if (refTrack.eta()>etaintervals[f]&&
        refTrack.eta()<etaintervals[f+1]) {
      if(isSigMatched){
        allSigRT_eta[counter][f]+=weight; 
       }else{
        allPURT_eta[counter][f]+=weight; 
      }   
      if(isRemoved){
        allRemovedRT_eta[counter][f]+=weight;
        if(isSigMatched){
          removedSigRT_eta[counter][f]+=weight;
        }else{
          removedPURT_eta[counter][f]+=weight; 
        }
      } // END if(isRemoved){
    } // END if(refTrack.eta()>etaintervals[counter][f]&&refTrack.eta()<etaintervals[counter][f+1){
  } // END for(unsigned int f=0; f<etaintervals.size()-1; f++){

  // vs pt
  for (unsigned int f=0; f<ptintervals.size()-1; f++){
    if (sqrt(refTrack.momentum().perp2())>ptintervals[f]&&
        sqrt(refTrack.momentum().perp2())<ptintervals[f+1]){
      if(isSigMatched){
        allSigRT_pt[counter][f]+=weight; 
      }else{
        allPURT_pt[counter][f]+=weight; 
      }   
      if(isRemoved){
        allRemovedRT_pt[counter][f]+=weight;
        if(isSigMatched){
          removedSigRT_pt[counter][f]+=weight; 
        }else{
	  removedPURT_pt[counter][f]+=weight; 
        }
      } // END if(isRemoved){
    } // END if(sqrt(refTrack.momentum().perp2())>ptintervals[counter][f]&&sqrt(refTrack.momentum().perp2())<ptintervals[counter][f+1]){
  } // END for(unsigned int f=0; f<ptintervals.size()-1; f++){

  // vs npu
  for (unsigned int f=0; f<vertcountintervals.size()-1; f++){
    if(npu>=vertcountintervals[f]&&
       npu<vertcountintervals[f+1]){
      if(isSigMatched){
        allSigRT_npu[counter][f]+=weight; 
      }else{
        allPURT_npu[counter][f]+=weight; 
      }   
      if(isRemoved){
        allRemovedRT_npu[counter][f]+=weight;
        if(isSigMatched){
          removedSigRT_npu[counter][f]+=weight;
        }else{
  	  removedPURT_npu[counter][f]+=weight; 
        }
      } // END if(isRemoved){
    } // END if(npu == vertcountintervals[counter][f]){
  } // END for(unsigned int f=0; f<vertcountintervals.size()-1; f++){

  //PileUp vs eta vs pt vs npu
  
  if ( isSigMatched ){

    //fake rate
    double pu_fr  = 1.;
    if (isRemoved) pu_fr = 0.;

    p_fakeratePileUp_val[counter]->Fill(trk_eta, trk_pt, npu, pu_fr, weight);
    p_fakeratePileUp_err[counter]->Fill(trk_etaerr, trk_pterr, npu, pu_fr, weight);
    p_fakeratePileUp_relerr[counter]->Fill(trk_etaerr, trk_ptrelerr, npu, pu_fr, weight);

  } else {

    //efficiency
    double pu_eff = 0.;
    if (isRemoved) pu_eff = 1.;

    p_efficiencyPileUp_val[counter]->Fill(trk_eta, trk_pt, npu, pu_eff, weight);
    p_efficiencyPileUp_err[counter]->Fill(trk_etaerr, trk_pterr, npu, pu_eff, weight);
    p_efficiencyPileUp_relerr[counter]->Fill(trk_etaerr, trk_ptrelerr, npu, pu_eff, weight);

  }    
   
}

void
TrackValidatorAlgos::fillFractionHistosFromVectors(int counter)
{

  fillFractionHisto(effic_eta[counter],assSignalTP_eta[counter],allSignalTP_eta[counter],"effic");
  fillFractionHisto(effic_pt[counter],assSignalTP_pt[counter],allSignalTP_pt[counter],"effic");
  fillFractionHisto(effic_npu[counter],assSignalTP_npu[counter],allSignalTP_npu[counter],"effic");

  fillFractionHisto(fakerate_eta[counter],assSignalRT_eta[counter],allRT_eta[counter],"fakerate");
  fillFractionHisto(fakerate_pt[counter],assSignalRT_pt[counter],allRT_pt[counter],"fakerate");
  fillFractionHisto(fakerate_npu[counter],assSignalRT_npu[counter],allRT_npu[counter],"fakerate");

  fillFractionHisto(PU_effic_eta[counter],removedPURT_eta[counter],allPURT_eta[counter],"effic");
  fillFractionHisto(PU_effic_pt[counter],removedPURT_pt[counter],allPURT_pt[counter],"effic");
  fillFractionHisto(PU_effic_npu[counter],removedPURT_npu[counter],allPURT_npu[counter],"effic");

  fillFractionHisto(PU_fakerate_1_eta[counter],removedSigRT_eta[counter],allRemovedRT_eta[counter],"effic");
  fillFractionHisto(PU_fakerate_1_pt[counter],removedSigRT_pt[counter],allRemovedRT_pt[counter],"effic");
  fillFractionHisto(PU_fakerate_1_npu[counter],removedSigRT_npu[counter],allRemovedRT_npu[counter],"effic");

  fillFractionHisto(PU_fakerate_2_eta[counter],removedSigRT_eta[counter],allSigRT_eta[counter],"effic");
  fillFractionHisto(PU_fakerate_2_pt[counter],removedSigRT_pt[counter],allSigRT_pt[counter],"effic");
  fillFractionHisto(PU_fakerate_2_npu[counter],removedSigRT_npu[counter],allSigRT_npu[counter],"effic");


  effic_npu_Contr[counter]->Divide(num_assoc_npu_Contr[counter],num_simul_tracks_npu_Contr[counter]);

  fakerate_npu_Contr_help[counter]->Divide(num_assoc2_npu_Contr[counter],num_reco_tracks_npu_Contr[counter]);

  for(int x_ite=1; x_ite<=nintVertcount; x_ite++){
    for(int y_ite=1; y_ite<=nintContribution; y_ite++){
      fakerate_npu_Contr[counter]->SetBinContent(x_ite,y_ite, 1. - fakerate_npu_Contr_help[counter]->GetBinContent(x_ite,y_ite));
      fakerate_npu_Contr[counter]->SetBinError(x_ite,y_ite,fakerate_npu_Contr_help[counter]->GetBinError(x_ite,y_ite));
    }
  }

}

void
TrackValidatorAlgos::fillHistosFromVectors(int counter)
{

  fillPlotFromVector(num_track_simul_eta[counter],allSignalTP_eta[counter]);
  fillPlotFromVector(num_assoc_eta[counter],assSignalTP_eta[counter]);

  fillPlotFromVector(num_track_simul_pt[counter],allSignalTP_pt[counter]);
  fillPlotFromVector(num_assoc_pt[counter],assSignalTP_pt[counter]);

  fillPlotFromVector(num_track_simul_npu[counter],allSignalTP_npu[counter]);
  fillPlotFromVector(num_assoc_npu[counter],assSignalTP_npu[counter]);


  fillPlotFromVector(num_track_reco_eta[counter],allRT_eta[counter]);
  fillPlotFromVector(num_assoc2_eta[counter],assSignalRT_eta[counter]);

  fillPlotFromVector(num_track_reco_pt[counter],allRT_pt[counter]);
  fillPlotFromVector(num_assoc2_pt[counter],assSignalRT_pt[counter]);

  fillPlotFromVector(num_track_reco_npu[counter],allRT_npu[counter]);
  fillPlotFromVector(num_assoc2_npu[counter],assSignalRT_npu[counter]);


  fillPlotFromVector(num_removed_reco_signal_eta[counter],removedSigRT_eta[counter]);
  fillPlotFromVector(num_removed_reco_eta[counter],allRemovedRT_eta[counter]);
  fillPlotFromVector(num_removed_reco_PU_eta[counter],removedPURT_eta[counter]);

  fillPlotFromVector(num_removed_reco_signal_pt[counter],removedSigRT_pt[counter]);
  fillPlotFromVector(num_removed_reco_pt[counter],allRemovedRT_pt[counter]);
  fillPlotFromVector(num_removed_reco_PU_pt[counter],removedPURT_pt[counter]);

  fillPlotFromVector(num_removed_reco_signal_npu[counter],removedSigRT_npu[counter]);
  fillPlotFromVector(num_removed_reco_npu[counter],allRemovedRT_npu[counter]);
  fillPlotFromVector(num_removed_reco_PU_npu[counter],removedPURT_npu[counter]);


  fillPlotFromVector(num_track_reco_signal_eta[counter],allSigRT_eta[counter]);
  fillPlotFromVector(num_track_reco_PU_eta[counter],allPURT_eta[counter]);

  fillPlotFromVector(num_track_reco_signal_pt[counter],allSigRT_pt[counter]);
  fillPlotFromVector(num_track_reco_PU_pt[counter],allPURT_pt[counter]);

  fillPlotFromVector(num_track_reco_signal_npu[counter],allSigRT_npu[counter]);
  fillPlotFromVector(num_track_reco_PU_npu[counter],allPURT_npu[counter]);


  fillPlotFromVector(num_reco_PU_eta[counter],allAssPURT_eta[counter]);

}

void
TrackValidatorAlgos::fillPlotFromVector(TH1F* histo,vector<double> vIn)
{

  for(unsigned ite=0; ite<vIn.size();ite++){
    histo->SetBinContent(ite+1,vIn.at(ite));
  } 

}

void
TrackValidatorAlgos::fillFractionHisto(TH1F* histo,vector<double> num,vector<double> denum,string type)
{

  for (unsigned int j=0; j<num.size(); j++){

    double val = 0.;
    double err = 0.;

    if (denum[j]!=0){
      if (type=="effic"){
	val = ((double) num[j])*1./((double) denum[j]);
        err = sqrt( val*(1-val)/(double) denum[j] );
      } else if (type=="fakerate"){
	val = 1-((double) num[j])*1./((double) denum[j]);
        err = sqrt( val*(1-val)/(double) denum[j] );
      } else if (type=="pileup"){
	val = ((double) num[j])*1./((double) denum[j]);
        err = sqrt( val*(1+val)/(double) denum[j] );
      } else return;
    }

    histo->SetBinContent(j+1, val);
    histo->SetBinError(j+1, err);

  }

}

bool
TrackValidatorAlgos::findRefTrack(const Track& refTrack,const Track& track)
{

	return (
	  refTrack.eta() == track.eta() &&
	  refTrack.phi() == track.phi() &&
	  refTrack.chi2() == track.chi2() &&
	  refTrack.ndof() == track.ndof() &&
	  refTrack.p() == track.p()
	);

}

double
TrackValidatorAlgos::getTrackWeight(const TrackingParticle* tp, const GenJetCollection* genJetColl)
{


  	double trackP = tp->p();

  	for(unsigned jet_ite=0; jet_ite<genJetColl->size(); jet_ite++){

    	  GenJetRef jetRef(genJetColl,jet_ite); 

    	  double jetP = jetRef->p();

 	  vector<const GenParticle*> constituents = jetRef->getGenConstituents();

	  for (unsigned const_ite=0; const_ite<constituents.size(); const_ite++) {

	    if(isGenPart(tp,constituents[const_ite])) return trackP*1./jetP;

  	  }
	
  	}

  	return 0.0001;

}

double
TrackValidatorAlgos::getTrackWeightReco(const Track& track, const PFJetCollection* jetColl)
{

  	double trackP = track.p();

  	for(unsigned jet_ite=0; jet_ite<jetColl->size(); jet_ite++){

    	  PFJetRef jetRef(jetColl,jet_ite); 

    	  double jetP = jetRef->p();

  	  for (unsigned daugth_ite=0; daugth_ite<jetRef->numberOfDaughters(); daugth_ite++){

    	    const PFCandidatePtr pfcand = jetRef->getPFConstituent(daugth_ite);
    	    TrackRef trackref = pfcand->trackRef();
    	    if((trackref.isAvailable()) && (trackref.isNonnull())){
              if(findRefTrack(track,*trackref)) return trackP*1./jetP;
    	    }
  	  }
	
  	}

  	return 0.0001;

}

bool
TrackValidatorAlgos::isGenPart(const TrackingParticle* tp, const GenParticle* gp)
{

	return((fabs(tp->momentum().eta() - gp->eta())<0.001) &&
	       (fabs(tp->momentum().phi() - gp->phi())<0.001) &&
	       (fabs(1.-(tp->pt() / gp->pt()))<0.01));

}