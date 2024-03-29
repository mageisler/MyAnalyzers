#ifndef TrackValidatorAlgos_h
#define TrackValidatorAlgos_h

// -*- C++ -*-
//
// Package:    TrackValidator
// Class:      TrackValidator
// 
/**\class TrackValidator TrackValidator.cc MGeisler/TrackValidator/src/TrackValidator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Fri Feb  3 13:57:40 CET 2012
// $Id: TrackValidatorAlgos.h,v 1.8 2012/07/10 08:34:08 mgeisler Exp $
//
//

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

#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// ROOT include files
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile3D.h>

using namespace std;
using namespace edm;
using namespace reco;

typedef vector<pair<TrackingParticleRef, double> > TpDoubV;

//
// class declaration
//

class TrackValidatorAlgos{
 public:

    TrackValidatorAlgos(const edm::ParameterSet&);

    void CreateIntervalVectors();

    void GetInputCollections(const Event&);

    void initialize(){setUpVectors();};

    void BookHistos(TFileDirectory); 

    void fill_independent_histos(int,int,int,unsigned);

    void fill_recoAssociated_simTrack_histos(int,TrackingParticle*,const Track*,int,unsigned*);

    void fill_simAssociated_recoTrack_histos(int,const Track&,bool,bool,int,TpDoubV);

    void fill_removedRecoTrack_histos(int,const Track&,bool,bool,int,TpDoubV);

    void fillFractionHistosFromVectors(int);

    void fillHistosFromVectors(int);

    bool findRefTrack(const Track&,const Track&);

 protected:
  //protected functions 

 private: 

  // private methods for internal usage
  void setUpVectors();

  void fillPlotFromVector(TH1F*,vector<double>);

  void fillFractionHisto(TH1F*,vector<double>,vector<double>,string);

  double getTrackWeight(const TrackingParticle*, const GenJetCollection*);

  double getTrackWeightReco(const Track&, const PFJetCollection*);

  bool isGenPart(const TrackingParticle*, const GenParticle*);


  //private data members  

  bool useLogpt_;

  bool useJetWeighting_;
  string genJetCollLabel_;
  string jetCollLabel_;

  Handle<GenJetCollection> genJetCollH;
  Handle<PFJetCollection> jetCollH;

  TrackingParticleSelector* generalTpSignalSelector;
  TrackingParticleSelector* generalTpSignalSelectorPt;
  TrackingParticleSelector* generalTpPUSelector;

  //parameters for the histograms

  double minEta, maxEta;  int nintEta;
  double minEtaErr, maxEtaErr;  int nintEtaErr;
  double minPt, maxPt;  int nintPt;
  double minPtErr, maxPtErr;  int nintPtErr;
  double minPtErrRel, maxPtErrRel;  int nintPtErrRel;
  double minTrackcount, maxTrackcount;  int nintTrackcount;
  double minContribution, maxContribution;  int nintContribution;
  double minTW, maxTW;  int nintTW;
  int minVertcount, maxVertcount, nintVertcount;


  // ###########
  // profiles 
  // ###########

  vector<TProfile3D*> p_efficiency_val; vector<TProfile3D*> p_fakerate_val; 
  vector<TProfile3D*> p_efficiencyPileUp_val; vector<TProfile3D*> p_fakeratePileUp_val;

  vector<TProfile3D*> p_fakerate_err; 
  vector<TProfile3D*> p_efficiencyPileUp_err; vector<TProfile3D*> p_fakeratePileUp_err; 

  vector<TProfile3D*> p_fakerate_relerr; 
  vector<TProfile3D*> p_efficiencyPileUp_relerr; vector<TProfile3D*> p_fakeratePileUp_relerr; 

  vector<TProfile*> p_purity_tw; 


  // ###########
  // histograms 
  // ###########

  // track collection

  vector<TH2F*> effic_npu_Contr; 
  vector<TH2F*> num_simul_tracks_npu_Contr; vector<TH2F*> num_assoc_npu_Contr;

  vector<TH2F*> fakerate_npu_Contr; vector<TH2F*> fakerate_npu_Contr_help;
  vector<TH2F*> num_reco_tracks_npu_Contr; vector<TH2F*> num_assoc2_npu_Contr;

  vector<TH1F*> PU_effic_eta; vector<TH1F*> PU_effic_pt; vector<TH1F*> PU_effic_npu;

  vector<TH1F*> PU_fakerate_1_eta; vector<TH1F*> PU_fakerate_1_pt; 
  vector<TH1F*> PU_fakerate_1_npu;
  vector<TH1F*> PU_fakerate_2_eta; vector<TH1F*> PU_fakerate_2_pt;  
  vector<TH1F*> PU_fakerate_2_npu;

  vector<TH1F*> effic_eta; vector<TH1F*> effic_pt; 
  vector<TH1F*> effic_npu;
  vector<TH1F*> fakerate_eta; vector<TH1F*> fakerate_pt; 
  vector<TH1F*> fakerate_npu;

  vector<TH1F*> num_simul_tracks; vector<TH1F*> num_track_simul_eta;
  vector<TH1F*> num_track_simul_pt; vector<TH1F*> num_track_simul_npu; 
  vector<TH1F*> num_simul_vertex; 
  vector<TH1F*> num_reco_tracks; vector<TH1F*> num_track_reco_eta;
  vector<TH1F*> num_track_reco_pt; vector<TH1F*> num_track_reco_npu;

  vector<TH1F*> num_removed_reco_signal_eta; vector<TH1F*> num_removed_reco_eta;
  vector<TH1F*> num_removed_reco_PU_eta; vector<TH1F*> num_reco_PU_eta;

  vector<TH1F*> num_removed_reco_signal_pt; vector<TH1F*> num_removed_reco_pt;
  vector<TH1F*> num_removed_reco_PU_pt;

  vector<TH1F*> num_removed_reco_signal_npu; vector<TH1F*> num_removed_reco_npu;
  vector<TH1F*> num_removed_reco_PU_npu;

  vector<TH1F*> num_assoc_eta; vector<TH1F*> num_assoc2_eta;
  vector<TH1F*> num_assoc_pt; vector<TH1F*> num_assoc2_pt;
  vector<TH1F*> num_assoc_npu; vector<TH1F*> num_assoc2_npu;

  vector<TH1F*> num_track_reco_PU_eta; vector<TH1F*> num_track_reco_signal_eta;
  vector<TH1F*> num_track_reco_PU_pt; vector<TH1F*> num_track_reco_signal_pt;
  vector<TH1F*> num_track_reco_PU_npu; vector<TH1F*> num_track_reco_signal_npu;

  vector<TH1F*> weights;
  vector<TH1F*> h_tweight;


  // ###########
  // vectors 
  // ###########

  vector<double> etaintervals;
  vector<double> ptintervals;
  vector<double> vertcountintervals;

  // track collection

  vector< vector<double> > allSignalTP_eta, allRT_eta;
  vector< vector<double> > allSignalTP_npu, allRT_npu;
  vector< vector<double> > allSignalTP_pt,  allRT_pt;

  vector< vector<double> > assSignalTP_eta, assSignalTP_npu, assSignalTP_pt;
  vector< vector<double> > assSignalRT_eta, assSignalRT_npu, assSignalRT_pt;

  vector< vector<double> > allSigRT_eta, allPURT_eta, allAssPURT_eta;
  vector< vector<double> > allSigRT_pt,  allPURT_pt;
  vector< vector<double> > allSigRT_npu, allPURT_npu;

  vector< vector<double> > allRemovedRT_eta, allRemovedRT_npu, allRemovedRT_pt;
  vector< vector<double> > removedSigRT_eta, removedSigRT_npu, removedSigRT_pt;
  vector< vector<double> > removedPURT_eta,  removedPURT_npu,  removedPURT_pt;

  vector<int> sim_tracks;

};

#endif