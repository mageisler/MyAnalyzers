#ifndef JetValidatorAlgos_h
#define JetValidatorAlgos_h

// -*- C++ -*-
//
// Package:    JetValidator
// Class:      JetValidator
// 
/**\class JetValidator JetValidator.cc MGeisler/JetValidator/src/JetValidator.cc
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Fri Feb  3 13:57:40 CET 2012
// $Id: TrackValidator.h,v 1.3 2012/03/15 15:44:15 mgeisler Exp $
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

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

// ROOT include files
#include <TH1F.h>
#include <TProfile3D.h>

//
// class declaration
//

class JetValidatorAlgos{
 public:

    JetValidatorAlgos(const edm::ParameterSet&);

    void CreateIntervalVectors();

    void initialize(){setUpVectors();};

    void BookHistos(TFileDirectory);  

    void fill_independent_histos(int, int, int, int);

    pat::JetRef get_best_matching_recoJet(reco::GenJetRef, edm::Handle<pat::JetCollection>);

    reco::GenJetRef get_best_matching_genJet(pat::JetRef, edm::Handle<reco::GenJetCollection>);

    void fill_recoAssociated_genJet_histos(int, reco::GenJetRef, pat::JetRef, int, unsigned);

    void fill_genAssociated_recoJet_histos(int, pat::JetRef, reco::GenJetRef, int);

    void fillHistosFromVectors(int);

 protected:
  //protected functions 

 private: 

  // private methods for internal usage
  void setUpVectors();

   reco::GenJet get_charged_gen_jet(reco::GenJetRef);

   pat::Jet get_charged_reco_jet(pat::JetRef);

  void fillPlotFromVector(TH1F*, std::vector<int>);

  //private data members  

  bool useLogpt_;
  unsigned MaxNumberOfjetsPerEvent_;

  //parameters for the histograms

  double minEta, maxEta;  int nintEta;
  double minPt, maxPt;  int nintPt;
  double minVertcount, maxVertcount;  int nintVertcount;
  double minJetcount, maxJetcount;  int nintJetcount;


  // ###########
  // vectors 
  // ###########

  std::vector<double> etaIntervals;
  std::vector<double> ptIntervals;
  std::vector<int> vertcountIntervals;

  std::vector< std::vector<int> > reco_jets_eta, matchedReco_jets_eta;
  std::vector< std::vector<int> > reco_jets_pt,  matchedReco_jets_pt;
  std::vector< std::vector<int> > reco_jets_npu, matchedReco_jets_npu;

  std::vector< std::vector<int> > gen_jets_eta, matchedGen_jets_eta;
  std::vector< std::vector<int> > gen_jets_pt,  matchedGen_jets_pt;
  std::vector< std::vector<int> > gen_jets_npu, matchedGen_jets_npu;


  // ###########
  // histograms 
  // ###########
 
  std::vector<TH1F*> h_reco_jets_num;
  std::vector<TH1F*> h_reco_jets_eta; std::vector<TH1F*> h_reco_jets_pt; std::vector<TH1F*> h_reco_jets_npu;
  std::vector<TH1F*> h_matchedReco_jets_eta; std::vector<TH1F*> h_matchedReco_jets_pt; std::vector<TH1F*> h_matchedReco_jets_npu;

  std::vector<TH1F*> h_gen_vertex_num; 
  std::vector<TH1F*> h_gen_jets_num;
  std::vector<TH1F*> h_gen_jets_eta; std::vector<TH1F*> h_gen_jets_pt; std::vector<TH1F*> h_gen_jets_npu;
  std::vector<TH1F*> h_matchedGen_jets_eta; std::vector<TH1F*> h_matchedGen_jets_pt; std::vector<TH1F*> h_matchedGen_jets_npu;

  std::vector<TProfile3D*> p_efficiency; std::vector<TProfile3D*> p_fakerate; 

  std::vector<TProfile3D*> p_ptResolution;         std::vector<TProfile3D*> p_ptResponse; 
  std::vector<TProfile3D*> p_ptResolutionCharged;  std::vector<TProfile3D*> p_ptResponseCharged;

  std::vector<TProfile3D*> p_deltaR; 

  std::vector<TProfile3D*> p_constResolution;         std::vector<TProfile3D*> p_constResponse; 
  std::vector<TProfile3D*> p_constResolutionCharged;  std::vector<TProfile3D*> p_constResponseCharged;

};

#endif