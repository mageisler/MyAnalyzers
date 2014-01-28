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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

// ROOT include files
#include <TProfile3D.h>


//
// class declaration
//

class ZDaughtersFinder : public edm::EDAnalyzer {
   public:
      explicit ZDaughtersFinder(const edm::ParameterSet&);
      ~ZDaughtersFinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      
      bool isInVector(int, std::vector<int>);
      
      void addDaughters(reco::GenParticle, reco::GenParticleCollection*, std::vector<int>);
      
      void isGenJet(HepMC::GenParticle, reco::GenJetCollection*, std::vector<reco::GenJetRef>*);
      
      bool isGenPart(HepMC::GenParticle, reco::GenParticle);

      pat::JetRef get_best_matching_recoJet(reco::GenJetRef, edm::Handle<pat::JetCollection>);
      
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------


      std::vector<edm::InputTag> input_recoTrackColls_;
      std::vector<edm::InputTag> input_recoJetColls_;
      
      std::vector<int> input_pdgIds_;
      
      edm::InputTag input_genParticles_;
      edm::InputTag input_genJets_;
      edm::InputTag input_trackingParticles_;
      
      edm::InputTag input_puSummary_;
      
      double input_maxDist_;
      double input_maxDeltaR_;

      TrackingParticleSelector* generalTpSignalSelector;

      //parameters for the histograms
 
      double minEta, maxEta;  int nintEta;
      double minPt, maxPt;  int nintPt;
      int minVertcount, maxVertcount, nintVertcount;


      // ###########
      // profiles 
      // ###########

      std::vector<TProfile3D*> p_trackEfficiency;
      std::vector<TProfile3D*> p_jetEfficiency; std::vector<TProfile3D*> p_jetPtResponse;
      
};