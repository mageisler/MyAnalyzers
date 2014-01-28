// -*- C++ -*-
//
// Package:    MuonSelection
// Class:      MuonSelection
// 
/**\class MuonSelection MuonSelection.cc MGeisler/MuonSelection/src/MuonSelection.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Sep  3 11:58:16 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"


//
// class declaration
//

class MuonSelection : public edm::EDProducer {
   public:
      explicit MuonSelection(const edm::ParameterSet&);
      ~MuonSelection();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      
      edm::InputTag muonLabel_;
      edm::InputTag vertexCollectionLabel_;
      
      int minLayers_;
      double minPt_;
      double maxEta_;
      
};