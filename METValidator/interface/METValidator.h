// -*- C++ -*-
//
// Package:    METValidator
// Class:      METValidator
// 
/**\class METValidator METValidator.cc MGeisler/METValidator/src/METValidator.cc

*/
//
// Original Author:  Matthias Geisler
//         Created:  Wed Feb  6 15:35:15 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

// ROOT include files
#include <TH1F.h>
#include <TH2F.h>

//
// class declaration
//

class METValidator : public edm::EDAnalyzer {
   public:
      explicit METValidator(const edm::ParameterSet&);
      ~METValidator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      std::vector<edm::InputTag> corrMETLabels_;
      std::vector<edm::InputTag> rawMETLabels_;

      edm::InputTag puLabel_;

      edm::InputTag wLabel_;

      edm::InputTag VertexCollectionLabel_;

  	  bool isData_;

  	  bool usePUInfo_;
  
      //parameters for the histograms

      double minMET, maxMET;  int nintMET;
      double minVertcount, maxVertcount;  int nintVertcount;

      //histograms

      std::vector<TH1F*> h_numVtx;

      std::vector<TH2F*> h_corrMET_Type0;
      std::vector<TH2F*> h_rawMET_Type0;

};