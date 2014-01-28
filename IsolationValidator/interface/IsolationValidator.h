// -*- C++ -*-
//
// Package:    IsolationValidator
// Class:      IsolationValidator
// 
/**\class IsolationValidator IsolationValidator.cc MGeisler/IsolationValidator/src/IsolationValidator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Wed Oct 16 10:31:55 CEST 2013
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

// ROOT include files
#include <TH1F.h>
#include <TH2F.h>

//
// class declaration
//

class IsolationValidator : public edm::EDAnalyzer {
  public:
    explicit IsolationValidator(const edm::ParameterSet&);
    ~IsolationValidator();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    std::vector<edm::InputTag> sufLabels_;

    edm::InputTag wLabel_;

    edm::InputTag VertexCollectionLabel_;

    bool isData_;
  
    //parameters for the histograms

    double minIso, maxIso;  int nintIso;
    double minVertcount, maxVertcount;  int nintVertcount;

    //histograms

    std::vector<TH1F*> h_numVtx;

    std::vector<TH2F*> h_phIso_numVtx, h_elIso_numVtx, h_muIso_numVtx;
      
};