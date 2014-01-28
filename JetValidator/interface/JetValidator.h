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
#include "MGeisler/JetValidator/interface/JetValidatorAlgos.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

//
// class declaration
//

class JetValidator : public edm::EDAnalyzer, public JetValidatorAlgos {

  public:
    explicit JetValidator(const edm::ParameterSet&);
    ~JetValidator();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&);

    virtual void endRun(edm::Run const&, edm::EventSetup const&);

  // ----------member data ---------------------------

  edm::InputTag genJetLabel_;
  std::vector<edm::InputTag> recoJetLabels_;

  edm::InputTag puLabel_;

};
