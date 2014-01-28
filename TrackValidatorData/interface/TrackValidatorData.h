// -*- C++ -*-
//
// Package:    TrackValidatorData
// Class:      TrackValidatorData
// 
/**\class TrackValidatorData TrackValidatorData.cc MGeisler/TrackValidatorData/src/TrackValidatorData.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Fri Feb  3 13:57:40 CET 2012
// $Id: TrackValidatorData.h,v 1.3 2012/03/15 15:44:15 mgeisler Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "MGeisler/TrackValidatorData/interface/TrackValidatorDataAlgos.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"


using namespace std;
using namespace edm;
// using namespace reco;

//
// class declaration
//

class TrackValidatorData : public EDAnalyzer, private TrackValidatorDataAlgos {

  public:
    explicit TrackValidatorData(const ParameterSet&);
    ~TrackValidatorData();

    static void fillDescriptions(ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const Event&, const EventSetup&);

  // ----------member data ---------------------------

  InputTag refMCLabel_;
  vector<InputTag> amLabels_;

  InputTag vcLabel_;

  InputTag wLabel_;

  bool isData_;
  bool isZMuMu_;
  bool ignoremissingtkcollection_;

};
