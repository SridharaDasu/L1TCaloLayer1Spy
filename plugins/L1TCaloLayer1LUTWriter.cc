// -*- C++ -*-
//
// Package:    L1Trigger/L1TCaloLayer1
// Class:      L1TCaloLayer1LUTWriter
// 
/**\class L1TCaloLayer1LUTWriter L1TCaloLayer1LUTWriter.cc L1Trigger/L1TCaloLayer1/plugins/L1TCaloLayer1LUTWriter.cc

   Description: This ED Analyzer compares output of CMS L1 Trigger Calo Layer-1 output (CaloTowers) from two sources

   Implementation:
   It is expected that we compare CaloTowers from the spy source to that of the emulator.  
   It can be used to compare any two CaloTower collections
*/
//
// Original Author:  Sridhara Dasu
//         Created:  Sun, 11 Oct 2015 08:14:01 GMT
//
//


// system include files
#include <iostream>
#include <iomanip>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"
#include "CondFormats/DataRecord/interface/L1RCTParametersRcd.h"

#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"

#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"

//
// class declaration
//

class L1TCaloLayer1LUTWriter : public edm::EDAnalyzer {
public:
  explicit L1TCaloLayer1LUTWriter(const edm::ParameterSet&);
  ~L1TCaloLayer1LUTWriter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  bool verbose;

};

L1TCaloLayer1LUTWriter::L1TCaloLayer1LUTWriter(const edm::ParameterSet& iConfig) :
  verbose(iConfig.getUntrackedParameter<bool>("verbose")) {}

L1TCaloLayer1LUTWriter::~L1TCaloLayer1LUTWriter() {}

//
// member functions
//

// ------------ method called for each event  ------------
void
L1TCaloLayer1LUTWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // get rct parameters - these should contain Laura Dodd's tower-level scalefactors (ET, eta)

  edm::ESHandle<L1RCTParameters> rctParameters;
  iSetup.get<L1RCTParametersRcd>().get(rctParameters);
  const L1RCTParameters* rctParameters_ = rctParameters.product();

  // get energy scale to convert input from ECAL - this should be linear with LSB = 0.5 GeV
  edm::ESHandle<L1CaloEcalScale> ecalScale;
  iSetup.get<L1CaloEcalScaleRcd>().get(ecalScale);
  const L1CaloEcalScale* e = ecalScale.product();
      
  // get energy scale to convert input from HCAL - this should be Landsberg's E to ET etc non-linear conversion factors
  edm::ESHandle<L1CaloHcalScale> hcalScale;
  iSetup.get<L1CaloHcalScaleRcd>().get(hcalScale);
  const L1CaloHcalScale* h = hcalScale.product();

  // Loop and write the ECAL LUT
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;
  std::cout << "Input  ECAL_01   ECAL_02   ECAL_03   ECAL_04   ECAL_05   ECAL-06   ECAL_07   ECAL-08   ECAL_09   ECAL_10   ECAL_11   ECAL_12   ECAL_13   ECAL_14   ECAL_15   ECAL_16   ECAL_17   ECAL-18   ECAL_19   ECAL_20   ECAL_21   ECAL_22   ECAL_23   ECAL_24   ECAL_25   ECAL_26   ECAL_27   ECAL_28" << std::endl;
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;
  for(uint32_t fb = 0; fb < 2; fb++) {
    for(uint32_t ecalInput = 0; ecalInput <= 0xFF; ecalInput++) {
      for(int absCaloEta = 1; absCaloEta <= 28; absCaloEta++) {
	double linearizedECalInput = e->et(ecalInput, absCaloEta, 1);
	if(linearizedECalInput != (e->et(ecalInput, absCaloEta, -1))) {
	  std::cerr << "L1TCaloLayer1LUTWriter - ecal scale factors are different for positive and negative eta ! :(" << std::endl;
	}
	// Use hcal = 0 to get ecal only energy but in RCT JetMET scale - should be 8-bit max
	uint32_t value = (rctParameters_->JetMETTPGSum(linearizedECalInput, 0, absCaloEta) / rctParameters_->jetMETLSB());
	if(value > 0xFF) {
	  value = 0xFF;
	}
	if(value == 0) {
	  value = (1 << 11);
	}
	else {
	  uint32_t et_log2 = ((uint32_t) log2(value)) & 0x7;
	  value |= (et_log2 << 12);
	}
	value |= (fb << 10);
	if(absCaloEta == 1) {
	  std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(5) << std::hex
		    << (ecalInput | (fb << 8)) << "   ";
	}
	std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(6) << std::hex
		  << value;
	if(absCaloEta == 28) std::cout << std::endl;
	else std::cout << "    ";
      }
    }
  }

  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;
  std::cout << "Input  HCAL_01   HCAL_02   HCAL_03   HCAL_04   HCAL_05   HCAL-06   HCAL_07   HCAL-08   HCAL_09   HCAL_10   HCAL_11   HCAL_12   HCAL_13   HCAL_14   HCAL_15   HCAL_16   HCAL_17   HCAL-18   HCAL_19   HCAL_20   HCAL_21   HCAL_22   HCAL_23   HCAL_24   HCAL_25   HCAL_26   HCAL_27   HCAL_28" << std::endl;
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;

  // Loop and write the HCAL LUT
  for(uint32_t fb = 0; fb < 2; fb++) {
    for(uint32_t hcalInput = 0; hcalInput <= 0xFF; hcalInput++) {
      for(int absCaloEta = 1; absCaloEta <= 28; absCaloEta++) {
	double linearizedHCalInput = h->et(hcalInput, absCaloEta, 1);
	if(linearizedHCalInput != (h->et(hcalInput, absCaloEta, -1))) {
	  std::cerr << "L1TCaloLayer1LUTWriter - hcal scale factors are different for positive and negative eta ! :(" << std::endl;
	}
	// Use ecal = 0 to get hcal only energy but in RCT JetMET scale - should be 8-bit max
	uint32_t value = (rctParameters_->JetMETTPGSum(0, linearizedHCalInput, absCaloEta) / rctParameters_->jetMETLSB());
	if(value > 0xFF) {
	  value = 0xFF;
	}
	if(value == 0) {
	  value = (1 << 11);
	}
	else {
	  uint32_t et_log2 = ((uint32_t) log2(value)) & 0x7;
	  value |= (et_log2 << 12);
	}
	value |= (fb << 10);
	if(absCaloEta == 1) {
	  std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(5) << std::hex
		    << (hcalInput | (fb << 8)) << "   ";
	}
	std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(6) << std::hex
		  << value;
	if(absCaloEta == 28) std::cout << std::endl;
	else std::cout << "    ";
      }
    }
  }
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1TCaloLayer1LUTWriter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TCaloLayer1LUTWriter::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TCaloLayer1LUTWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TCaloLayer1LUTWriter);
