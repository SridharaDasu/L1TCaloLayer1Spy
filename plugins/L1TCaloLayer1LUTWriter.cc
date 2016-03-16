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

#include "L1Trigger/L1TCaloLayer1/src/L1TCaloLayer1FetchLUTs.hh"

//
// class declaration
//

class L1TCaloLayer1LUTWriter : public edm::EDAnalyzer {
public:
  explicit L1TCaloLayer1LUTWriter(const edm::ParameterSet&);
  ~L1TCaloLayer1LUTWriter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:

  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------

  std::vector< std::vector< std::vector < uint32_t > > > ecalLUT;
  std::vector< std::vector< std::vector < uint32_t > > > hcalLUT;
  std::vector< std::vector< uint32_t > > hfLUT;

  bool verbose;

};

L1TCaloLayer1LUTWriter::L1TCaloLayer1LUTWriter(const edm::ParameterSet& iConfig) :
  ecalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hcalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hfLUT(12, std::vector< uint32_t >(256)),
  verbose(iConfig.getUntrackedParameter<bool>("verbose")) {}

L1TCaloLayer1LUTWriter::~L1TCaloLayer1LUTWriter() {}

//
// member functions
//

// ------------ method called for each event  ------------
void
L1TCaloLayer1LUTWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(!L1TCaloLayer1FetchLUTs(iSetup, ecalLUT, hcalLUT, hfLUT)) {
    std::cerr << "beginRun: failed to fetch LUTS" << std::endl;
  }

  // Loop and write the ECAL LUT
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;
  std::cout << "Input  ECAL_01   ECAL_02   ECAL_03   ECAL_04   ECAL_05   ECAL-06   ECAL_07   ECAL-08   ECAL_09   ECAL_10   ECAL_11   ECAL_12   ECAL_13   ECAL_14   ECAL_15   ECAL_16   ECAL_17   ECAL-18   ECAL_19   ECAL_20   ECAL_21   ECAL_22   ECAL_23   ECAL_24   ECAL_25   ECAL_26   ECAL_27   ECAL_28" << std::endl;
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;
  for(uint32_t fb = 0; fb < 2; fb++) {
    for(uint32_t ecalInput = 0; ecalInput <= 0xFF; ecalInput++) {
      for(int absCaloEta = 1; absCaloEta <= 28; absCaloEta++) {
	if(absCaloEta == 1) {
	  std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(5) << std::hex
		    << (ecalInput | (fb << 8)) << "   ";
	}
	std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(6) << std::hex
		  << ecalLUT[absCaloEta - 1][fb][ecalInput];
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
	if(absCaloEta == 1) {
	  std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(5) << std::hex
		    << (hcalInput | (fb << 8)) << "   ";
	}
	std::cout << std::showbase << std::internal << std::setfill('0') << std::setw(6) << std::hex
		  << hcalLUT[absCaloEta - 1][fb][hcalInput];
	if(absCaloEta == 28) std::cout << std::endl;
	else std::cout << "    ";
      }
    }
  }
  std::cout << "============================================================================================================================================================================================================================================================================================" << std::endl;

}

// ------------ method called once each job just before starting event loop  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::beginJob()
  {
  }
*/

// ------------ method called once each job just after ending the event loop  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::endJob() 
  {
  }
*/
// ------------ method called when starting to processes a run  ------------
/*
  void 
  L1TCaloLayer1LUTWriter::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
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
