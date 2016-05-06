// -*- C++ -*-
//
// Package:    L1Trigger/L1TCaloLayer1Spy
// Class:      L1TCaloLayer1Spy
// 
/**\class L1TCaloLayer1Spy L1TCaloLayer1Spy.cc L1Trigger/L1TCaloLayer1Spy/plugins/L1TCaloLayer1Spy.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Fri, 09 Oct 2015 10:22:41 GMT
//
//


// system include files
#include <memory>
#include <stdexcept>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

// CTP7 TCP/IP Client 

#include <UCT2016Layer1.hh>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveSample.h"

#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveSample.h"

#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

using namespace l1t;

//
// class declaration
//

class L1TCaloLayer1Spy : public edm::one::EDProducer<> {
   public:
      explicit L1TCaloLayer1Spy(const edm::ParameterSet&);
      ~L1TCaloLayer1Spy();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  std::string setupString;

  uint32_t selectedBXNumber;

  bool verbose;

  uint32_t eventNumber;

  std::unique_ptr<UCT2016Layer1> layer1;

  std::vector< std::vector< std::vector< uint32_t> > > negativeEtaInputData;
  std::vector< std::vector< std::vector< uint32_t> > > positiveEtaInputData;
  std::vector< std::vector< std::vector< uint32_t> > > negativeEtaOutputData;
  std::vector< std::vector< std::vector< uint32_t> > > positiveEtaOutputData;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
L1TCaloLayer1Spy::L1TCaloLayer1Spy(const edm::ParameterSet& iConfig) :
  setupString(iConfig.getUntrackedParameter<std::string>("setupString")),
  selectedBXNumber(iConfig.getUntrackedParameter<uint32_t>("SelectedBXNumber")),
  verbose(iConfig.getUntrackedParameter<bool>("verbose")),
  eventNumber(0),
  layer1(new UCT2016Layer1(setupString)),
  negativeEtaInputData(nLayer1Cards, std::vector< std::vector<uint32_t> >(nInputLinks, std::vector<uint32_t>(nInputWordsPerCapturePerLink))),
  positiveEtaInputData(nLayer1Cards, std::vector< std::vector<uint32_t> >(nInputLinks, std::vector<uint32_t>(nInputWordsPerCapturePerLink))),
  negativeEtaOutputData(nLayer1Cards, std::vector< std::vector<uint32_t> >(nTMTLinks, std::vector<uint32_t>(nOutputWordsPerCapturePerLink))),
  positiveEtaOutputData(nLayer1Cards, std::vector< std::vector<uint32_t> >(nTMTLinks, std::vector<uint32_t>(nOutputWordsPerCapturePerLink)))
{
  produces<HcalTrigPrimDigiCollection>();
  produces<EcalTrigPrimDigiCollection>();
  produces<CaloTowerBxCollection>();
}


L1TCaloLayer1Spy::~L1TCaloLayer1Spy()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1TCaloLayer1Spy::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  static std::string theRunConfiguration;

  using namespace edm;

  std::auto_ptr<EcalTrigPrimDigiCollection> ecalTPGs(new EcalTrigPrimDigiCollection);
  std::auto_ptr<HcalTrigPrimDigiCollection> hcalTPGs(new HcalTrigPrimDigiCollection);
  std::auto_ptr<CaloTowerBxCollection> towersColl (new CaloTowerBxCollection);
  
  // Determine if we need to take action getting data from layer1, and if so do!
  if((eventNumber % nInputEventsPerCapture) == 0) {

    std::cout << "Calling getRunMode() " << std::endl;
    std::cout.flush();

    // Skip processing if layer-1 is in inappropriate mode
    // getRunMode also verifies that the mode is useable for captures
    UCT2016Layer1CTP7::RunMode mode;
    if(!layer1->getRunMode(mode)) {
      std::cerr << "L1TCaloLayer1Spy: Inappropriate mode set for Layer-1; Set mode = " << mode << std::endl;
      return;
    }
    
    // Get configuration for record
    // getConfiguration ensures that the entire layer-1 system is in same configuration
    
    std::cout << "Calling getConfiguration() " << std::endl;
    std::cout.flush();

    std::string configuration;
    if(!layer1->getConfiguration(configuration)) {
      std::cerr << "L1TCaloLayer1Spy: Could not read configuration from CTP7 " << std::endl;
      return;
    }
    if(theRunConfiguration == "") theRunConfiguration = configuration;
    if(configuration != theRunConfiguration) {
      std::cerr << "L1TCaloLayer1Spy: Layer1 configuration changed midstream; Will stop processing from now on"
		<< "; configuration  = " << configuration << " != " << theRunConfiguration << std::endl;
      return;
    }

    std::cout << "getNextCapture() " << std::endl;
    std::cout.flush();

    UCT2016Layer1CTP7::CaptureMode captureMode;
    UCT2016Layer1CTP7::CaptureStatus captureStatus;
    uint32_t nCaptured = layer1->getNextCapture(captureMode, selectedBXNumber, captureStatus, 
						negativeEtaInputData, positiveEtaInputData,
						negativeEtaOutputData, positiveEtaOutputData);
    if(nCaptured != nInputEventsPerCapture) {
      std::cerr << "L1TCaloLayer1Spy: Layer1 could not make a capture" << std::endl;
      return;
    }

  }

  // Trigger system processes one BX per "event", and the
  // users of this products should expect the hit BX data only.
  // Use data for the event from the layer1 buffers incrementing 
  // to pick the next BX

  int theBX = 0;

  // Produce input data
  // Loop over all cards
  for(uint32_t cardPhi = 0; cardPhi < 18; cardPhi++) {
    // Loop over both sides
    for(int cardSide = -1; cardSide <= 1; cardSide+=2) {
      uint32_t offset = (eventNumber % nInputEventsPerCapture) * 4;
      // Loop over all input links
      for(uint32_t link = 0; link < 14; link++) {
	// Each event eats four words of input buffers
	uint32_t *ecalLinkData = 0;
	uint32_t *hcalLinkData = 0;
	if(cardSide == 1) {
	  ecalLinkData = &(positiveEtaInputData[cardPhi][2*link].data())[offset];
	  hcalLinkData = &(positiveEtaInputData[cardPhi][2*link+1].data())[offset];
	}
	else {
	  ecalLinkData = &(negativeEtaInputData[cardPhi][2*link].data())[offset];
	  hcalLinkData = &(negativeEtaInputData[cardPhi][2*link+1].data())[offset];
	}
	// Bottom eight bits of the third word contain ECAL finegrain feature bits
	// Store them for later access in the loop
	uint8_t ecalFBits = (ecalLinkData[2] & 0xFF);
	// The third 32-bit word + bottom 16 bits of the fourth word make up 
	// 6-bit feature word for each of the eight towers.  They are stitched
	// together in a 64-bit word here to be pealed of as needed later in loops
	uint64_t hcalFBits = hcalLinkData[2];
	hcalFBits |= (((uint64_t) (hcalLinkData[3] & 0xFFFF)) << 32);
	// Process all Eta in a link
	for(uint32_t dEta = 0; dEta < 2; dEta++) {
	  uint32_t ecalDataWord = ecalLinkData[dEta];
	  uint32_t hcalDataWord = hcalLinkData[dEta];
	  // Process all Phi in a link
	  for(uint32_t dPhi = 0; dPhi < 4; dPhi++) {
	    // Determine tower data and location in (caloEta, caloPhi)
	    int absCaloEta = link*2 + dEta + 1;
	    int caloEta = cardSide * absCaloEta;
	    int caloPhi = cardPhi * 4 + dPhi - 1;
	    if(caloPhi <= 0)caloPhi += 72;
	    // Make ECALTriggerPrimitive
	    uint32_t em = (ecalDataWord >> (dPhi * 8)) & (0xFF);
	    bool efb = ((ecalFBits & (0x1 << (dEta * 4 + dPhi))) != 0);
	    uint16_t towerDatum = em;
	    if(efb) towerDatum |= 0x0100;
	    EcalTriggerPrimitiveSample sample(towerDatum); 
	    EcalSubdetector ecalTriggerTower = EcalSubdetector::EcalTriggerTower;
	    EcalTrigTowerDetId id(cardSide, ecalTriggerTower, absCaloEta, caloPhi);
	    EcalTriggerPrimitiveDigi etpg(id);
	    etpg.setSize(1);
	    etpg.setSample(0, sample);
	    ecalTPGs->push_back(etpg);
	    // Make HCALTriggerPrimitive
	    uint32_t hd = (hcalDataWord >> (dPhi * 8)) & (0xFF);
	    uint8_t hfb = (hcalFBits >> (dEta * 4 + dPhi)*6) & 0x3F;
	    towerDatum = (hd + (hfb << 8));
	    HcalTriggerPrimitiveSample hSample(towerDatum); 
	    HcalTrigTowerDetId hid(caloEta, caloPhi);
	    HcalTriggerPrimitiveDigi htpg(hid);
	    htpg.setSize(1);
	    htpg.setSample(0, hSample);
	    hcalTPGs->push_back(htpg);
	  }
	}
      }
      // HF Data
      for(uint32_t abPhi=0; abPhi<2; ++abPhi) {
        uint32_t *hfData;
        if(cardSide == 1) {
          hfData = &(positiveEtaInputData[cardPhi][28+abPhi].data())[offset];
        } else {
          hfData = &(negativeEtaInputData[cardPhi][28+abPhi].data())[offset];
        }
        for(uint32_t hfEta=0; hfEta<12; ++hfEta) {
          if ( abPhi == 0 && hfEta == 11 ) continue;
          if ( abPhi == 1 && hfEta == 10 ) continue;
          size_t word = hfEta/4;
          size_t shift = 8 * (hfEta%4);
          if ( hfEta > 10 ) shift = 16;
          uint32_t et = (hfData[word]>>shift) & 0xff;

          int iEta = cardSide * (hfEta+30);
          int iPhi = 1 + cardPhi*4 + abPhi*2;
          if(std::abs(iEta) == 41) iPhi -= 2; // Last two HF are 3, 7, 11, ...
          iPhi = (iPhi+69)%72 + 1; // iPhi -= 2 mod 72
          std::cout << iEta << ", " << iPhi << ", " << et << std::endl;
          HcalTriggerPrimitiveSample sample(et); 
          HcalTrigTowerDetId id(iEta, iPhi);
          id.setVersion(1); // To not process these 1x1 HF TPGs with RCT
          HcalTriggerPrimitiveDigi tpg(id);
          tpg.setSize(1);
          tpg.setSample(0, sample);
          hcalTPGs->push_back(tpg);
        }
      }

      // Output data
      // Make caloTower collection just for Barrel and Endcap for the moment
      uint32_t nHeader = 1;
      uint32_t nBEDataWords = 28;
      for(uint32_t tEta = 0; tEta < nBEDataWords; tEta++) {
	for(uint32_t dPhi = 0; dPhi < 4; dPhi++) {
	  uint32_t outputLink = (eventNumber % nTMTCards)*2 + (dPhi/2);
	  //	  uint32_t offset = ((eventNumber % nOutputEventsPerCapturePerLinkPair) / nTMTCards) * nOutputEventWords + nHeader + tEta;
	  uint32_t offset = (eventNumber/nTMTCards) * nOutputEventWords + nHeader + tEta;
	  uint16_t dataWord;
	  if(cardSide == -1) {
	    dataWord = (((negativeEtaOutputData[cardPhi][outputLink].data())[offset]) >> ((dPhi%2) * 16)) & 0xFFFF;
	  }
	  else {
	    dataWord = (((positiveEtaOutputData[cardPhi][outputLink].data())[offset]) >> ((dPhi%2) * 16)) & 0xFFFF;
	  }
	  CaloTower caloTower;
	  caloTower.setHwPt(dataWord & 0x1FF);               // Bits 0-8 of the 16-bit word per the interface protocol document
	  caloTower.setHwEtRatio((dataWord >> 9) & 0x7);     // Bits 9-11 of the 16-bit word per the interface protocol document
	  caloTower.setHwQual((dataWord >> 12) & 0xF);       // Bits 12-15 of the 16-bit word per the interface protocol document
	  // Determine tower data and location in (caloEta, caloPhi)
	  int absCaloEta = tEta + 1;
	  int caloEta = cardSide * absCaloEta;
	  int caloPhi = cardPhi * 4 + dPhi - 1;
	  if(caloPhi <= 0)caloPhi += 72;
	  caloTower.setHwEta(caloEta);
	  caloTower.setHwPhi(caloPhi);
	  // Push the tower in
	  towersColl->push_back(theBX, caloTower);
	}
      }
    }
  }

  iEvent.put(ecalTPGs);
  iEvent.put(hcalTPGs);
  iEvent.put(towersColl);

  eventNumber++;

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1TCaloLayer1Spy::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TCaloLayer1Spy::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TCaloLayer1Spy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
L1TCaloLayer1Spy::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TCaloLayer1Spy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TCaloLayer1Spy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TCaloLayer1Spy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TCaloLayer1Spy);
