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

#include <UCT2016Layer1CTP7.hh>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

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

class L1TCaloLayer1Spy : public edm::EDProducer {
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

  std::string phiMapFile;

  uint32_t selectedOrbitNumber;
  uint32_t selectedBXNumber;

  bool verbose;

  std::vector<UCT2016Layer1CTP7*> cards;

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
  phiMapFile(iConfig.getUntrackedParameter<std::string>("Layer1PhiMapXMLFile")),
  selectedOrbitNumber(iConfig.getUntrackedParameter<uint32_t>("SelectedOrbitNumber")),
  selectedBXNumber(iConfig.getUntrackedParameter<uint32_t>("SelectedBXNumber")),
  verbose(iConfig.getUntrackedParameter<bool>("verbose")) 
{
  produces<HcalTrigPrimDigiCollection>();
  produces<EcalTrigPrimDigiCollection>();
  produces<CaloTowerBxCollection>();
  for(int phi = 0; phi <= 17; phi++) {
    std::cout << "Connecting to phi=" << phi << std::endl;
    try {
      cards.push_back(new UCT2016Layer1CTP7(phi, phiMapFile));
    }
    catch (std::runtime_error &e) {
      std::cout << "Failed connecting to phi=" << phi << "; " << e.what() << std::endl;
    }
  }
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

  using namespace edm;

  std::auto_ptr<EcalTrigPrimDigiCollection> ecalTPGs(new EcalTrigPrimDigiCollection);
  std::auto_ptr<HcalTrigPrimDigiCollection> hcalTPGs(new HcalTrigPrimDigiCollection);
  std::auto_ptr<CaloTowerBxCollection> towersColl (new CaloTowerBxCollection);

  // Trigger system processes one BX per "event", and the
  // users of this products should expect the hit BX data only.

  int theBX = 0;

  // FIXME: For testing purpose use some random data :)

  for(int caloPhi = 1; caloPhi <= 72; caloPhi++) {
    for(int caloEta = -28; caloEta <= 28; caloEta++) {
      if(caloEta != 0) {
	// Make random ECALTrigPrimitives
	uint32_t em = (random() & 0xFF);
	bool efb = (random() & 0x1);
	uint16_t towerDatum = em;
	if(efb) towerDatum |= 0x0100;
	EcalTriggerPrimitiveSample sample(towerDatum); 
	int iEta = abs(caloEta);
	int zSide = caloEta / iEta;
	EcalSubdetector ecalTriggerTower = EcalSubdetector::EcalTriggerTower;
	EcalTrigTowerDetId id(zSide, ecalTriggerTower, iEta, caloPhi);
	EcalTriggerPrimitiveDigi etpg(id);
	etpg.setSize(1);
	etpg.setSample(0, sample);
	ecalTPGs->push_back(etpg);
	// Make random HCALTrigPrimitives
	uint32_t hd = (random() & 0xFF);
	bool hfb = (random() & 0x1);
	towerDatum = em;
	if(hfb) towerDatum |= 0x0100;
	HcalTriggerPrimitiveSample hSample(towerDatum); 
	HcalTrigTowerDetId hid(caloEta, caloPhi);
	HcalTriggerPrimitiveDigi htpg(hid);
	htpg.setSize(1);
	htpg.setSample(0, hSample);
	hcalTPGs->push_back(htpg);
	// Make caloTowers
	CaloTower caloTower;
	uint32_t et = em + hd;
	uint32_t er = 0;
	uint32_t qualBits = 0;
	if(em == 0 || hd == 0) qualBits |= 0x01;
	if(hd == 0 || em >= hd) qualBits |= 0x10;
	if(em == hd) {
	  er = 0;
	}
	else if(em != 0 && em > hd) {
	  er = (uint32_t) log2(((double) hd) / ((double) em));
	}
	else if(hd != 0 && hd > em) {
	  er = (uint32_t) log2(((double) em) / ((double) hd));
	}
	if(er > 0x7) er = 7;
	caloTower.setHwPt(et);               // Bits 0-8 of the 16-bit word per the interface protocol document
	caloTower.setHwEtRatio(er);          // Bits 9-11 of the 16-bit word per the interface protocol document
	caloTower.setHwQual(qualBits);       // Bits 12-15 of the 16-bit word per the interface protocol document
	caloTower.setHwEta(caloEta);
	caloTower.setHwPhi(caloPhi);
	caloTower.setHwEtEm(em);             // This is provided as a courtesy - not available to hardware
	caloTower.setHwEtHad(hd);            // This is provided as a courtesy - not available to hardware
	towersColl->push_back(theBX, caloTower);
      }
    }
  }

  iEvent.put(ecalTPGs);
  iEvent.put(hcalTPGs);
  iEvent.put(towersColl);

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
