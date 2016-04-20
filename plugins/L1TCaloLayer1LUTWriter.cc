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

#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include <openssl/md5.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"

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

  bool writeXMLParam(std::string id, std::string type, std::string body);
  bool writeSWATCHVector(std::string id, std::vector<int> vect);
  bool writeSWATCHVector(std::string id, std::vector<double> vect);
  bool writeSWATCHTableRow(std::vector<uint32_t> vect);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  // Wrapper for xmllib error codes
  // returnCode < 0 if error
  inline bool rcWrap(int rc) const {
    if ( rc < 0 ) {
      edm::LogError("L1TCaloLayer1LUTWriter") << "Error while processing an xmllib command :<";
      return false;
    }
    return true;
  };

  // ----------member data ---------------------------

  bool useLSB;
  bool useCalib;
  bool useECALLUT;
  bool useHCALLUT;
  bool useHFLUT;

  std::vector< std::vector< std::vector < uint32_t > > > ecalLUT;
  std::vector< std::vector< std::vector < uint32_t > > > hcalLUT;
  std::vector< std::vector< uint32_t > > hfLUT;

  bool verbose;
  xmlTextWriterPtr writer_;
};

L1TCaloLayer1LUTWriter::L1TCaloLayer1LUTWriter(const edm::ParameterSet& iConfig) :
  useLSB(iConfig.getParameter<bool>("useLSB")),
  useCalib(iConfig.getParameter<bool>("useCalib")),
  useECALLUT(iConfig.getParameter<bool>("useECALLUT")),
  useHCALLUT(iConfig.getParameter<bool>("useHCALLUT")),
  useHFLUT(iConfig.getParameter<bool>("useHFLUT")),
  ecalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hcalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hfLUT(12, std::vector< uint32_t >(256)),
  verbose(iConfig.getUntrackedParameter<bool>("verbose"))
{
  std::string fileName = iConfig.getParameter<std::string>("fileName");

  writer_ = xmlNewTextWriterFilename(fileName.c_str(), 0);
  if (writer_ == NULL) {
    edm::LogError("L1TCaloLayer1LUTWriter") << ("testXmlwriterFilename: Error creating the xml writer");
  }
  else {
    // See these links for why this is bad:
    // https://www.w3.org/TR/2008/REC-xml-20081126/#sec-white-space
    // http://usingxml.com/Basics/XmlSpace
    xmlTextWriterSetIndent(writer_, 1);
    // If curious about all this BAD_CAST stuff:
    // http://xmlsoft.org/html/libxml-xmlstring.html#BAD_CAST
    xmlTextWriterSetIndentString(writer_, BAD_CAST "  ");
  }
}


L1TCaloLayer1LUTWriter::~L1TCaloLayer1LUTWriter()
{
  xmlFreeTextWriter(writer_);
}


//
// member functions
//

bool
L1TCaloLayer1LUTWriter::writeXMLParam(std::string id, std::string type, std::string body)
{
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "param")) ) return false;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "id", BAD_CAST id.c_str())) ) return false;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "type", BAD_CAST type.c_str())) ) return false;
  if ( !rcWrap(xmlTextWriterWriteString(writer_, BAD_CAST body.c_str())) ) return false;
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return false;

  // Success!
  return true;
}


bool
L1TCaloLayer1LUTWriter::writeSWATCHTableRow(std::vector<uint32_t> vect)
{
  std::stringstream output;
  for(auto it=vect.begin(); it!=vect.end(); ++it) {
    output << std::showbase << std::internal << std::setfill('0') << std::setw(6) << std::hex << *it;
    if ( it != vect.end()-1 ) {
      output << ", ";
    }
  }

  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "row", BAD_CAST output.str().c_str())) ) return false;
  return true;
}

bool
L1TCaloLayer1LUTWriter::writeSWATCHVector(std::string id, std::vector<int> vect)
{
  std::stringstream output;
  for(auto it=vect.begin(); it!=vect.end(); ++it) {
    output << *it;
    if ( it != vect.end()-1 ) {
      output << ", ";
    }
  }
  return writeXMLParam(id, "vector:int", output.str());
}

bool
L1TCaloLayer1LUTWriter::writeSWATCHVector(std::string id, std::vector<double> vect)
{
  std::stringstream output;
  for(auto it=vect.begin(); it!=vect.end(); ++it) {
    output << (float) *it;
    if ( it != vect.end()-1 ) {
      output << ", ";
    }
  }
  return writeXMLParam(id, "vector:float", output.str());
}

// ------------ method called for each event  ------------
void
L1TCaloLayer1LUTWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Can't bail out of constructor, so bail here
  if ( writer_ == NULL ) return;

  // CaloParams contains all persisted parameters for Layer 1
  edm::ESHandle<l1t::CaloParams> paramsHandle;
  iSetup.get<L1TCaloParamsRcd>().get(paramsHandle);
  if ( paramsHandle.product() == nullptr ) {
    edm::LogError("L1TCaloLayer1LUTWriter") << "Missing CaloParams object! Check Global Tag, etc.";
    return;
  }
  l1t::CaloParamsHelper caloParams(*paramsHandle.product());

  // Helper function translates CaloParams into actual LUT vectors
  if(!L1TCaloLayer1FetchLUTs(iSetup, ecalLUT, hcalLUT, hfLUT, useLSB, useCalib, useECALLUT, useHCALLUT, useHFLUT)) {
    edm::LogError("L1TCaloLayer1LUTWriter") << "Failed to fetch LUTs";
    return;
  }

  if ( !rcWrap(xmlTextWriterStartDocument(writer_, NULL, NULL, NULL)) ) return;

  // Root node <algo>
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "algo")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "id", BAD_CAST "calol1")) ) return;


  // SWATCH magic for all cards
  // If we ever have different LUTs for different phi, this would need to change
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "context")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "id", BAD_CAST "processors")) ) return;

  // LUT generation parameters
  // This is not needed for SWATCH
  // but necessary for O2O, given the offline format
  // (i.e. what we are reading right now)
  if ( !writeSWATCHVector("layer1ECalScaleETBins", caloParams.layer1ECalScaleETBins()) ) return;
  if ( !writeSWATCHVector("layer1ECalScaleFactors", caloParams.layer1ECalScaleFactors()) ) return;
  if ( !writeSWATCHVector("layer1HCalScaleETBins", caloParams.layer1HCalScaleETBins()) ) return;
  if ( !writeSWATCHVector("layer1HCalScaleFactors", caloParams.layer1HCalScaleFactors()) ) return;
  if ( !writeSWATCHVector("layer1HFScaleETBins", caloParams.layer1HFScaleETBins()) ) return;
  if ( !writeSWATCHVector("layer1HFScaleFactors", caloParams.layer1HFScaleFactors()) ) return;
  if ( !writeXMLParam("towerLsbSum", "float", std::to_string(caloParams.towerLsbSum())) ) return;
  if ( !writeXMLParam("useLSB", "bool", (useLSB) ? "true":"false") ) return;
  if ( !writeXMLParam("useCalib", "bool", (useCalib) ? "true":"false") ) return;
  if ( !writeXMLParam("useECALLUT", "bool", (useECALLUT) ? "true":"false") ) return;
  if ( !writeXMLParam("useHCALLUT", "bool", (useHCALLUT) ? "true":"false") ) return;
  if ( !writeXMLParam("useHFLUT", "bool", (useHFLUT) ? "true":"false") ) return;

  // We will checksum the LUT contents and put it at the end
  MD5_CTX md5context;
  MD5_Init(&md5context);

  // Loop and write the ECAL LUT
  // <param id="ECALLUT" type="table">
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "param")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "id", BAD_CAST "ECALLUT")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "type", BAD_CAST "table")) ) return;

  // <columns>
  const char * ecal_columns{"Input, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28"};
  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "columns", BAD_CAST ecal_columns)) ) return;

  // <types>
  const char * ecal_types{"uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint"};
  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "types", BAD_CAST ecal_types)) ) return;

  // <rows>
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "rows")) ) return;

  for(uint32_t fb = 0; fb < 2; fb++) {
    for(uint32_t ecalInput = 0; ecalInput <= 0xFF; ecalInput++) {
      std::vector<uint32_t> row;
      uint32_t fullInput = (fb << 8) | ecalInput;
      row.push_back(fullInput);
      for(int iEta=1; iEta<=28; ++iEta) {
        row.push_back(ecalLUT[iEta-1][fb][ecalInput]);
      }
      MD5_Update(&md5context, &row[1], (row.size()-1)*sizeof(uint32_t));
      if ( !writeSWATCHTableRow(row) ) return;
    }
  }

  // </rows>
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return;
  // </param>
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return;


  // Loop and write the HCAL LUT
  // <param id="HCALLUT" type="table">
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "param")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "id", BAD_CAST "HCALLUT")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "type", BAD_CAST "table")) ) return;

  // <columns>
  const char * hcal_columns{"Input, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28"};
  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "columns", BAD_CAST hcal_columns)) ) return;

  // <types>
  const char * hcal_types{"uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint"};
  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "types", BAD_CAST hcal_types)) ) return;

  // <rows>
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "rows")) ) return;

  for(uint32_t fb = 0; fb < 2; fb++) {
    for(uint32_t hcalInput = 0; hcalInput <= 0xFF; hcalInput++) {
      std::vector<uint32_t> row;
      uint32_t fullInput = (fb << 8) | hcalInput;
      row.push_back(fullInput);
      for(int iEta=1; iEta<=28; ++iEta) {
        row.push_back(hcalLUT[iEta-1][fb][hcalInput]);
      }
      MD5_Update(&md5context, &row[1], (row.size()-1)*sizeof(uint32_t));
      if ( !writeSWATCHTableRow(row) ) return;
    }
  }

  // </rows>
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return;
  // </param>
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return;


  // Loop and write the HF LUT
  // <param id="HFLUT" type="table">
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "param")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "id", BAD_CAST "HFLUT")) ) return;
  if ( !rcWrap(xmlTextWriterWriteAttribute(writer_, BAD_CAST "type", BAD_CAST "table")) ) return;

  // <columns>
  const char * hf_columns{"Input, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41"};
  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "columns", BAD_CAST hf_columns)) ) return;

  // <types>
  const char * hf_types{"uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint, uint"};
  if ( !rcWrap(xmlTextWriterWriteElement(writer_, BAD_CAST "types", BAD_CAST hf_types)) ) return;

  // <rows>
  if ( !rcWrap(xmlTextWriterStartElement(writer_, BAD_CAST "rows")) ) return;

  for(uint32_t fb = 0; fb < 4; fb++) {
    for(uint32_t hfInput = 0; hfInput <= 0xFF; hfInput++) {
      std::vector<uint32_t> row;
      uint32_t fullInput = (fb << 8) | hfInput;
      row.push_back(fullInput);
      for(int hfEta=0; hfEta<12; ++hfEta) {
        uint32_t output = hfLUT[hfEta][hfInput];
        // HF LUT in emulator does not currently handle
        // feature bits, instead emulator passes them
        // unaltered. So this is what we have hardware do
        output |= (fb << 8);
        row.push_back(output);
      }
      MD5_Update(&md5context, &row[1], (row.size()-1)*sizeof(uint32_t));
      if ( !writeSWATCHTableRow(row) ) return;
    }
  }

  // </rows>
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return;
  // </param>
  if ( !rcWrap(xmlTextWriterEndElement(writer_)) ) return;


  // Now to write the checksum
  unsigned char checksum[MD5_DIGEST_LENGTH];
  MD5_Final(checksum, &md5context);
  std::stringstream checksumString;
  for(size_t i=0; i<MD5_DIGEST_LENGTH; ++i) {
    checksumString << std::hex << static_cast<unsigned int>(checksum[i]);
  }
  if ( !writeXMLParam("md5checksum", "string", checksumString.str()) ) return;


  // Closes all open elements recursively for us
  if ( !rcWrap(xmlTextWriterEndDocument(writer_)) ) return;
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
