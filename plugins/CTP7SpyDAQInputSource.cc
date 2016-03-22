// N. Smith <nick.smith@cern.ch>
// For reading CTP7 ethernet captures directly into CMSSW framework
//

#include <sys/time.h>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "DataFormats/Provenance/interface/BranchID.h"
#include "DataFormats/Provenance/interface/BranchType.h"
#include "DataFormats/Provenance/interface/LuminosityBlockAuxiliary.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "DataFormats/Provenance/interface/ProcessHistoryID.h"
#include "DataFormats/Provenance/interface/ProcessHistoryRegistry.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"
#include "DataFormats/Provenance/interface/ProductRegistry.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Framework/interface/InputSourceDescription.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Sources/interface/RawInputSource.h"
#include "FWCore/Utilities/interface/GetPassID.h"
#include "FWCore/Utilities/interface/TypeWithDict.h"
#include "FWCore/Version/interface/GetReleaseVersion.h"


class CTP7SpyDAQInputSource: public edm::RawInputSource {

  public:
    explicit CTP7SpyDAQInputSource(edm::ParameterSet const&,edm::InputSourceDescription const&);
    virtual ~CTP7SpyDAQInputSource();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  protected:
    virtual bool checkNextEvent() override;
    virtual void read(edm::EventPrincipal& eventPrincipal) override;

  private:
    virtual void preForkReleaseResources() override;
    virtual void postForkReacquireResources(std::shared_ptr<edm::multicore::MessageReceiverForSource>) override;
    edm::Timestamp getTimestamp();
    edm::BranchDescription getBranchDescription(const edm::TypeWithDict& type, std::string productInstanceName = "") const;

    edm::ProcessHistoryID processHistoryID_;
    edm::BranchDescription ecalDigiBranchDescr_;
    edm::BranchDescription hcalDigiBranchDescr_;
    edm::BranchDescription caloTowerBranchDescr_;

    edm::RunNumber_t runNumber_;

    unsigned long long iCurrentLumiSection_;
    unsigned long long iCurrentEvent_;
    unsigned long long maxEvents_;
};

CTP7SpyDAQInputSource::CTP7SpyDAQInputSource(edm::ParameterSet const& pset,
                                             edm::InputSourceDescription const& desc) :
  edm::RawInputSource(pset, desc),
  processHistoryID_(),
  iCurrentEvent_(1),
  maxEvents_(pset.getUntrackedParameter<unsigned long long>("maxEvents", 10))
{
  // We have to make up provenance for the 3 collections that we save
  // This is entirely cribbed from DaqProvenanceHelper
  // I'd love to just use that class but we want to save 3 collections
  // not just one fedRawDataCollection :<
  edm::ProcessHistory ph;
  edm::ParameterSet processParameterSet_;

  std::string const& moduleLabel = "source";
  std::string const& processName = "SPY";
  std::string const& moduleName = "CTP7SpyDAQInputSource";

  processParameterSet_.addParameter("@all_sources", std::vector<std::string>({moduleLabel}));

  typedef std::vector<std::string> vstring;
  vstring empty;

  edm::ParameterSet triggerPaths;
  triggerPaths.addParameter<vstring>("@trigger_paths", empty);
  processParameterSet_.addParameter<edm::ParameterSet>("@trigger_paths", triggerPaths);

  edm::ParameterSet pseudoInput;
  pseudoInput.addParameter<std::string>("@module_edm_type", "Source");
  pseudoInput.addParameter<std::string>("@module_label", moduleLabel);
  pseudoInput.addParameter<std::string>("@module_type", moduleName);
  processParameterSet_.addParameter<edm::ParameterSet>(moduleLabel, pseudoInput);

  processParameterSet_.addParameter<vstring>("@all_esmodules", empty);
  processParameterSet_.addParameter<vstring>("@all_esprefers", empty);
  processParameterSet_.addParameter<vstring>("@all_essources", empty);
  processParameterSet_.addParameter<vstring>("@all_loopers", empty);
  processParameterSet_.addParameter<vstring>("@all_modules", empty);
  processParameterSet_.addParameter<vstring>("@end_paths", empty);
  processParameterSet_.addParameter<vstring>("@paths", empty);
  processParameterSet_.addParameter<std::string>("@process_name", processName);
  // Now we register the process parameter set.
  processParameterSet_.registerIt();
  ph.emplace_back(processName, processParameterSet_.id(), edm::getReleaseVersion(), edm::getPassID());
  processHistoryID_ = ph.setProcessHistoryID();
  processHistoryRegistryForUpdate().registerProcessHistory(ph);

  // Now set up branch descriptions for each product

  ecalDigiBranchDescr_ = getBranchDescription(edm::TypeWithDict(typeid(EcalTrigPrimDigiCollection)), "ecalDigis");
  productRegistryUpdate().copyProduct(ecalDigiBranchDescr_);

  hcalDigiBranchDescr_ = getBranchDescription(edm::TypeWithDict(typeid(HcalTrigPrimDigiCollection)), "hcalDigis");
  productRegistryUpdate().copyProduct(hcalDigiBranchDescr_);

  caloTowerBranchDescr_ = getBranchDescription(edm::TypeWithDict(typeid(l1t::CaloTowerBxCollection)), "caloTowers");
  productRegistryUpdate().copyProduct(caloTowerBranchDescr_);


  // TODO: Attach to Layer1 and get run number?
  runNumber_ = 1L;
  setNewRun();
  setRunAuxiliary(new edm::RunAuxiliary(runNumber_, edm::Timestamp::beginOfTime(), edm::Timestamp::invalidTimestamp()));

  runAuxiliary()->setProcessHistoryID(processHistoryID_);
}

CTP7SpyDAQInputSource::~CTP7SpyDAQInputSource()
{
}

void CTP7SpyDAQInputSource::read(edm::EventPrincipal& eventPrincipal)
{
  std::cout << "Reading... " << std::endl;
  edm::Timestamp tstamp = getTimestamp();

  auto eventID_ = edm::EventID(runNumber_, iCurrentLumiSection_, iCurrentEvent_);
  // FYI:
  // enum   ExperimentType {
  //     Undefined = 0, PhysicsTrigger = 1, CalibrationTrigger = 2, RandomTrigger = 3,
  //     Reserved = 4, TracedEvent = 5, TestTrigger = 6, ErrorTrigger = 15
  //     }
  edm::EventAuxiliary aux(eventID_, processGUID(), tstamp, true, edm::EventAuxiliary::TestTrigger);
  aux.setProcessHistoryID(processHistoryID_);
  std::cout << "Making event... " << std::endl;
  makeEvent(eventPrincipal, aux);

  std::cout << "Making collections... " << std::endl;
  std::unique_ptr<EcalTrigPrimDigiCollection> ecalDigis(new EcalTrigPrimDigiCollection);
  std::unique_ptr<HcalTrigPrimDigiCollection> hcalDigis(new HcalTrigPrimDigiCollection);
  std::unique_ptr<l1t::CaloTowerBxCollection> caloTowers(new l1t::CaloTowerBxCollection);

  //
  // Fill!
  //
  uint32_t cardSide{1}, absCaloEta{3}, caloPhi{5};
  uint32_t towerDatum{0x02};
  EcalTriggerPrimitiveSample sample(towerDatum); 
  EcalTrigTowerDetId id(cardSide, EcalSubdetector::EcalTriggerTower, absCaloEta, caloPhi);
  EcalTriggerPrimitiveDigi etpg(id);
  etpg.setSize(1);
  etpg.setSample(0, sample);
  ecalDigis->push_back(etpg);

  //
  // Put in our new event
  //
  std::cout << "Putting... " << std::endl;
  std::unique_ptr<edm::WrapperBase> edp(new edm::Wrapper<EcalTrigPrimDigiCollection>(std::move(ecalDigis)));
  eventPrincipal.put(ecalDigiBranchDescr_, std::move(edp), edm::ProductProvenance(ecalDigiBranchDescr_.branchID()));

  std::cout << "Putting... " << std::endl;
  std::unique_ptr<edm::WrapperBase> hdp(new edm::Wrapper<HcalTrigPrimDigiCollection>(std::move(hcalDigis)));
  eventPrincipal.put(hcalDigiBranchDescr_, std::move(hdp), edm::ProductProvenance(hcalDigiBranchDescr_.branchID()));

  std::cout << "Putting... " << std::endl;
  std::unique_ptr<edm::WrapperBase> tdp(new edm::Wrapper<l1t::CaloTowerBxCollection>(std::move(caloTowers)));
  eventPrincipal.put(caloTowerBranchDescr_, std::move(tdp), edm::ProductProvenance(caloTowerBranchDescr_.branchID()));

  iCurrentEvent_++;
}

bool CTP7SpyDAQInputSource::checkNextEvent()
{
  // TODO: Possible to guess lumi based on info from CTP7?
  iCurrentLumiSection_ = 1;
  if (!luminosityBlockAuxiliary() || luminosityBlockAuxiliary()->luminosityBlock() != iCurrentLumiSection_) {
    resetLuminosityBlockAuxiliary();
    edm::LuminosityBlockAuxiliary* lumiBlockAuxiliary =
      new edm::LuminosityBlockAuxiliary(
        runAuxiliary()->run(),
        iCurrentLumiSection_, getTimestamp(),
        edm::Timestamp::invalidTimestamp()
      );
    setLuminosityBlockAuxiliary(lumiBlockAuxiliary);
    luminosityBlockAuxiliary()->setProcessHistoryID(processHistoryID_);
  }

  // TODO: Return false when no events left (run stopped, capture failed, etc.)
  if ( maxEvents_ != 0 && iCurrentEvent_ > maxEvents_ ) return false;

  setEventCached();
  return true;
}

// I would be lying if I said I knew what these were for...
void CTP7SpyDAQInputSource::preForkReleaseResources()
{}

void CTP7SpyDAQInputSource::postForkReacquireResources(std::shared_ptr<edm::multicore::MessageReceiverForSource>)
{
  setRunAuxiliary( new edm::RunAuxiliary(runNumber_, edm::Timestamp::beginOfTime(), edm::Timestamp::invalidTimestamp()));
}

void
CTP7SpyDAQInputSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Input Source for reading CTP7 card ethernet readout");
  desc.addUntracked<unsigned long long> ("maxEvents", 10)->setComment("How many events before we bail out");
  desc.setAllowAnything();
  descriptions.add("source", desc);
}

edm::Timestamp
CTP7SpyDAQInputSource::getTimestamp() {
  edm::TimeValue_t time;
  timeval stv;
  gettimeofday(&stv,0);
  time = stv.tv_sec;
  time = (time << 32) + stv.tv_usec;
  edm::Timestamp tstamp(time);
  return tstamp;
}

edm::BranchDescription
CTP7SpyDAQInputSource::getBranchDescription(const edm::TypeWithDict& type, std::string productInstanceName) const {
  return edm::BranchDescription(edm::InEvent
                  , "source"
                  , "SPY"
                  , type.userClassName()
                  , type.friendlyClassName()
                  , productInstanceName
                  , "CTP7SpyDAQInputSource"
                  , edm::ParameterSetID()
                  , type
                  , false);
}

#include "FWCore/Framework/interface/InputSourceMacros.h"
DEFINE_FWK_INPUT_SOURCE(CTP7SpyDAQInputSource);
