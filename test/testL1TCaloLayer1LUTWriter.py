import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process("L1TCaloLayer1LUTWriter",eras.Run2_2016)

options = VarParsing()
options.register('runNumber', 260627, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.parseArguments()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_HLT_v6', '')

process.source = cms.Source('EmptySource',
    firstRun = cms.untracked.uint32(options.runNumber)
)

# Writes LUT for the only event to be processed - ignores data itself.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load('L1Trigger.L1TCaloLayer1Spy.l1tCaloLayer1LUTWriter_cfi')

process.p = cms.Path(process.l1tCaloLayer1LUTWriter)

process.schedule = cms.Schedule(process.p)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW,L1TEventSetupForHF1x1TPs 

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAW(process)

#call to customisation function L1TEventSetupForHF1x1TPs imported from L1Trigger.Configuration.customiseReEmul
process = L1TEventSetupForHF1x1TPs(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
#from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAODEMU 

#call to customisation function L1NtupleAODEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
#process = L1NtupleAODEMU(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseUtils
#from L1Trigger.Configuration.customiseUtils import L1TTurnOffUnpackStage2GtGmtAndCalo 

#call to customisation function L1TTurnOffUnpackStage2GtGmtAndCalo imported from L1Trigger.Configuration.customiseUtils
#process = L1TTurnOffUnpackStage2GtGmtAndCalo(process)

# End of customisation functions

