import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("L1TCaloLayer1LUTWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('L1Trigger.L1TCaloLayer1Spy.l1tCaloLayer1LUTWriter_cfi')

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '74X_dataRun2_Express_v1'
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource', 'GlobalTag')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/data/dasu/0EFD41DE-866B-E511-9644-02163E0143CE.root')
)

# Writes LUT for the only event to be processed - ignores data itself.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.p = cms.Path(process.l1tCaloLayer1LUTWriter)
