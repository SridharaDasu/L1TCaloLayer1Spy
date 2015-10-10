import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloLayer1Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('L1Trigger.L1TCaloLayer1Spy.l1tCaloLayer1SpyDigis_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("EmptySource")

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('l1tCaloLayer1Spy.root'),
    outputCommands = cms.untracked.vstring('keep *')
)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '74X_dataRun2_Express_v1'

process.p = cms.Path(process.l1tCaloLayer1SpyDigis)

process.e = cms.EndPath(process.out)
