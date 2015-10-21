import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloLayer1Spy")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('L1Trigger.L1TCaloLayer1Spy.l1tCaloLayer1SpyDigis_cfi')

# Put multiples of 162 - output data for eighteen BXs are available for each capture
# One event is created for each capture.  Putting non-multiples of 162 just means
# that some of the events captured are "wasted".

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1620) )

process.source = cms.Source("EmptySource")

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/data/dasu/l1tCaloLayer1Spy.root'),
    outputCommands = cms.untracked.vstring('keep *')
)

process.p = cms.Path(process.l1tCaloLayer1SpyDigis)

process.e = cms.EndPath(process.out)
