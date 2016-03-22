import FWCore.ParameterSet.Config as cms

process = cms.Process("RUNSPY")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("CTP7SpyDAQInputSource",
    maxEvents = cms.untracked.uint64(10),
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(process.content)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('spyOutput.root'),
)

process.e = cms.EndPath(process.out)
