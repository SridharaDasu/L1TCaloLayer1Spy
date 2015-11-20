import FWCore.ParameterSet.Config as cms

l1tCaloLayer1LUTWriter = cms.EDAnalyzer('L1TCaloLayer1LUTWriter',
                                        verbose = cms.untracked.bool(False)
                                        )
