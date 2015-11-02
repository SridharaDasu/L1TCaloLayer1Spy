import FWCore.ParameterSet.Config as cms

l1tCaloLayer1SpyDigis = cms.EDProducer('L1TCaloLayer1Spy',
                                       setupString = cms.untracked.string("ethernet:CTP7phiMap.xml"),
                                       SelectedBXNumber = cms.untracked.uint32(1000),
                                       verbose = cms.untracked.bool(False)
                                       )
