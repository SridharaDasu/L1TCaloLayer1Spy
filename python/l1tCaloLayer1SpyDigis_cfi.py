import FWCore.ParameterSet.Config as cms

l1tCaloLayer1SpyDigis = cms.EDProducer('L1TCaloLayer1Spy',
                                       Layer1PhiMapXMLFile = cms.untracked.string("CTP7phiMap.xml"),
                                       SelectedOrbitNumber = cms.untracked.uint32(1000),
                                       SelectedBXNumber = cms.untracked.uint32(1000),
                                       verbose = cms.untracked.bool(False)
                                       )
