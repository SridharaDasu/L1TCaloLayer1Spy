import FWCore.ParameterSet.Config as cms

l1tCaloLayer1SpyDigis = cms.EDProducer('L1TCaloLayer1Spy',
                                       Layer1PhiMapXMLFile = cms.untracked.string("CTP7phiMap.xml"),
                                       verbose = cms.bool(False)
                                       )
