import FWCore.ParameterSet.Config as cms

l1tCaloLayer1LUTWriter = cms.EDAnalyzer('L1TCaloLayer1LUTWriter',
    verbose = cms.untracked.bool(False),
    fileName = cms.string("luts.xml"),
    useLSB = cms.bool(True),
    useCalib = cms.bool(True),
    useECALLUT = cms.bool(True),
    useHCALLUT = cms.bool(True),
    useHFLUT = cms.bool(True),
    firmwareVersion = cms.int32(1),
)
