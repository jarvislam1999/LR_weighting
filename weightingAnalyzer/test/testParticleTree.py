import FWCore.ParameterSet.Config as cms

process = cms.Process("testParticle")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/p/pthapa/private/miniAOD_2016/CMSSW_8_0_21/src/GenStudy/weightingAnalyzer/test/RootFilesWeighting/RL13_Con_110_M500_Lamda10000_Test_Tree_Structure_MINIAOD.root')
)

process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("prunedGenParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
)

process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("prunedGenParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex  = cms.untracked.bool(False)
)

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.printTree1*process.printTree2)
process.outpath = cms.EndPath(process.printEventNumber)
process.MessageLogger.destinations = cms.untracked.vstring('cout','cerr')
