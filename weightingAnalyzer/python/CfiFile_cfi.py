import FWCore.ParameterSet.Config as cms

weighting_Analyzer = cms.EDAnalyzer('weightingAnalyzer',
                      prunedGenParticlesTag=cms.InputTag("prunedGenParticles"),
                      genInfoProduct=cms.InputTag("generator"),
                      genEventInfoProduct=cms.InputTag("generator"),
                      decayParticlePID=cms.int32(13),
                      debug=cms.int32(0),
                                    
)


