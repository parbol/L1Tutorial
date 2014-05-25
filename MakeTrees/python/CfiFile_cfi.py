import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('MakeTrees',
  outputFileName        = cms.string("tree.root")
)
