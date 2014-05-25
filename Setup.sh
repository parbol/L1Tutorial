#!/bin/bash

#####
##### Runing UCT2015 + HLT (Setup for CMSSW_7_0_0)
#####


#Step 1. Create release area and get the code
cmsrel CMSSW_7_0_0_pre9
cd CMSSW_7_0_0_pre9/src/
cmsenv

git cms-addpkg DataFormats/L1CaloTrigger
git cms-addpkg L1TriggerConfig/L1ScalesProducers
git cms-addpkg L1Trigger/RegionalCaloTrigger     
git clone https://github.com/uwcms/UCT2015.git L1Trigger/UCT2015
cd L1Trigger/UCT2015/
git checkout uct2015Core
cd ../..
scramv1 b -j8


#Step 2. Getting a copy of the HLT menu
edmConfigFromDB --cff --configName /users/tropiano/2013/iterative_tracking/PVconstraint/hltMenu8e33/V2 --nopaths > setup_cff.py
hltGetConfiguration /users/tropiano/2013/iterative_tracking/PVconstraint/hltMenu8e33/V2 --full --offline --mc --unprescale --process TEST --globaltag auto:hltonline > hlt.py


#Step 3. Add the following lines to hlt.py (just below the definition of the process name)
#process.load("setup_cff")
#process.load("L1Trigger.UCT2015.emulationMC_cfi")
#process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")


#Step 4. Replace
#process.HLTL1UnpackerSequence = cms.Sequence( process.gtUCTDigis + process.gctUCTDigis  + process.hltL1GtObjectMap + process.l1extraParticlesUCT )
#by
#process.HLTL1UnpackerSequence=cms.Sequence(process.emulationSequence *  process.uct2015L1Extra)


#Step 5. Replace
#hltL1extraParticles by l1extraParticlesUCT
#hltGctDigis by gctUCTDigis
#hltGtDigis by gtUCTDigis 


#Step 6. Replace the L1seed module you are interested in:

#process.hltL1sL1SingleJet36 = cms.EDFilter( "HLTLevel1GTSeed",
#    L1SeedsLogicalExpression = cms.string( "L1_SingleJet36" ),
#    saveTags = cms.bool( True ),
#    L1MuonCollectionTag = cms.InputTag( "l1extraParticlesUCT" ),
#    L1UseL1TriggerObjectMaps = cms.bool( True ),
#    L1UseAliasesForSeeding = cms.bool( True ),
#    L1GtReadoutRecordTag = cms.InputTag( "gtUCTDigis" ),
#    L1CollectionsTag = cms.InputTag( "l1extraParticlesUCT" ),
#    L1NrBxInEvent = cms.int32( 3 ),
#    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
#    L1TechTriggerSeeding = cms.bool( False )
#)

#by

#process.hltL1sL1SingleJet36 = cms.EDFilter( "HLTLevel1Jet",
#    inputTag = cms.InputTag( "l1extraParticlesUCT","Central" ),
#    MinE=cms.double(-1),
#    MinPt=cms.double(16),
#    MinN=cms.int32(1),
#    MinMass=cms.double(-1),
#    MaxEta=cms.double(-1), 
#    triggerType=cms.int32(-84)
#)



