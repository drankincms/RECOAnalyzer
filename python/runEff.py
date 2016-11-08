import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_ICHEP16JEC_v0"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_8_1_0_pre9/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/FA5F7E7E-7055-E611-8095-0025905B856E.root'
    )
)

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

process.ak4PFchsCorrectedJets   = cms.EDProducer('CorrectedPFJetProducer',
    src         = cms.InputTag('ak4PFJetsCHS'),
    correctors  = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')
)

process.load("PhysicsTools/JetMCAlgos/HadronAndPartonSelector_cfi")  
process.load("PhysicsTools/JetMCAlgos/AK4PFJetsMCFlavourInfos_cfi")
#process.ak4JetFlavourInfos.jets = cms.InputTag("ak4PFchsCorrectedJets")
process.ak4JetFlavourInfos.jets = cms.InputTag("ak4PFJetsCHS")

process.check = cms.EDAnalyzer("RECOAnalyzer",
    vertices = cms.InputTag("offlinePrimaryVertices"),
    #jets = cms.InputTag("ak4PFchsCorrectedJets"),
    jets = cms.InputTag("ak4PFJetsCHS"),
    btags = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    flavorinfo = cms.InputTag("ak4JetFlavourInfos"),
    genparticles = cms.InputTag("genParticles"),
    electrons = cms.InputTag("gedGsfElectrons"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("particleFlowEGamma"),
    muons = cms.InputTag("muons"),
    debug = cms.bool(True)
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("histo.root"),
    closeFileFast = cms.untracked.bool(True)
)

#process.p = cms.Path(process.ak4PFCHSL1FastL2L3CorrectorChain + process.ak4PFchsCorrectedJets + process.selectedHadronsAndPartons + process.ak4JetFlavourInfos + process.check)
process.p = cms.Path(process.selectedHadronsAndPartons + process.ak4JetFlavourInfos + process.check)
