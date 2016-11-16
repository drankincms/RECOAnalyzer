import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_ICHEP16JEC_v0"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/02C87E40-6B55-E611-9C28-0025905A48C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/0650DB43-6B55-E611-B9F5-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/067773B7-6A55-E611-B4CE-0CC47A4C8F10.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/0850F677-6A55-E611-97B7-0025905B85DC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/08813C41-6B55-E611-A3B4-0CC47A4C8EB6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/0A07A379-6C55-E611-912E-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/0A90A46B-6855-E611-8523-0CC47A4D762E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/0C2C9319-6A55-E611-BFD9-0025905B8592.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/108CACCD-6A55-E611-89F7-0025905A612C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/16A68198-6955-E611-BD77-0025905B85DC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/1E7AD0FD-6A55-E611-ABF5-0CC47A78A4BA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/205390F6-6855-E611-840B-0CC47A78A42C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/240A8AB5-6D55-E611-B4F0-0025905B85DC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/2C388CB1-6B55-E611-96F2-0CC47A78A3D8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/34758B8D-6955-E611-97FA-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/3A5BCECB-6A55-E611-B02C-0CC47A7452D8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/3C6C7B6C-6955-E611-8CED-0CC47A78A3D8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/3E63A615-6A55-E611-8F94-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/56219144-6B55-E611-97A9-0025905B856E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/5A294E7D-6955-E611-93D0-0CC47A78A3D8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/64D4B573-6A55-E611-9B62-0025905A48C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/66CB5746-6B55-E611-947B-0025905B85CA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/6AC3389B-6B55-E611-AACA-0CC47A4D764A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/6C4A43A8-6B55-E611-B377-0CC47A4C8F10.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/76EF7F72-6A55-E611-A980-0CC47A4C8F10.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/7EF1A135-7155-E611-9535-0025905A612C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/82421FB4-6B55-E611-9275-0CC47A78A45A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/841D3270-6A55-E611-B871-0CC47A4D764A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/8CE59941-6B55-E611-AE57-0025905A60BE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/964050BD-6A55-E611-8106-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/967B3542-6B55-E611-8984-0CC47A78A440.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/9AAC2316-6A55-E611-9919-0CC47A4D762E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/9EE869CE-7155-E611-9ADE-0025905B8592.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/AE527CF5-6855-E611-B462-0025905A48C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/B2598F95-6955-E611-B083-0025905A612C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/B2CE3BCD-6A55-E611-A022-0CC47A4D764A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/B4B6B941-6B55-E611-BCDF-0CC47A78A42C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/B6A0C8F8-6855-E611-842B-0025905B856E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/BCC7AC54-6A55-E611-B6F0-0CC47A78A3D8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/D615FF19-7255-E611-ADAA-0CC47A78A45A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/DAB9FD47-6B55-E611-9414-0025905B85DC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/DCE4857B-6955-E611-BBEB-0CC47A78A45A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/E4DBBAFE-6B55-E611-BB1C-0CC47A78A3D8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/EC112F43-6B55-E611-936F-0025905B8592.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/F07CD074-6A55-E611-A2D9-0025905B856E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU35-v1/10000/F8F5C215-6A55-E611-B333-0CC47A4D76B8.root'
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
    met = cms.InputTag("pfMet"),
    genparticles = cms.InputTag("genParticles"),
    electrons = cms.InputTag("gedGsfElectrons"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("particleFlowEGamma"),
    debug = cms.bool(False)
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("histo_ttbar_pu35.root"),
    closeFileFast = cms.untracked.bool(True)
)

#process.p = cms.Path(process.ak4PFCHSL1FastL2L3CorrectorChain + process.ak4PFchsCorrectedJets + process.check)
process.p = cms.Path(process.check)
