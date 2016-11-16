import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_ICHEP16JEC_v0"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/02DFF30C-1C56-E611-88C0-0025905A60B6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/046F571B-1856-E611-988B-0CC47A4C8E66.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/087CF19F-1A56-E611-A310-003048FFD72C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/0A589E89-1956-E611-BE61-0CC47A78A2F6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/0E9546A3-1B56-E611-8477-0CC47A7452D0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/10D95652-1B56-E611-9140-0CC47A78A468.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/10E26870-1956-E611-AE01-0025905B85DE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/1C081C70-1B56-E611-8FA8-0025905B8568.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/24C6246C-1956-E611-8C7B-0025905A60D6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/28DAACE3-1A56-E611-8BFD-0025905A6088.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/2AAC2C9A-1A56-E611-BD9F-0025905B85C6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/2AAFE0DC-1B56-E611-BF23-0025905B85A2.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/2AC93E4D-1A56-E611-91AD-0025905A6092.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/305D6B4A-1A56-E611-8A36-0025905B8594.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/306A6616-1956-E611-ABAD-0025905B85EE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/30DEACC8-1A56-E611-8DC4-0025905A6090.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/346E82C2-1B56-E611-8F91-0CC47A78A3F4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/34FA4D0B-1C56-E611-AC79-0025905B8590.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/36119547-1B56-E611-B8B5-0025905A6088.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/36A12D51-1B56-E611-A1DC-0CC47A4D76BE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/3A499DBF-1856-E611-8809-0CC47A4C8F0A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/3AE33404-1B56-E611-966C-0CC47A4D76AA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/3EB99A56-1A56-E611-996B-0025905A60B0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/4017D265-1956-E611-8D2D-0025905B859A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/44DAD8DC-1A56-E611-A5E0-0025905AA9F0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/468D3732-1B56-E611-924F-0025905A6076.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/48BE86DA-1856-E611-B05E-0025905A60BC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/4E3FDBA3-1956-E611-9CE7-0025905B8612.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/52199A18-1C56-E611-B86E-0025905A60B8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/52A4948A-1A56-E611-A6CA-0025905A60DE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/52FC1731-1B56-E611-B854-0025905A607E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/54F4234F-1A56-E611-A246-0CC47A4D761A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/581A0C41-1B56-E611-82D8-0025905A60D0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/58DC2632-1B56-E611-9ACD-0CC47A78A2F6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/5CA9579A-1A56-E611-89B5-0025905A60EE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/5CEF08B1-1B56-E611-A4C8-0CC47A4C8E96.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/5EF5945B-1956-E611-A9F6-0025905A48FC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/60734201-1956-E611-A072-0025905B857E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/628BB106-1C56-E611-9799-0025905B85C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/6417DB96-1A56-E611-AACE-003048FFD79E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/6821E309-1956-E611-BE5C-0025905A48D0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/7082F1E9-1B56-E611-8F43-0CC47A4C8E28.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/70B59E2D-1B56-E611-AC85-0025905B85B8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/74ED6511-1B56-E611-8E7E-0CC47A745298.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/781E1E47-1B56-E611-96D6-0CC47A4C8EE2.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/78D38C7D-1B56-E611-BFD7-0025905B859A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/7A096145-1A56-E611-AF8D-003048FFD7A4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/7C0A84A7-1956-E611-AC10-0025905B85E8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/7CBD8117-1C56-E611-A912-0025905A48C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/7E6FDB49-1856-E611-AF70-0CC47A78A3F4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/8095AC88-1756-E611-8855-0025905B8600.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/86FF7786-1A56-E611-9337-0025905A60D2.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/8A7C0F48-1856-E611-8E68-0CC47A4D768E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/8CC4F598-1956-E611-8CE5-0CC47A4D7686.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/904DA5A7-1956-E611-9BFD-0025905B85AA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/9200F597-1A56-E611-94A0-0CC47A78A41C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/922D6FE6-1A56-E611-B3BC-0025905A48C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/92A6016E-1956-E611-BC69-0025905B858E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/98A6600B-1C56-E611-BE3A-0025905A4964.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/98A9300B-1C56-E611-ADFA-0025905A48BA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/9A6C790C-1C56-E611-BC03-0025905A60BC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/9E3156A9-1A56-E611-9FED-0025905A60D6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/A2168CDF-1A56-E611-A635-0025905A60FE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/A24F72DE-1A56-E611-98B4-0025905A609A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/A47CA347-1A56-E611-A867-0CC47A7452D0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/A6C5D73B-1A56-E611-BF1A-0025905A60FE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/AA6DFD9B-1A56-E611-9606-0CC47A4D7690.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/B4775342-1956-E611-A046-0CC47A4C8ED8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/B8C2F243-1A56-E611-B7AE-0CC47A4D76C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/BC782D12-1B56-E611-9178-0CC47A78A33E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/BE18603D-1B56-E611-868A-0CC47A745294.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/C0C5B656-1A56-E611-A21E-0025905B8600.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/C80645CC-1856-E611-9686-0CC47A4C8EE8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/C82443DF-1B56-E611-8DB2-0025905A609A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/CADCC2E4-1856-E611-A2F0-0025905A606A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/CC5623D6-1856-E611-82FD-0025905A611E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/D042CEC4-1B56-E611-8925-0CC47A4C8E2A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/D082E492-1A56-E611-8D18-0025905A497A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/D20053BD-1A56-E611-BDA4-0CC47A4C8E2E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/D206AA0E-1C56-E611-B3C4-0025905B8580.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/D2786166-1856-E611-AD94-0CC47A4C8E1C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/D287DB2A-1B56-E611-AC8E-003048FFD7AA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/DE600232-1756-E611-9ED7-0025905B85DE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/DEAD849B-1A56-E611-B5E8-0CC47A78A33E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/DEE3FB08-1956-E611-B8F1-0025905B8564.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/E804448B-1956-E611-B905-00248C55CC97.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/EAE123CB-1A56-E611-91B8-0025905A6118.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/EC5F3CDD-1A56-E611-9FAF-0025905A60AA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/F2A5BEC7-1B56-E611-BDED-0CC47A4C8E7E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU140-v1/10000/FEDF8C8C-1A56-E611-BFE3-0025905B860C.root'
    )
)

process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

process.ak4PFchsCorrectedJets   = cms.EDProducer('CorrectedPFJetProducer',
    src         = cms.InputTag('ak4PFJetsCHS'),
    correctors  = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')
)

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
    fileName = cms.string("histo_ttbar_pu140.root"),
    closeFileFast = cms.untracked.bool(True)
)

#process.p = cms.Path(process.ak4PFCHSL1FastL2L3CorrectorChain + process.ak4PFchsCorrectedJets + process.check)
process.p = cms.Path(process.check)
