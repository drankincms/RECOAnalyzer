import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_ICHEP16JEC_v0"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/024C76E8-6857-E611-B63D-0025905A612A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/02C39237-6357-E611-9051-0025905B85AA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/0A249723-4157-E611-9B57-0CC47A4D75F2.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/0E14F83E-6C57-E611-9C95-0CC47A4D7618.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/0ECD89B2-4457-E611-819D-0CC47A4D765E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/10F0F383-2057-E611-8752-0CC47A78A3F4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/16264747-2857-E611-998D-0025905A60EE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/1C95321E-4157-E611-8D7C-0CC47A4D7662.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/1CA35366-3357-E611-9562-0CC47A4D7614.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/1CB4FA0B-2B57-E611-A1D8-0CC47A4C8F2C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/1EB23B68-3357-E611-A392-0CC47A4D76BE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/1EEE33B9-3A57-E611-8B8C-0CC47A4C8E28.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/225487E2-2157-E611-BE3E-0025905A605E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/244E84C3-6F57-E611-A417-0025905B8590.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/262B4F34-2857-E611-9828-0025905A6082.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/2A568E99-2C57-E611-B1AB-0025905B85EE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/2AC56122-7057-E611-BB2B-0025905A48E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/2EAAF142-2C57-E611-BF1F-0CC47A4D761A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/2EB9A349-6957-E611-A44F-0CC47A4D76D2.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/3022E021-3057-E611-9AFE-0025905B85C6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/34BCC3FC-3C57-E611-A3FD-0025905B857E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/3661D070-4757-E611-8734-0025905B8606.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/366628C6-4457-E611-A2C3-0025905A60F8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/3803DA6A-4A57-E611-8876-0CC47A4D7638.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/38CACB6B-4A57-E611-80CE-0CC47A4D769C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/3ABF08B7-4357-E611-85C8-0CC47A78A33E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/3CD8280C-6357-E611-94CF-0025905A610C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4073958A-2E57-E611-AE8C-0CC47A78A4A0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/40DE5E64-3957-E611-B656-0CC47A4D75F6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/44C1E910-3457-E611-9780-0CC47A4C8F30.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/485DEA6B-3857-E611-8A70-0CC47A4C8E3C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4A0FF08C-6C57-E611-A87B-0CC47A4D762E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4A7B0E35-4957-E611-8DA5-0025905B85C0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4C03F866-7157-E611-9EF6-0025905B8564.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4C37648A-3057-E611-B7FD-0CC47A4D75EC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4C931299-6257-E611-B5F0-0025905A60A8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/4CC8E8C1-2457-E611-B3D5-0CC47A4D765E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/5009673C-7357-E611-B160-0025905A60C6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/5248E0C3-3157-E611-8C7A-0CC47A4C8F1C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/5270BA73-7557-E611-B138-0025905A48EC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/580F7436-2C57-E611-B8C9-0CC47A4D7664.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/58BF53DD-3657-E611-82FC-0025905A48E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/5A99CED9-2457-E611-98B4-0CC47A4C8EB6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/604BFAFC-3B57-E611-98A6-0025905B857E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/64E5D090-2957-E611-B774-0025905A497A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/68090693-3D57-E611-B72C-0025905A48FC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/68E6FBBC-2257-E611-912E-0025905B8596.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/6A00D768-3957-E611-87DA-0CC47A4D765E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/6A2AA240-4957-E611-BAA1-0025905B85DA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/6A5A9D3A-2E57-E611-9C32-0025905B85B6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/6E3D2CE4-6857-E611-A283-0CC47A4D76A2.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/6E4F904C-3A57-E611-9218-0CC47A78A3B4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/7035ED0E-2957-E611-A117-0CC47A4D7632.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/7418EACA-4257-E611-8A18-0CC47A4D7674.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/78190B98-6357-E611-AA09-0025905A610C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/7AD160EC-6857-E611-8E81-0CC47A4D767C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/7E6C01E5-3F57-E611-A7AF-0025905B858A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/805D017E-6D57-E611-93B4-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/806517CA-3A57-E611-8197-0CC47A4C8F10.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/84110BF8-2757-E611-8B5B-0CC47A745298.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/8CB37882-6C57-E611-9DB2-0CC47A78A436.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/8E54A04E-3557-E611-BCF7-0CC47A4D76CC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/901C0F4C-2C57-E611-AFAC-0025905B85DE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/98AF2ED8-3F57-E611-A919-0CC47A4D7632.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/9CC8EA88-3357-E611-9C80-0025905A613C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/9CE563C3-3A57-E611-BC89-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/A6375AD2-4257-E611-8255-0CC47A4D7644.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/A6B9AF8A-2B57-E611-9502-0CC47A4C8EE8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/A829B61B-3B57-E611-9467-0025905B857E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/AA3B3411-3B57-E611-A626-0CC47A78A456.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/AC69583C-4857-E611-95A0-0CC47A4C8ECE.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/AE720405-2857-E611-9978-0025905A48EC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/B059C47C-6B57-E611-8147-0CC47A4D7638.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/B493DE44-4657-E611-AA9A-0025905B860C.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/BA770802-7057-E611-84A5-0CC47A4D7662.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/C69AC578-2657-E611-8A69-0025905B858E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/C8D2EC73-2957-E611-B6F0-0CC47A4C8E1E.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/D07516B8-6F57-E611-BC3E-0CC47A4C8E98.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/D629E11A-2957-E611-89EF-0CC47A78A4B8.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/DC235FBD-2157-E611-B869-0CC47A78A3B4.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/DCCF2EFA-2757-E611-B884-0CC47A78A2EC.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/DEB9C381-2957-E611-B590-0025905B85C6.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/E04B7E03-4557-E611-8041-0025905A60E0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/E208A83A-4857-E611-9AD4-0CC47A4D7606.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/E4853A9A-2257-E611-809D-0CC47A4C8E8A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/E8829A8B-3957-E611-819E-0025905B85AA.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/E8C89E84-3C57-E611-8342-0025905A497A.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/F49F4B3B-2857-E611-A8E8-0CC47A4D76D0.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/F81BAA21-2757-E611-94F1-0025905A6084.root',
       '/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_81X_mcRun2_asymptotic_v2_2023tiltedPU200-v1/10000/FA0D822B-7457-E611-A37F-0CC47A4C8F26.root'
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
    fileName = cms.string("histo_ttbar_pu200.root"),
    closeFileFast = cms.untracked.bool(True)
)

#process.p = cms.Path(process.ak4PFCHSL1FastL2L3CorrectorChain + process.ak4PFchsCorrectedJets + process.check)
process.p = cms.Path(process.check)
