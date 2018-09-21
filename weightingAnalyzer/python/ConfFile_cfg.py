import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(# Input file to be analyzed.
#        'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/CITo2Mu_M300_CUETP8M1_Lam22TeVConLL_13TeV_Pythia8_Corrected-v3/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/62DF9A66-1B3B-E711-AFFB-3417EBE5291B.root',

#        'file:/uscms/home/qlam/nobackup/CMSSW_9_4_0/src/MC17miniAOD.root'

#        'file:/afs/cern.ch/work/s/szaleski/private/prakashTest/CMSSW_8_0_21/src/prundedGen/weightingAnalyzer/python/step_MINIAODSIM_bumpHunt_99.root'
#         'root://cmsxrootd.fnal.gov//store/group/lpcci2dileptons/CITo2Mu_GENSIM/MuMu_GENSIM/170609_183422/0000/EXO-RunIISummer15GS-09252_1.root',
#        'root://cmsxrootd.fnal.gov///store/user/szaleski/CITo2E_GENSIM_Lam10/EE_GENSIM_LLConM2000/170814_202655/0000/EXO-RunIISummer15GS-09252_20.root'

        #Contact interaction MC17 CIto2Mu M300to800 16 TeV ConLR 

#        'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/3049DF8B-E866-E811-A2FC-3417EBE6449E.root',

#       'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/1CFF8E58-E866-E811-AC46-5065F38182E1.root',

#       'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/5CDB7456-E866-E811-8E3A-008CFAEBDEEC.root',

#       'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/266ADD49-E866-E811-B54D-7CD30ACDCACA.root',

#       'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/B4722E1A-0B69-E811-A80E-E0071B7A8560.root',

#       'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/F603BF8A-436B-E811-8A9E-0025905B8564.root',

        # Drell_Yan MC17 MC17 DYtoMu M300to800

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M300to800_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/2C6E8724-4472-E811-86C0-A0369FD1EF00.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M300to800_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/68B576BC-4372-E811-A262-AC1F6B8DD2A0.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M300to800_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/B23E8205-4472-E811-B77C-A4BF01125E26.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M300to800_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/8EA9ED25-4472-E811-AA3F-20CF300E9EB7.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M300to800_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/26822A14-4472-E811-9B5E-801844DEF358.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M300to800_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/00000/B2F984A4-DA8D-E811-B46F-509A4C858BCB.root',

# CI MC17 CiTo2E M2000toInf 16 TeV ConLR
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/F4AA40CD-9163-E811-B28B-0025905C3E68.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/10000/5643F087-2264-E811-8451-24BE05CEFB41.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/98216DB8-9163-E811-9EF0-90B11C4442A1.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/381E3BB2-9163-E811-9F8D-24BE05CEBD61.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/B06A1F7F-EC65-E811-BBD6-002481DE49B6.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/7E342B16-3F65-E811-A9B9-0CC47A4C8ECA.root',

# MC17 DYtoMu M2000toInf 16 TeV ConLR

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M2000toInf_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/28ED272D-6663-E811-A4DC-F01FAFE0F396.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M2000toInf_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/A02EDB0B-4863-E811-AE2D-F01FAFE5FB02.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M2000toInf_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/DCA85C6D-7263-E811-8013-E0071B7AE500.root',

# CI MC17 CItoMu M2000toInf 16 TeV ConLR
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/10000/889C9871-F163-E811-B225-FA163EEFE78F.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/60F39FB9-E762-E811-B0D5-FA163E8A1979.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/14BE8087-E762-E811-9EC5-44A84225C4EB.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/F2945632-D562-E811-B050-B8CA3A70A5E8.root',

# CI MC17 CIToMu M2000toInf 16TeV DesLR
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/1402547E-4B5A-E811-A9B5-0025905C54DA.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/8A33307D-8958-E811-995B-0025905C3E66.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/50B12152-4B5A-E811-9570-0CC47A4D9A72.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/C844520E-625A-E811-9B78-0CC47A1E0DC2.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/1810C325-4B5A-E811-A75A-24B6FDFF17A7.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/B2A4F669-905A-E811-8E15-90B11C2AA16C.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/6EE4B94F-4B5A-E811-858C-0CC47A4D76D6.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M2000toInf_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/30000/06C39215-FE5A-E811-8857-D8D385AF883C.root',


# MC17 CIToMu M1300to2000 16TeV ConLR

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M1300to2000_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/CA743409-3D98-E811-943C-0CC47A6C063E.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M1300to2000_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/66B4BD17-3D98-E811-85A0-0026B92779F1.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M1300to2000_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/A81F5005-3D98-E811-A6B4-0CC47A7E6A6C.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M1300to2000_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/12E06E94-3C98-E811-8045-5065F37D50E2.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M1300to2000_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/40000/607642B8-D1A9-E811-B02A-0CC47A4D76C0.root',
 
# MC17 DYToMu M1300to2000 16TeV ConLR
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M1300to2000_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/0A0D1B5B-5664-E811-A38C-B4969109F628.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M1300to2000_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/962146FB-A364-E811-8084-008CFAF554E6.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M1300to2000_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/EAD43262-6264-E811-8BA7-E0071B7AC760.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M1300to2000_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/42974AB2-4E64-E811-BA96-90B11C4434EB.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M1300to2000_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/1AE57097-2066-E811-AA6F-008CFAF3543C.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M1300to2000_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/F69C9367-9F65-E811-BF52-509A4C63A055.root',


# MC17 CIToE M800to1300 16TeV ConLR

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M800to1300_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/AA6FC652-ED5E-E811-8D59-842B2B5C2299.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M800to1300_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/268B152B-ED5E-E811-B251-0CC47A78A468.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M800to1300_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/9C9A2E3C-ED5E-E811-A4C9-0CC47A57D136.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M800to1300_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/000E3752-ED5E-E811-8924-008CFA56D770.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M800to1300_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/D0E6D5A3-9A63-E811-9980-E0071B7B2320.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2E_M800to1300_CP5_Lam16TeVConLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/926A146C-6263-E811-9DE9-EC0D9A8222DE.root',

# MC17 DYToMu M800to1300 16TeV ConLR

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M800to1300_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/E2BBA9BB-6363-E811-885C-1866DA85DFF0.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M800to1300_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/F48E4C43-E563-E811-B272-003048FFD76C.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M800to1300_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/A2BD46B7-6363-E811-860C-EC0D9A82261E.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_DYTo2Mu_M800to1300_CP5_Pythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/20000/78E74850-CA61-E811-A1A9-E0071B73B6B0.root',


# MC17 CIToE M300to800 16TeV ConLR

# MC17 CIToMu M300to800 16TeV DesLR

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/C82F4DAF-C759-E811-A91A-D48564593F64.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/B8B53FB5-C759-E811-8671-002590D9D8AE.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/B453CBA3-B859-E811-8181-5065F3815221.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/4E38299E-2359-E811-ADCA-0CC47AFB7DA8.root',

#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/MC17_CITo2Mu_M300to800_CP5_Lam16TeVDesLRPythia8_v1/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/B4A908AC-5760-E811-A39D-0242AC130002.root',

)
                            )


from GenStudy.weightingAnalyzer.CfiFile_cfi import *
process.weightingAnalyzer=weighting_Analyzer.clone()

process.TFileService = cms.Service("TFileService",


 

                                   fileName = cms.string("Pythia8_Sep12_DYToMu_16TeV_CIM1300to2000_ConLR.root")#  


)


process.p = cms.Path(process.weightingAnalyzer)

