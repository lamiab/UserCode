import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIZ")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = "Zmumu_Tree_pPbData_Prompt_210658-210737_ntrk.root"
#import os,commands
#def getDirectoryList(path):
#    cmd  = 'ls %s/ ' % (path)
#    file = ["file:%s/%s" % (path,i) for i in commands.getoutput(cmd).split("\n")]
#    return file

#options.inputFiles = getDirectoryList("/afs/cern.ch/work/m/mironov/public/regitSkim147mub/skim")

options.inputFiles = "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_10_1_M16.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_11_1_LTK.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_12_1_sPu.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_13_1_rue.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_16_1_hXV.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_17_1_dmi.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_18_1_oO4.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_19_1_kH7.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_1_1_MA9.root",  "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_21_1_WK2.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_22_1_cR1.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_23_1_cEv.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_24_1_xBV.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_25_1_lkq.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_2_1_hSA.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_3_1_WJu.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_4_1_O5Y.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_5_1_PpA.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_6_1_oUv.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_7_1_XHP.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_8_1_Q89.root", "/store/caf/user/tdahms/Data2013/pPb/PromptSkims/Runs_210658-210737/onia2MuMuPAT_9_1_KVk.root"


options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag ='GR_P_V43D::All' # the same tag that have been used for the skim

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowersPlusTrunc"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("pACentrality")
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    )
)

process.hltMult100DblMu3 = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_PAPixelTrackMultiplicity100_L2DoubleMu3_v*"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)

process.hiz = cms.EDAnalyzer('HiZAnalyzer',
                                srcMuon = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig = cms.InputTag("patMuonsWithoutTrigger"),
                                src = cms.InputTag("onia2MuMuPatTrkTrk"),
                                genParticles = cms.InputTag("generator"),
                                primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
				srcCentrality = cms.InputTag("pACentrality"),

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                               
                                #--
                                maxAbsZ = cms.double(24.0),
                               
                                centralityRanges = cms.vdouble(10,20,30,40,50,100),
       
                                applyCuts = cms.bool(True),           
                     
                                removeSignalEvents = cms.untracked.bool(False),
                                removeTrueMuons = cms.untracked.bool(False),
                                storeSameSign = cms.untracked.bool(True),
                               
                                #-- Gen Details
                                oniaPDG = cms.int32(23),
                                isHI = cms.untracked.bool(True),
                                isMC = cms.untracked.bool(False),
                                isPromptMC = cms.untracked.bool(True),

                                #-- Histogram configuration
                                fillTree = cms.bool(True),
                                minimumFlag = cms.bool(True),
                                fillSingleMuons = cms.bool(True),
                                histFileName = cms.string(options.outputFile),       
                               
                                #--
                                NumberOfTriggers = cms.uint32(7),
                                )

#process.p = cms.Path(process.hltMult100DblMu3 * process.hiz)
process.p = cms.Path(process.hiz)
