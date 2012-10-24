import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIZ")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = "Zmumu_Tree_RegitData_147mub.root"
import os,commands
def getDirectoryList(path):
    cmd  = 'ls %s/ ' % (path)
    file = ["file:%s/%s" % (path,i) for i in commands.getoutput(cmd).split("\n")]
    return file

options.inputFiles = getDirectoryList("/afs/cern.ch/work/m/mironov/public/regitSkim147mub/skim")

options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag ='GR_R_44_V10::All' # the same tag that have been used for the skim
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"), #HFhits for prompt reco
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    )
)

process.hltDoubleMuOpen = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_HIL1DoubleMuOpen"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)

process.hiz = cms.EDAnalyzer('HiZAnalyzer',
                                srcMuon = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig = cms.InputTag("patMuonsWithoutTrigger"),
                                src = cms.InputTag("onia2MuMuPatGlbGlb"),
                                #src = cms.InputTag("onia2MuMuPatStaSta"),
                                genParticles = cms.InputTag("hiGenParticles"),
                                primaryVertexTag = cms.InputTag("hiSelectedVertex"),

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                                useRapidity = cms.bool(True),
                               
                                #--
                                maxAbsZ = cms.double(24.0),
                               
                                centralityRanges = cms.vdouble(10,20,30,40,50,100),
       
                                applyCuts = cms.bool(False),           
                     
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
                                NumberOfTriggers = cms.uint32(9),
                                )


#process.p = cms.Path(process.hltDoubleMuOpen + process.hiz)
process.p = cms.Path(process.hiz)
