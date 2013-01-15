import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIZ")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = "Zmumu_Tree_pPbData_test.root"
#import os,commands
#def getDirectoryList(path):
#    cmd  = 'ls %s/ ' % (path)
#    file = ["file:%s/%s" % (path,i) for i in commands.getoutput(cmd).split("\n")]
#    return file

#options.inputFiles = getDirectoryList("/afs/cern.ch/work/m/mironov/public/regitSkim147mub/skim")

options.inputFiles = "file:/afs/cern.ch/user/t/tdahms/public/ForAnna/onia2MuMuPAT.root"

options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

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

process.hiz = cms.EDAnalyzer('HiZAnalyzer',
                                srcMuon = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig = cms.InputTag("patMuonsWithoutTrigger"),
                                src = cms.InputTag("onia2MuMuPatTrkTrk"),
                                genParticles = cms.InputTag("generator"),
                                primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),

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

process.p = cms.Path(process.hiz)
