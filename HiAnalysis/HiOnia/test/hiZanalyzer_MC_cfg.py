import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIZ")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

#setup any defaults you want
options.outputFile = "Zmumu_Tree_RegitMC.root"
#import os,commands
#def getDirectoryList(path):
#    cmd  = 'ls %s/ ' % (path)
#    file = ["file:%s/%s" % (path,i) for i in commands.getoutput(cmd).split("\n")]
#    return file

#options.inputFiles = getDirectoryList("/afs/cern.ch/work/a/azsigmon/CMSSW_4_4_2_patch5/src/HiSkim/HiOnia2MuMu/test/SkimRegit_tau/root/")

options.inputFiles = "test.root"

options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'STARTHI44_V12::All'

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Drum"), # different from data!
    centralitySrc = cms.InputTag("hiCentrality")
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.hltDoubleMuOpen = cms.EDFilter("HLTHighLevel",
                                       TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                       HLTPaths = cms.vstring("HLT_HIL1DoubleMu0_HighQ_v1"),
                                       eventSetupPathsKey = cms.string(''),
                                       andOr = cms.bool(True),
                                       throw = cms.bool(False)
                                       )

process.hiz = cms.EDAnalyzer('HiZAnalyzer',
                                srcMuon = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig = cms.InputTag("patMuonsWithoutTrigger"),
                                src = cms.InputTag("onia2MuMuPatGlbGlb"),
                                genParticles = cms.InputTag("hiGenParticles"),
                                primaryVertexTag = cms.InputTag("hiSelectedVertex"),

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                                useRapidity = cms.bool(True),
                                
                                #--
                                maxAbsZ = cms.double(24.0),
                                
                                centralityRanges = cms.vdouble(10,20,40,60,100),
	
                                applyCuts = cms.bool(False),
                      
                                removeSignalEvents = cms.untracked.bool(False),
                                removeTrueMuons = cms.untracked.bool(False),
                                storeSameSign = cms.untracked.bool(True),
                                
                                #-- Gen Details
                                oniaPDG = cms.int32(23),
                                isHI = cms.untracked.bool(True),
                                isMC = cms.untracked.bool(True),
                                isPromptMC = cms.untracked.bool(True),

                                #-- Histogram configuration
                                fillTree = cms.bool(True),
                                fillSingleMuons = cms.bool(True),
                                histFileName = cms.string(options.outputFile),
                                
                                #--
                                NumberOfTriggers = cms.uint32(9),
                                )


#process.p = cms.Path(process.hltDoubleMuOpen + process.hiz)
process.p = cms.Path(process.hiz)
