#!/usr/bin/env python
# @(#)root/tmva $Id: TMVAClassification.py 38475 2011-03-17 10:46:00Z evt $
# ------------------------------------------------------------------------------ #
# Project      : TMVA - a Root-integrated toolkit for multivariate data analysis #
# Package      : TMVA                                                            #
# Python script: TMVAClassification.py                                           #
#                                                                                #
# This python script provides examples for the training and testing of all the   #
# TMVA classifiers through PyROOT.                                               #
#                                                                                #
# The Application works similarly, please see:                                   #
#    TMVA/macros/TMVAClassificationApplication.C                                 #
# For regression, see:                                                           #
#    TMVA/macros/TMVARegression.C                                                #
#    TMVA/macros/TMVARegressionpplication.C                                      #
# and translate to python as done here.                                          #
#                                                                                #
# The methods to be used can be switched on and off via the prompt command, for  #
# example:                                                                       #
#                                                                                #
#    python TMVAClassification.py --methods Fisher,Likelihood                    #
#                                                                                #
# The output file "TMVA.root" can be analysed with the use of dedicated          #
# macros (simply say: root -l <../macros/macro.C>), which can be conveniently    #
# invoked through a GUI that will appear at the end of the run of this macro.    #
#                                                                                #
# for help type "python TMVAClassification.py --help"                            #
# ------------------------------------------------------------------------------ #

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import ROOT
# --------------------------------------------

# Default settings for command line arguments
DEFAULT_OUTFNAME         = "20141208.10.14_ForwardPileupJets.Pythia8_mc12_8TeV_JZ2W_NTUP_COMMON_TowersV4.forwardnew.TMVAtrained.EtaGT2p4_pt20to30.root"
DEFAULT_INFNAME          = "/atlas/output/pnef/20141208.10.14_ForwardPileupJets.Pythia8_mc12_8TeV_JZ2W_NTUP_COMMON_TowersV4.forwardnew.root"
DEFAULT_TREESIG          = "PUJetsTreeFiller/EventTree"
DEFAULT_TREEBKG          = "PUJetsTreeFiller/EventTree"
DEFAULT_METHODS          = "kNN100trim,likelihood"
DEFAULT_NEVENTS_TEST_S   =0
DEFAULT_NEVENTS_TEST_B   =0
DEFAULT_NEVENTS_TRAIN_S  =0
DEFAULT_NEVENTS_TRAIN_B  =0

VARIABLES             = ["CaloTowersWidth", "BackToBackJVT"]
SELECTION             = "fabs(Jeta)>2.4 && Jpt>20 && Jpt<30 && CaloTowersWidth>=0 && CaloTowersWidth<2 "
SPECTATORS            = ["NPV", "Mu" ,"Jpt", "Jeta", "isPileup", "isHardScatter","EventWeight","isQCDPileup3","NClus"]


# Print usage help ------------------------------------------------------------------
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: all methods)"
    print "  -i | --inputfile  : name of input ROOT file (default: '%s')" % DEFAULT_INFNAME
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % DEFAULT_OUTFNAME
    print "  -t | --inputtrees : input ROOT Trees for signal and background (default: '%s %s')" \
          % (DEFAULT_TREESIG, DEFAULT_TREEBKG)
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

# Main routine
def main():

    try:
        # retrive command line options
        shortopts  = "m:i:t:o:vh?"
        longopts   = ["methods=", "inputfile=", "inputtrees=", "outputfile=", "verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    infname     = DEFAULT_INFNAME
    treeNameSig = DEFAULT_TREESIG
    treeNameBkg = DEFAULT_TREEBKG
    outfname    = DEFAULT_OUTFNAME
    methods     = DEFAULT_METHODS
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods = a
        elif o in ("-i", "--inputfile"):
            infname = a
        elif o in ("-o", "--outputfile"):
            outfname = a
        elif o in ("-t", "--inputtrees"):
            a.strip()
            trees = a.rsplit( ' ' )
            trees.sort()
            trees.reverse()
            if len(trees)-trees.count('') != 2:
                print "ERROR: need to give two trees (each one for signal and background)"
                print trees
                sys.exit(1)
            treeNameSig = trees[0]
            treeNameBkg = trees[1]
        elif o in ("-v", "--verbose"):
            verbose = True

    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()

    # Import ROOT classes
    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    
    # check ROOT version, give alarm if 5.18 
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)."
        sys.exit(1)
    
    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    gROOT.SetMacroPath( "./" )
    gROOT.Macro       ( "./TMVAlogon.C" )    
    gROOT.LoadMacro   ( "./TMVAGui.C" )
    
    # Import TMVA classes from ROOT
    from ROOT import TMVA

    # Output file
    outputFile = TFile( outfname, 'RECREATE' )
    
    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
#    factory = TMVA.Factory( "TMVAClassification", outputFile, 
#                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )
    jobname =DEFAULT_OUTFNAME
    factory = TMVA.Factory( jobname.replace(".root", ""), outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" ) # pascal

    # Set verbosity
    factory.SetVerbose( verbose )


    
    # Define the input variables that shall be used for the classifier training
    # note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    for var in VARIABLES:
        factory.AddVariable( var, 'F')

    # You can add so-called "Spectator variables", which are not used in the MVA training, 
    for spect in SPECTATORS:
        factory.AddSpectator( spect, spect)

    # Apply additional cuts on the signal and background sample. 
    mycutBkg = TCut( SELECTION+"&&isPileup") 
    mycutSig = TCut( SELECTION+"&&isHardScatter" ) 

    # open file
    input = TFile.Open( infname )

    # Get the signal and background trees for training
    signal      = input.Get( treeNameSig )
    background  = input.Get( treeNameBkg )

    # Global event weights (see below for setting event-wise weights)
    signalWeight     = 1.0
    backgroundWeight = 1.0

    # ====== register trees ====================================================
    factory.AddSignalTree    ( signal,     signalWeight     )
    factory.AddBackgroundTree( background, backgroundWeight )

    # To give different trees for training and testing, do as follows:
    #    factory.AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" )
    #    factory.AddSignalTree( signalTestTree,     signalTestWeight,  "Test" )
    
    # Set individual event weights (the variables must exist in the original TTree)
    #    for signal    : factory.SetSignalWeightExpression    ("weight1*weight2");
    #    for background: factory.SetBackgroundWeightExpression("weight1*weight2");

    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples
    TrainingAndTestTreeStr= "nTrain_Signal="+str(DEFAULT_NEVENTS_TRAIN_S)+\
                            ":nTrain_Background="+str(DEFAULT_NEVENTS_TRAIN_B)+\
                            ":nTest_Signal="+str(DEFAULT_NEVENTS_TEST_S)+\
                            ":nTest_Background="+str(DEFAULT_NEVENTS_TEST_B)+\
                            ":SplitMode=Random:NormMode=EqualNumEvents:!V"
    factory.PrepareTrainingAndTestTree( mycutSig, mycutBkg, TrainingAndTestTreeStr)

    # --------------------------------------------------------------------------------------------------

    # ---- Book MVA methods
    #
    if "kNN100trim" in mlist:
        factory.BookMethod( TMVA.Types.kKNN, "KNN100trim",
                "!V:H:nkNN=100:ScaleFrac=0.8:UseKernel=F:UseWeight=F:Trim=True:BalanceDepth=6" )
    
    if "likelihood" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "Likelihood","H:!V:" );
    

    if "BDT" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDT",
                           "!H:NTrees=850:nEventsMin=150:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" )

    
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs. 
    # Train MVAs
    factory.TrainAllMethods()
    
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % outfname
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros    
    gROOT.ProcessLine( "TMVAGui(\"%s\")" % outfname )
    
    # keep the ROOT thread running
    gApplication.Run() 

# ----------------------------------------------------------


if __name__ == "__main__":
    main()
