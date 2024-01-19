#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
//#include <G4_Bbc.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_QA.C>

#include <Trkr_Diagnostics.C>
#include <G4_User.C>
#include <QA.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <../src/InttEventDisplay.h>

R__LOAD_LIBRARY(libintteventdisplay.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)

InttEventDisplay*inttEventDisplay;
Fun4AllServer * se;// = Fun4AllServer::instance();

void Fun4All_InttEventDisplay(string inputfilename=/*"/sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/AnaTutorial/macro/dst_intt_cosmic_run25184.root",*/
	      "/sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/AnaTutorial/macro/dst_intt_beam_run20868.root",
	      int minclus =10,int maxclus=50, bool savePictures=false,bool windowview = true,
string saveDirectory="/sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/InttEventDisplay/picture/run20868/"){
  
  const char*inputfile = inputfilename.c_str();
  se = Fun4AllServer::instance();
  se->Verbosity(0);
  const string &outputFile = "G4sPHENIX";
  string outputroot = outputFile;
  const int nEvents = 0;

  Fun4AllInputManager*in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(inputfile);
  se->registerInputManager(in);

  
  //======================
  // What to run
  //======================
  //Enable::MICROMEGAS = true;
  recoConsts *rc = recoConsts::instance();
  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);
  
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", CDB::global_tag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);

  TrackingInit();
  
  Tracking_Reco();
  
  inttEventDisplay = new InttEventDisplay("inttEventDisplay", outputroot + "_inttEventDisplay.root");
  
  //bool windowview = false;
  inttEventDisplay->setMinJetPt(10.);
  inttEventDisplay->Verbosity(0);
  inttEventDisplay->analyzeTracks(true);
  inttEventDisplay->analyzeClusters(true);
  inttEventDisplay->analyzeJets(false);
  inttEventDisplay->analyzeTruth(false);
  inttEventDisplay->useTruthClusters(false);
  inttEventDisplay->setMinNClus(minclus);
  inttEventDisplay->setMaxNClus(maxclus);
  inttEventDisplay->savePictures(savePictures);
  inttEventDisplay->setSaveDirectory(saveDirectory);
  inttEventDisplay->openWindowViewer(windowview);

  se->registerSubsystem(inttEventDisplay);
  se->run(nEvents);
  
  if(windowview==false){
    se->End();
    std::cout << "All done" << std::endl;
    delete se;
    //--  if (Enable::PRODUCTION)
    //--  {
    //--    Production_MoveOutput();
    //--  }
    
    gSystem->Exit(0);
  }
  
  return 0;
}
