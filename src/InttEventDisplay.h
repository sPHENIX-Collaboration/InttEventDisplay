#ifndef INTTEVENTDISPLAY_H__
#define INTTEVENTDISPLAY_H__

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <Acts/Definitions/Algebra.hpp>
#include <TEveTrack.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGLAnnotation.h>
#include <TGLViewer.h>
#include <TEveViewer.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>
#include <TEveProjectionManager.h>
#include <TEveWindow.h>

#include <vector>
#include <string>

/// Class declarations for use in the analysis module
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TTree;
class TH1;
class TH2;
class PHCompositeNode;
class RawClusterContainer;
class RawCluster;
class SvtxTrackMap;
//class JetMap;
class GlobalVertex;
class PHHepMCGenEventMap;
class JetEvalStack;
class JetRecoEval;
class SvtxTrackEval;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class CaloTriggerInfo;
class JetTruthEval;
class SvtxEvalStack;
class JetEvalStack;

class CylinderGeomIntt;
class TrkrClusterContainer;
class ActsGeometry;
class TrkrHitSetContainer;
class PHField;
class PHFieldUtility;
class PHFieldConfig;
class PHG4Reco;
class PHGeomUtility;
//class PHG4GDMLUtility;

class InttG4HitRead;

/////
class TCanvas;
class TMarker;
class TGeoVolume;
class TGeoNode;
/////


/// Definition of this analysis module class
class InttEventDisplay : public SubsysReco
{
 public:
  /// Constructor
  InttEventDisplay(const std::string &name = "InttEventDisplay"
		   /*const std::string &fname = "InttEventDisplay.root"*/);

  // Destructor
  virtual ~InttEventDisplay();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);

  /// Set the minimum jet pT to cut on
  void setMinJetPt(float minjetpt) { m_minjetpt = minjetpt; }

  /// Set the minimum cluster pT to cut on
  void setMinClusPt(float mincluspt) { m_mincluspt = mincluspt; }

  /// Set things to analyze
  void analyzeTracks(bool analyzeTracks) { m_analyzeTracks = analyzeTracks; }
  void analyzeClusters(bool analyzeClusters) { m_analyzeClusters = analyzeClusters; }
  void analyzeJets(bool analyzeJets) { m_analyzeJets = analyzeJets; }
  void analyzeTruth(bool analyzeTruth) { m_analyzeTruth = analyzeTruth; }
  void useTruthClusters(bool useTruthClusters){ m_useTruthClusters = useTruthClusters; }

  //Set the minimum and max number of cluster to cut on
  void setMinNClus(int mincluster)  { m_mincluster = mincluster; }
  void setMaxNClus(int maxcluster)  { m_maxcluster = maxcluster; }

  //save or not save picture
  void savePictures(bool savePictures) { m_savePictures = savePictures; }

  //setting directory data saved
  void setSaveDirectory(std::string saveDirectory) {m_saveDirectory = saveDirectory;}

  void openWindowViewer(bool windowview) {m_windowView = windowview;}

//////

  //input vector to TEveTrack
  void DrawTracks();

  //input vector to TEvePoint
  void DrawHits();
  void DrawVertex();
  void DrawClusters();

  //do other Draw functions
  void Drawall();

  void PrintEidNclusters(TGLViewer*vi);
  void ShowScale(TGLViewer*vi);
  void Loadgeom();

  //get intt geometry from default gGeoManager
  void GetMatrices(std::string name, TGeoVolume *world, TGeoNode *node);
  void Extractinttgeom();

  //display event only one kind viewer
  void Display_3D();
  void Display_rphi();
  void Display_rhoz(TGLViewer::ECameraType Camera = TGLViewer::kCameraOrthoZOY);

  //display event 3D, rphi and rhoz projection at once
  void display();
  
  //make viewer
  void RPhiViewer(TGLViewer*vi);
  void RhoZViewer(TGLViewer*vi,TGLViewer::ECameraType Camera = TGLViewer::kCameraOrthoZOY);
  void ThreeDViewer(TGLViewer*vi);

  //make tab and input viewer
  void Make_2DViewerTab();
  void Make_3DViewerTab();

  //setting camera angle and background color of each viewer
  void Setting_3dGLViewer(TGLViewer*vi);
  void Setting_rphiGLViewer(TGLViewer*vi);
  void Setting_rhozGLViewer(TGLViewer*vi,TGLViewer::ECameraType Camera);

  void SnapShot();

  void Terminate();
  void VectorClear();

  //TEveManager*gEve=0;

 private:
  TCanvas *m_c1;

  std::vector<Acts::Vector3> m_clusters;
  std::vector<Acts::Vector3> m_hits;
  std::vector<Acts::Vector3> m_tracks;
  std::vector<Acts::Vector3> m_vertex;
 
  TEvePointSet * m_ps = nullptr;
  TEvePointSet * m_psv = nullptr;
  TEvePointSet * m_psh = nullptr;
  TEveTrackList * m_list = nullptr;

  TEveViewer *rphiev = 0;
  TEveViewer *rhozev = 0;
  TEveViewer * evev = 0;

  TEveProjectionManager * rphimng;
  TEveProjectionManager * rhozmng;

  TEveGeoTopNode * inttgeom;
  
  //TGLViewer * rphiv = nullptr;
  //TGLViewer * rhozv =  nullptr;

  TGLAnnotation * an = nullptr;

  TEveWindowSlot *slot2D = 0;
  TEveWindowPack *pack2D = 0;

  TEveWindowSlot *slot3D = 0;
  TEveWindowPack *pack3D = 0;
  
  Double_t camcenter[3] = {0., 0., 0.};

  int ievt =0;
  int k=0; //numbering TGeoVolume 

  //int m_mincluster;
  
  
//////

 private:
  /// String to contain the outfile name containing the trees
  std::string m_outfilename;

  /// Fun4All Histogram Manager tool
  Fun4AllHistoManager *m_hm;

  /// A float for cutting on jet pt
  float m_minjetpt;

  /// A float for cutting on cluster pt
  float m_mincluspt;

  /// A boolean for running over tracks
  bool m_analyzeTracks;

  /// A boolean for running over clusters
  bool m_analyzeClusters;

  /// A boolean for running over jets
  bool m_analyzeJets;

  /// A boolean for collecting hepmc information
  bool m_analyzeTruth;

  bool m_useTruthClusters;

  // A int for cutting on number of cluster
  int m_mincluster;
  int m_maxcluster;

  // A boolean for saving picture
  bool m_savePictures;

  // A string for setting directory data saved
  std::string m_saveDirectory;

  // A boolean for viewer window open or not open
  bool m_windowView;

  /// TFile to hold the following TTrees and histograms
  TFile *m_outfile;
  TTree *m_clustertree;
  TTree *m_tracktree;
  TTree *m_hepmctree;
  TTree *m_truthtree;
  TTree *m_recojettree;
  TTree *m_truthjettree;
  TH1 *m_phi_h;
  TH2 *m_eta_phi_h;


  SvtxEvalStack *m_svtxEvalStack = nullptr;
  JetEvalStack *m_jetEvalStack = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterMap = nullptr;
  TrkrHitSetContainer *m_hitMap = nullptr;
  CylinderGeomIntt* m_geom = nullptr;
  PHField *fieldmap = nullptr;
  PHFieldConfig *default_config = nullptr;

  /// Methods for grabbing the data
  void getTracks(PHCompositeNode *topNode);
  void getTruthJets(PHCompositeNode *topNode);
  void getReconstructedJets(PHCompositeNode *topNode);
  void getEMCalClusters(PHCompositeNode *topNode);
  void getHEPMCTruth(PHCompositeNode *topNode);
  void getPHG4Truth(PHCompositeNode *topNode);
  void getNode(PHCompositeNode *topNode);

  std::vector<Acts::Vector3> writeInttClusters(PHCompositeNode */*topNode*/);
  std::vector<Acts::Vector3> writeInttHits(PHCompositeNode */*topNode*/);
  std::vector<Acts::Vector3> writeInttTracks(PHCompositeNode *topNode);
  std::vector<Acts::Vector3> writeInttVertex(PHCompositeNode *topNode);
  //std::vector<double> writeMagnetField(PHCompositeNode *topNode);
  
  void initializeVariables();
  void initializeTrees();

  /**
   * Make variables for the relevant trees
   */

  /// HEPMC Tree variables
  int m_partid1;
  int m_partid2;
  double m_x1;
  double m_x2;
  int m_mpi;
  int m_process_id;
  double m_truthenergy;
  double m_trutheta;
  double m_truthphi;
  double m_truthpx;
  double m_truthpy;
  double m_truthpz;
  double m_truthpt;
  double m_truthp;
  int m_numparticlesinevent;
  int m_truthpid;

  /// Track variables
  double m_tr_px;
  double m_tr_py;
  double m_tr_pz;
  double m_tr_p;
  double m_tr_pt;
  double m_tr_phi;
  double m_tr_eta;
  int m_charge;
  double m_chisq;
  int m_ndf;
  double m_dca;
  double m_tr_x;
  double m_tr_y;
  double m_tr_z;
  int m_truth_is_primary;
  double m_truthtrackpx;
  double m_truthtrackpy;
  double m_truthtrackpz;
  double m_truthtrackp;
  double m_truthtracke;
  double m_truthtrackpt;
  double m_truthtrackphi;
  double m_truthtracketa;
  int m_truthtrackpid;

  /// Reconstructed jet variables
  double m_recojetpt;
  int m_recojetid;
  double m_recojetpx;
  double m_recojetpy;
  double m_recojetpz;
  double m_recojetphi;
  double m_recojetp;
  double m_recojetenergy;
  double m_recojeteta;
  int m_truthjetid;
  double m_truthjetp;
  double m_truthjetphi;
  double m_truthjeteta;
  double m_truthjetpt;
  double m_truthjetenergy;
  double m_truthjetpx;
  double m_truthjetpy;
  double m_truthjetpz;
  double m_dR;

  /// Cluster variables
  double m_clusenergy;
  double m_cluseta;
  double m_clustheta;
  double m_cluspt;
  double m_clusphi;
  double m_cluspx;
  double m_cluspy;
  double m_cluspz;
  double m_E_4x4;

 //Geometry variables
  const static unsigned int m_nInttLayers = 4;
  /// Search window for phi to match intt clusters in cm
  double m_rPhiSearchWin = 0.1;

};

#endif
