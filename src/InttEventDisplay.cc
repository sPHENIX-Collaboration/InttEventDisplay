#include "InttEventDisplay.h"

/// Cluster/Calorimeter includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
//#include <calobase/RawTower.h>
//#include <calobase/RawTowerContainer.h>
//#include <calobase/RawTowerGeom.h>
//#include <calobase/RawTowerGeomContainer.h>
//#include <calotrigger/CaloTriggerInfo.h>

/// Tracking includes
//#include <globalvertex/GlobalVertex.h>
//#include <globalvertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/InttDefs.h>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <trackreco/PHActsSiliconSeeding.h>

/// Truth evaluation includes
//#include <g4eval/JetEvalStack.h>
#include <g4eval/SvtxEvalStack.h>

/// HEPMC truth includeshitsetkey
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//#include <HepMC/GenEvent.h>
//#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

//#include <phhepmc/PHHepMCGenEvent.h>
//#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

/// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGTab.h>

//TEve
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveViewer.h>
#include <TEveGeoNode.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveViewer.h>
#include <TEveWindow.h>
#include <TEveBrowser.h>
#include <TEveBrowser.h>
#include <TEveGedEditor.h>

//TGeo
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoOverlap.h>

//TGLViewer
#include <TGLViewer.h>
#include <TGLOrthoCamera.h>
#include <TGLAnnotation.h>
#include <TGLCameraOverlay.h>
#include <TGLEmbeddedViewer.h>

/// C++ includes
#include <cassert>
#include <sstream>
#include <string>
#include <typeinfo>

// Geometory includes
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <intt/CylinderGeomIntt.h>

//////////////
#include <TCanvas.h>
//#include <TRandom.h>
//#include <TMarker.h>
//#include <TGraph.h>
//#include <TGraph2D.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <phfield/PHField.h>
#include <phfield/PHField3DCartesian.h>
#include <phfield/PHFieldConfigv1.h>
#include <phfield/PHFieldUtility.h>
#include <phool/PHIODataNode.h>
//#include <g4main/PHG4Reco.h>
//#include <phgeom/PHGeomUtility.h>
//////////////

using namespace std;

/**
 * Constructor of module
 */
InttEventDisplay::InttEventDisplay(const std::string &name, const std::string &filename)
    : SubsysReco(name), m_outfilename(filename), m_hm(nullptr), m_minjetpt(5.0), m_mincluspt(0.25), m_analyzeTracks(true), m_analyzeClusters(true), m_analyzeJets(true), m_analyzeTruth(false)
{
    /// Initialize variables and trees so we don't accidentally access
    /// memory that was never allocated
    initializeVariables();
    initializeTrees();

    //
    m_c1 = NULL;
}

/**
 * Destructor of module
 */
InttEventDisplay::~InttEventDisplay()
{
    delete m_hm;
    delete m_hepmctree;
    delete m_truthjettree;
    delete m_recojettree;
    delete m_tracktree;

    ///
    if (m_c1 != NULL)
        delete m_c1;
}

/**
 * Initialize the module and prepare looping over events
 */
int InttEventDisplay::Init(PHCompositeNode * /*topNode*/)
{
    if (Verbosity() > 5)
    {
        cout << "Beginning Init in InttEventDisplay" << endl;
    }

    m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");

    m_phi_h = new TH1D("phi_h", ";Counts;#phi [rad]", 50, -6, 6);
    m_eta_phi_h = new TH2F("phi_eta_h", ";#eta;#phi [rad]", 10, -1, 1, 50, -6, 6);

    return 0;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int InttEventDisplay::process_event(PHCompositeNode *topNode)
{
    cout << "InttEvt::process evt : " << ievt++ << endl;

    PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

    CylinderGeomIntt *geom = dynamic_cast<CylinderGeomIntt *>(geom_container->GetLayerGeom(3));
    if (geom == NULL)
        cout << "No CylinderGeomIntt" << endl;
    if (Verbosity() > 5)
    {
        cout << "Beginning process_event in AnaTutorial" << endl;
    }

    // Get some Nodes. INTT geometry + hit container
    getNode(topNode);

    m_clusters = writeInttClusters(topNode);
    cout << "number of clusters is " << m_clusters.size() << endl;
    
    m_hits = writeInttHits(topNode);
    cout << "number of hits is " << m_hits.size() << endl;

    m_tracks = writeInttTracks(topNode);
    cout << "number of tracks is " << m_tracks.size() << endl;

    m_vertex = writeInttVertex(topNode);
    cout << "number of vertex is " << m_vertex.size() << endl;

    if(m_clusters.size()>=3){
      InttEventDisplay::display();
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // cout tracks
    /*
    for (auto itr = m_tracks.begin(); itr != m_tracks.end(); ++itr)
    {
      auto track = itr[0];
      cout << "--------------" << endl;
      // cout <<"itr = "  <<*itr << endl;
      cout << "px=" << track[0] << "   py=" << track[1] << "   pz=" << track[2] << endl;
      cout << "--------------" << endl<< endl<< endl<< endl<< endl;
    }
    */

    // cout hit positions
    /*
    for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
      auto cluster = itr[0];
      cout <<"--------------"<<endl;
      cout <<"itr = "  <<*itr << endl;
      cout<<"cluster[0]="<<cluster[0]<<"   cluster[1]="<<cluster[1]<<"   cluster[2]="<<cluster[2]<<endl;
      //cout <<"--------------"<<endl<<endl<<endl<<endl<<endl;
    }
    */

    // cout << "number of points is " << m_clusters.size() << endl<< endl<< endl<< endl;

    return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int InttEventDisplay::End(PHCompositeNode * /*topNode*/)
{
    if (Verbosity() > 1)
    {
        cout << "Ending AnaTutorial analysis package" << endl;
    }

    /*
    /// Change to the outfile
    m_outfile->cd();

    /// If we analyzed the tracks, write the tree out
    if (m_analyzeTracks)
        m_tracktree->Write();

    /// If we analyzed the jets, write them out
    if (m_analyzeJets)
    {
        m_truthjettree->Write();
        m_recojettree->Write();
    }

    /// If we analyzed the truth particles, write them out
    if (m_analyzeTruth)
    {
        m_hepmctree->Write();
        m_truthtree->Write();
    }

    /// If we analyzed the clusters, write them out
    if (m_analyzeClusters)
    {
        m_clustertree->Write();
    }

    /// Write out any other histograms
    m_phi_h->Write();
    m_eta_phi_h->Write();

    /// Write and close the outfile
    m_outfile->Write();
    m_outfile->Close();

    /// Clean up pointers and associated histos/trees in TFile
    delete m_outfile;

    if (Verbosity() > 1)
    {
        cout << "Finished AnaTutorial analysis package" << endl;
    }
    */

    return 0;
}

/**
 * This method gets the tracks as reconstructed from the tracker. It also
 * compares the reconstructed track to its truth track counterpart as determined
 * by the
 */
void InttEventDisplay::getTracks(PHCompositeNode *topNode)
{
    /// SVTX tracks node
    SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    if (!trackmap)
    {
        cout << PHWHERE
             << "SvtxTrackMap node is missing, can't collect tracks"
             << endl;
        return;
    }

    /// EvalStack for truth track matching
    if (!m_svtxEvalStack)
    {
        m_svtxEvalStack = new SvtxEvalStack(topNode);
        m_svtxEvalStack->set_verbosity(Verbosity());
    }

    m_svtxEvalStack->next_event(topNode);

    /// Get the track evaluator
    SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();

    /// Get the range for primary tracks
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    int count1 = 0;
    if (Verbosity() > 1)
    {
        cout << "Get the SVTX tracks" << endl;
    }
    for (SvtxTrackMap::Iter iter = trackmap->begin();
         iter != trackmap->end();
         ++iter)
    {
        SvtxTrack *track = iter->second;

        /// Get the reconstructed track info
        m_tr_px = track->get_px();
        m_tr_py = track->get_py();
        m_tr_pz = track->get_pz();
        m_tr_p = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py + m_tr_pz * m_tr_pz);

        m_tr_pt = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py);

        // Make some cuts on the track to clean up sample

        if (m_tr_pt < 0.5)
            continue;

        m_tr_phi = track->get_phi();
        m_tr_eta = track->get_eta();

        m_charge = track->get_charge();
        m_chisq = track->get_chisq();
        m_ndf = track->get_ndf();
        m_dca = track->get_dca();
        m_tr_x = track->get_x();
        m_tr_y = track->get_y();
        m_tr_z = track->get_z();

        /// Get truth track info that matches this reconstructed track
        PHG4Particle *truthtrack = trackeval->max_truth_particle_by_nclusters(track);
        m_truth_is_primary = truthinfo->is_primary(truthtrack);

        m_truthtrackpx = truthtrack->get_px();
        m_truthtrackpy = truthtrack->get_py();
        m_truthtrackpz = truthtrack->get_pz();
        m_truthtrackp = sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy + m_truthtrackpz * m_truthtrackpz);

        m_truthtracke = truthtrack->get_e();

        m_truthtrackpt = sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy);
        m_truthtrackphi = atan(m_truthtrackpy / m_truthtrackpx);
        m_truthtracketa = atanh(m_truthtrackpz / m_truthtrackp);
        m_truthtrackpid = truthtrack->get_pid();

        count1 = count1 + 1;
        m_tracktree->Fill();
    }
    cout << "get track loop =" << count1 << endl;
}

/**
 * This function puts all of the tree branch assignments in one place so as to not
 * clutter up the AnaTutorial::Init function.
 */
void InttEventDisplay::initializeTrees()
{
    m_recojettree = new TTree("jettree", "A tree with reconstructed jets");
    m_recojettree->Branch("m_recojetpt", &m_recojetpt, "m_recojetpt/D");
    m_recojettree->Branch("m_recojetid", &m_recojetid, "m_recojetid/I");
    m_recojettree->Branch("m_recojetpx", &m_recojetpx, "m_recojetpx/D");
    m_recojettree->Branch("m_recojetpy", &m_recojetpy, "m_recojetpy/D");
    m_recojettree->Branch("m_recojetpz", &m_recojetpz, "m_recojetpz/D");
    m_recojettree->Branch("m_recojetphi", &m_recojetphi, "m_recojetphi/D");
    m_recojettree->Branch("m_recojeteta", &m_recojeteta, "m_recojeteta/D");
    m_recojettree->Branch("m_recojetenergy", &m_recojetenergy, "m_recojetenergy/D");
    m_recojettree->Branch("m_truthjetid", &m_truthjetid, "m_truthjetid/I");
    m_recojettree->Branch("m_truthjetp", &m_truthjetp, "m_truthjetp/D");
    m_recojettree->Branch("m_truthjetphi", &m_truthjetphi, "m_truthjetphi/D");
    m_recojettree->Branch("m_truthjeteta", &m_truthjeteta, "m_truthjeteta/D");
    m_recojettree->Branch("m_truthjetpt", &m_truthjetpt, "m_truthjetpt/D");
    m_recojettree->Branch("m_truthjetenergy", &m_truthjetenergy, "m_truthjetenergy/D");
    m_recojettree->Branch("m_truthjetpx", &m_truthjetpx, "m_truthjetpx/D");
    m_recojettree->Branch("m_truthjetpy", &m_truthjetpy, "m_truthjyetpy/D");
    m_recojettree->Branch("m_truthjetpz", &m_truthjetpz, "m_truthjetpz/D");
    m_recojettree->Branch("m_dR", &m_dR, "m_dR/D");

    m_truthjettree = new TTree("truthjettree", "A tree with truth jets");
    m_truthjettree->Branch("m_truthjetid", &m_truthjetid, "m_truthjetid/I");
    m_truthjettree->Branch("m_truthjetp", &m_truthjetp, "m_truthjetp/D");
    m_truthjettree->Branch("m_truthjetphi", &m_truthjetphi, "m_truthjetphi/D");
    m_truthjettree->Branch("m_truthjeteta", &m_truthjeteta, "m_truthjeteta/D");
    m_truthjettree->Branch("m_truthjetpt", &m_truthjetpt, "m_truthjetpt/D");
    m_truthjettree->Branch("m_truthjetenergy", &m_truthjetenergy, "m_truthjetenergy/D");
    m_truthjettree->Branch("m_truthjetpx", &m_truthjetpx, "m_truthjetpx/D");
    m_truthjettree->Branch("m_truthjetpy", &m_truthjetpy, "m_truthjetpy/D");
    m_truthjettree->Branch("m_truthjetpz", &m_truthjetpz, "m_truthjetpz/D");
    m_truthjettree->Branch("m_dR", &m_dR, "m_dR/D");
    m_truthjettree->Branch("m_recojetpt", &m_recojetpt, "m_recojetpt/D");
    m_truthjettree->Branch("m_recojetid", &m_recojetid, "m_recojetid/I");
    m_truthjettree->Branch("m_recojetpx", &m_recojetpx, "m_recojetpx/D");
    m_truthjettree->Branch("m_recojetpy", &m_recojetpy, "m_recojetpy/D");
    m_truthjettree->Branch("m_recojetpz", &m_recojetpz, "m_recojetpz/D");
    m_truthjettree->Branch("m_recojetphi", &m_recojetphi, "m_recojetphi/D");
    m_truthjettree->Branch("m_recojeteta", &m_recojeteta, "m_recojeteta/D");
    m_truthjettree->Branch("m_recojetenergy", &m_recojetenergy, "m_recojetenergy/D");
    m_tracktree = new TTree("tracktree", "A tree with svtx tracks");
    m_tracktree->Branch("m_tr_px", &m_tr_px, "m_tr_px/D");
    m_tracktree->Branch("m_tr_py", &m_tr_py, "m_tr_py/D");
    m_tracktree->Branch("m_tr_pz", &m_tr_pz, "m_tr_pz/D");
    m_tracktree->Branch("m_tr_p", &m_tr_p, "m_tr_p/D");
    m_tracktree->Branch("m_tr_pt", &m_tr_pt, "m_tr_pt/D");
    m_tracktree->Branch("m_tr_phi", &m_tr_phi, "m_tr_phi/D");
    m_tracktree->Branch("m_tr_eta", &m_tr_eta, "m_tr_eta/D");
    m_tracktree->Branch("m_charge", &m_charge, "m_charge/I");
    m_tracktree->Branch("m_chisq", &m_chisq, "m_chisq/D");
    m_tracktree->Branch("m_ndf", &m_ndf, "m_ndf/I");
    m_tracktree->Branch("m_dca", &m_dca, "m_dca/D");
    m_tracktree->Branch("m_tr_x", &m_tr_x, "m_tr_x/D");
    m_tracktree->Branch("m_tr_y", &m_tr_y, "m_tr_y/D");
    m_tracktree->Branch("m_tr_z", &m_tr_z, "m_tr_z/D");
    m_tracktree->Branch("m_truth_is_primary", &m_truth_is_primary, "m_truth_is_primary/I");
    m_tracktree->Branch("m_truthtrackpx", &m_truthtrackpx, "m_truthtrackpx/D");
    m_tracktree->Branch("m_truthtrackpy", &m_truthtrackpy, "m_truthtrackpy/D");
    m_tracktree->Branch("m_truthtrackpz", &m_truthtrackpz, "m_truthtrackpz/D");
    m_tracktree->Branch("m_truthtrackp", &m_truthtrackp, "m_truthtrackp/D");
    m_tracktree->Branch("m_truthtracke", &m_truthtracke, "m_truthtracke/D");
    m_tracktree->Branch("m_truthtrackpt", &m_truthtrackpt, "m_truthtrackpt/D");
    m_tracktree->Branch("m_truthtrackphi", &m_truthtrackphi, "m_truthtrackphi/D");
    m_tracktree->Branch("m_truthtracketa", &m_truthtracketa, "m_truthtracketa/D");
    m_tracktree->Branch("m_truthtrackpid", &m_truthtrackpid, "m_truthtrackpid/I");

    m_hepmctree = new TTree("hepmctree", "A tree with hepmc truth particles");
    m_hepmctree->Branch("m_partid1", &m_partid1, "m_partid1/I");
    m_hepmctree->Branch("m_partid2", &m_partid2, "m_partid2/I");
    m_hepmctree->Branch("m_x1", &m_x1, "m_x1/D");
    m_hepmctree->Branch("m_x2", &m_x2, "m_x2/D");
    m_hepmctree->Branch("m_mpi", &m_mpi, "m_mpi/I");
    m_hepmctree->Branch("m_process_id", &m_process_id, "m_process_id/I");
    m_hepmctree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
    m_hepmctree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
    m_hepmctree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
    m_hepmctree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
    m_hepmctree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
    m_hepmctree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
    m_hepmctree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
    m_hepmctree->Branch("m_numparticlesinevent", &m_numparticlesinevent, "m_numparticlesinevent/I");
    m_hepmctree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

    m_truthtree = new TTree("truthg4tree", "A tree with truth g4 particles");
    m_truthtree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
    m_truthtree->Branch("m_truthp", &m_truthp, "m_truthp/D");
    m_truthtree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
    m_truthtree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
    m_truthtree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
    m_truthtree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
    m_truthtree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
    m_truthtree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
    m_truthtree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

    m_clustertree = new TTree("clustertree", "A tree with emcal clusters");
    m_clustertree->Branch("m_clusenergy", &m_clusenergy, "m_clusenergy/D");
    m_clustertree->Branch("m_cluseta", &m_cluseta, "m_cluseta/D");
    m_clustertree->Branch("m_clustheta", &m_clustheta, "m_clustheta/D");
    m_clustertree->Branch("m_cluspt", &m_cluspt, "m_cluspt/D");
    m_clustertree->Branch("m_clusphi", &m_clusphi, "m_clusphi/D");
    m_clustertree->Branch("m_cluspx", &m_cluspx, "m_cluspx/D");
    m_clustertree->Branch("m_cluspy", &m_cluspy, "m_cluspy/D");
    m_clustertree->Branch("m_cluspz", &m_cluspz, "m_cluspz/D");
    m_clustertree->Branch("m_E_4x4", &m_E_4x4, "m_E_4x4/D");
}

/**
 * This function initializes all of the member variables in this class so that there
 * are no variables that might not be set before e.g. writing them to the output
 * trees.
 */
void InttEventDisplay::initializeVariables()
{
    m_outfile = new TFile();
    m_phi_h = new TH1F();
    m_eta_phi_h = new TH2F();

    m_partid1 = -99;
    m_partid2 = -99;
    m_x1 = -99;
    m_x2 = -99;
    m_mpi = -99;
    m_process_id = -99;
    m_truthenergy = -99;
    m_trutheta = -99;
    m_truthphi = -99;
    m_truthp = -99;
    m_truthpx = -99;
    m_truthpy = -99;
    m_truthpz = -99;
    m_truthpt = -99;
    m_numparticlesinevent = -99;
    m_truthpid = -99;

    m_tr_px = -99;
    m_tr_py = -99;
    m_tr_pz = -99;
    m_tr_p = -99;
    m_tr_pt = -99;
    m_tr_phi = -99;
    m_tr_eta = -99;
    m_charge = -99;
    m_chisq = -99;
    m_ndf = -99;
    m_dca = -99;
    m_tr_x = -99;
    m_tr_y = -99;
    m_tr_z = -99;
    m_truth_is_primary = -99;
    m_truthtrackpx = -99;
    m_truthtrackpy = -99;
    m_truthtrackpz = -99;
    m_truthtrackp = -99;
    m_truthtracke = -99;
    m_truthtrackpt = -99;
    m_truthtrackphi = -99;
    m_truthtracketa = -99;
    m_truthtrackpid = -99;

    m_recojetpt = -99;
    m_recojetid = -99;
    m_recojetpx = -99;
    m_recojetpy = -99;
    m_recojetpz = -99;
    m_recojetphi = -99;
    m_recojetp = -99;
    m_recojetenergy = -99;
    m_recojeteta = -99;
    m_truthjetid = -99;
    m_truthjetp = -99;
    m_truthjetphi = -99;
    m_truthjeteta = -99;
    m_truthjetpt = -99;
    m_truthjetenergy = -99;
    m_truthjetpx = -99;
    m_truthjetpy = -99;
    m_truthjetpz = -99;
    m_dR = -99;
}

//////////////////////////
//
//    my analysis
//
//////////////////////////

void InttEventDisplay::getNode(PHCompositeNode *topNode)
{
    m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!m_tGeometry)
    {
        std::cout << PHWHERE << "No ActsGeometry on node tree. Bailing."
                  << std::endl;
        return;
    }
    else
        cout << "ActGeometry on node tree is successfully loaded" << endl;

    if (m_useTruthClusters)
    {
        m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
                                                                "TRKR_CLUSTER_TRUTH");
        cout << "TRKR_CLUSTER_TRUTH" << endl;
    }
    else
    {
        m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
                                                                "TRKR_CLUSTER");
        cout << "TRKR_CLUSTER" << endl;
    }

    if (!m_clusterMap)
    {
        cout << PHWHERE
             << "TrkrClusterContainer node is missing."
             << endl;
        // return ABORTEVENT;
        return;
    }
    else
        cout << "TrkrClusterContainr is successfully loaded" << endl;

    return;
}

std::vector<Acts::Vector3> InttEventDisplay ::writeInttHits(PHCompositeNode *topNode)
{
    std::vector<Acts::Vector3> hits;

    PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
    if (!g4hit)
    {
        std::cout << "Could not locate g4 hit node G4HIT_INTT " << std::endl;
        exit(1);
    }
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
    for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
        Acts::Vector3 glob;
        // std::cout << "x: " << hiter->second->get_avg_x() << std::endl;
        glob(0) = hiter->second->get_avg_x();
        // std::cout << "y: " << hiter->second->get_avg_y() << std::endl;
        glob(1) = hiter->second->get_avg_y();
        // std::cout << "z: " << hiter->second->get_avg_z() << std::endl;
        glob(2) = hiter->second->get_avg_z();
        hits.push_back(glob);
    }
    std::cout << "InttG4HitRead::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;

    return hits;
}

std::vector<Acts::Vector3> InttEventDisplay ::writeInttClusters(PHCompositeNode * /*topNode*/)
{
    std::vector<Acts::Vector3> clusters;
    for (unsigned int inttlayer = 0; inttlayer < m_nInttLayers; inttlayer++)
    {
        for (const auto &hitsetkey : m_clusterMap->getHitSetKeys(TrkrDefs::TrkrId::inttId, inttlayer + 3))
        {
            auto range = m_clusterMap->getClusters(hitsetkey);

            for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
            {
                const auto cluskey = clusIter->first;
                const auto cluster = clusIter->second;
                const auto globalPos = m_tGeometry->getGlobalPosition(cluskey, cluster);
                clusters.push_back(globalPos);
            }
        }
    }

    return clusters;
}

std::vector<Acts::Vector3> InttEventDisplay ::writeInttTracks(PHCompositeNode *topNode)
{
    std::vector<Acts::Vector3> tracks;

    // SVTX tracks node
    SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    /// EvalStack for truth track matching
    if (!m_svtxEvalStack)
    {
        m_svtxEvalStack = new SvtxEvalStack(topNode);
        m_svtxEvalStack->set_verbosity(Verbosity());
    }

    m_svtxEvalStack->next_event(topNode);

    Acts::Vector3 globt;

    if (Verbosity() > 1)
    {
        cout << "Get the SVTX tracks" << endl;
    }

    for (SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter)
    {
        SvtxTrack *track = iter->second;

        // Get the reconstructed track info
        m_tr_px = track->get_px();
        globt(0) = m_tr_px;
        // cout << "px =" << m_tr_px << endl;
        m_tr_py = track->get_py();
        // cout << "py =" << m_tr_py << endl;
        globt(1) = m_tr_py;
        m_tr_pz = track->get_pz();
        // cout << "pz =" << m_tr_pz << endl;
        globt(2) = m_tr_pz;
        tracks.push_back(globt);
    }
    return tracks;
}

std::vector<Acts::Vector3> InttEventDisplay ::writeInttVertex(PHCompositeNode *topNode)
{
    std::vector<Acts::Vector3> vertex;

    // SVTX Vertex node
    SvtxVertexMap *vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

    /// EvalStack for truth track matching
    if (!m_svtxEvalStack)
    {
        m_svtxEvalStack = new SvtxEvalStack(topNode);
        m_svtxEvalStack->set_verbosity(Verbosity());
    }

    m_svtxEvalStack->next_event(topNode);

    Acts::Vector3 globv;

    if (Verbosity() > 1)
    {
        cout << "Get the SVTX tracks" << endl;
    }
    for (SvtxVertexMap::Iter iter = vertexmap->begin(); iter != vertexmap->end(); ++iter)
    {
        SvtxVertex *vertexs = iter->second;

        // Get the reconstructed track info
        m_tr_x = vertexs->get_x();
        globv(0) = m_tr_x;
        cout << "x =" << m_tr_x << endl;
        m_tr_y = vertexs->get_y();
        cout << "y =" << m_tr_y << endl;
        globv(1) = m_tr_y;
        m_tr_z = vertexs->get_z();
        cout << "z =" << m_tr_z << endl;
        globv(2) = m_tr_z;
        vertex.push_back(globv);
    }
    return vertex;
}

void InttEventDisplay::Display_3D()
{
    TEveManager::Terminate();
    TEveManager::Create();

    InttEventDisplay::Drawall();
    InttEventDisplay::Loadgeom();

    TGLViewer*v = gEve->GetDefaultGLViewer();
    InttEventDisplay::ThreeDViewer(v);
}

void InttEventDisplay::Display_rphi()
{
    TEveManager::Terminate();
    TEveManager::Create();

    InttEventDisplay::Drawall();
    InttEventDisplay::Loadgeom();

    TGLViewer*vi = gEve->GetDefaultGLViewer();
    InttEventDisplay::RPhiViewer(vi);
}

void InttEventDisplay::Display_rhoz(TGLViewer::ECameraType Camera)
{
    TEveManager::Terminate();
    TEveManager::Create();

    InttEventDisplay::Drawall();
    InttEventDisplay::Loadgeom();

    TGLViewer*vi = gEve->GetDefaultGLViewer();
    InttEventDisplay::RhoZViewer(vi,Camera);
}

void InttEventDisplay::DrawTracks()
{
    int numtrack = m_tracks.size();
    cout << "number of tracks = " << numtrack << endl;
    m_list = new TEveTrackList();

    if (numtrack == 0)
    {
        cout << "No track data" << endl;
        m_list = nullptr;
        return;
    }

    TEveTrackPropagator *prop = m_list->GetPropagator();
    TEveTrack *track[numtrack];
    TEveRecTrackD *rc = new TEveRecTrackD();
    int counter = 0;

    prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 0.));

    m_list->SetElementName(Form("%s, constB", m_list->GetElementName()));
    m_list->SetLineColor(kMagenta);
    auto vertex = m_vertex[0];
    cout << "vertex x =" << vertex[0] << "vertex y =" << vertex[1] << "vertex z =" << vertex[2] << endl;
    for (auto itr = m_tracks.begin(); itr != m_tracks.end(); ++itr)
    {
        auto moment = itr[0];
        rc->fV.Set(vertex[0], vertex[1], vertex[2]);
        rc->fP.Set(moment[0], moment[1], moment[2]);
        rc->fSign = 1;
        track[counter] = new TEveTrack(rc, prop);
        m_list->AddElement(track[counter]);
        track[counter]->SetLineColor(m_list->GetLineColor());
        track[counter]->MakeTrack();
        counter++;
    }

    cout << "counter =" << counter << endl;
}

void InttEventDisplay::DrawClusters()
{
    int npoints = m_clusters.size();
    cout << "Nclusters = " << npoints << endl;
    m_ps = new TEvePointSet(npoints);

    if (npoints == 0)
    {
        cout << "No Cluster data" << endl;
        m_ps = nullptr;
        return;
    }

    m_ps = new TEvePointSet(npoints);
    m_ps->SetOwnIds(kTRUE);

    int counter = 0;
    for (auto itr = m_clusters.begin(); itr != m_clusters.end(); ++itr)
    {
        auto cluster = itr[0];
        m_ps->SetNextPoint(cluster[0], cluster[1], cluster[2]);
        // cout<<"xyz = "<<cluster[0]<<","<< cluster[1]<<","<< cluster[2]<<endl;
        m_ps->SetPointId(new TNamed(Form("Point %d", counter), ""));
        counter++;
    }

    m_ps->SetMarkerColor(2);
    m_ps->SetMarkerSize(1.2);
    m_ps->SetMarkerStyle(4);
}

void InttEventDisplay::DrawVertex()
{

    int nvertex = m_vertex.size();
    cout << "nvertex = " << nvertex << endl;
    // m_psv = new TEvePointSet(nvertex);

    if (nvertex == 0)
    {
        cout << "No vertex data" << endl;
        m_psv = nullptr;
        return;
    }

    m_psv = new TEvePointSet(nvertex);
    m_psv->SetOwnIds(kTRUE);

    int counter = 0;
    for (auto itr = m_vertex.begin(); itr != m_vertex.end(); ++itr)
    {
        auto cluster = itr[0];
        m_psv->SetNextPoint(cluster[0], cluster[1], cluster[2]);
        m_psv->SetPointId(new TNamed(Form("Point %d", counter), ""));
        counter++;
    }

    m_psv->SetMarkerColor(4);
    m_psv->SetMarkerSize(1.2);
    m_psv->SetMarkerStyle(4);
}

void InttEventDisplay::DrawHits()
{
    int nhits = m_hits.size();
    cout << "nhits = " << nhits << endl;
    m_psh = new TEvePointSet(nhits);

    if (nhits == 0)
    {
        cout << "No hit data" << endl;
        m_psh = nullptr;
        return;
    }

    m_psh->SetOwnIds(kTRUE);

    int counter = 0;
    for (auto itr = m_hits.begin(); itr != m_hits.end(); ++itr)
    {
        auto hit = itr[0];
        m_psh->SetNextPoint(hit[0], hit[1], hit[2]);
        m_psh->SetPointId(new TNamed(Form("Point %d", counter), ""));
        counter++;
    }

    m_psh->SetMarkerColor(5);
    m_psh->SetMarkerSize(1.2);
    m_psh->SetMarkerStyle(4);
}

void InttEventDisplay::GetMatrices(string name, TGeoVolume *world, TGeoNode *node)
{
    if (!node)
    {
        return;
    }

    TGeoVolume *volume = node->GetVolume();

    if (volume->GetName() == name)
    {

        TGeoMatrix *matrix = node->GetMatrix();

        if (volume && matrix)
        {
            world->AddNode(volume, k, matrix);
            k++;
        }
    }
    for (Int_t i = 0; i < node->GetNdaughters(); ++i)
    {
        TGeoNode *daughterNode = node->GetDaughter(i);
        GetMatrices(name, world, daughterNode);
    }
}

void InttEventDisplay::Extractinttgeom()
{
    // get GeoMatrix from nodes
    TGeoNode *topNode = gGeoManager->GetTopNode();
    TGeoVolume *intt = gGeoManager->MakeBox("BOX", gGeoManager->GetTopVolume()->GetMedium(), 1, 1, 1);
    for (int inttlayer = 0; inttlayer < 4; inttlayer++)
    {
        for (int itype = 0; itype < 2; itype++)
        {
            string name = "ladder_" + to_string(inttlayer) + "_" + to_string(itype);
            GetMatrices(name, intt, topNode);
        }
    }

    gGeoManager->SetTopVolume(intt);
    intt->SetVisibility(kFALSE);
    gGeoManager->CloseGeometry();
}

void InttEventDisplay::PrintEidNclusters(TGLViewer*vi)
{
    string legend = "eid :" + to_string(ievt - 1) + ", Nclusters: " + to_string(m_clusters.size());
    an = new TGLAnnotation(vi, legend.c_str(), 0.05, 0.95);
    an->SetTextSize(0.05);
    an->SetTextColor(kBlack);
}

void InttEventDisplay::ShowScale(TGLViewer*vi)
{
    // scale
    TGLCameraOverlay *glovlay = vi->GetCameraOverlay();
    glovlay->SetShowOrthographic(1);
}

void InttEventDisplay::Drawall()
{

    InttEventDisplay::DrawHits();
    if (m_psh != nullptr)
    {
        gEve->AddElement(m_psh);
    }

    InttEventDisplay::DrawClusters();
    if (m_ps != nullptr)
    {
        gEve->AddElement(m_ps);
    }

    InttEventDisplay::DrawTracks();
    if (m_list != nullptr)
    {
        gEve->AddElement(m_list);
    }

    InttEventDisplay::DrawVertex();
    if (m_psv != nullptr)
    {
        gEve->AddElement(m_psv);
    }
}

void InttEventDisplay::Setting_3dGLViewer(TGLViewer*vi){
  //show xyz axes
  vi->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

  //set background color
  vi->SetClearColor(kWhite);

  //set initial camera angle
  vi->CurrentCamera().RotateRad(0, -3.14 / 2);
  vi->SetPerspectiveCamera(TGLViewer::kCameraPerspXOZ, 2, 30, &camcenter[0], 0, 0);
}

void InttEventDisplay::Setting_rphiGLViewer(TGLViewer*vi){
  //set background color
  vi->SetClearColor(kWhite);

  //set initial camera
  vi->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 9, 100, &camcenter[0], 0., 0.);
  vi->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
}

void InttEventDisplay::Setting_rhozGLViewer(TGLViewer*vi,TGLViewer::ECameraType Camera){
  //set background color
  vi->SetClearColor(kWhite);
  
  //set initial camera
  vi->SetOrthoCamera(Camera, 4, 100, &camcenter[0], 0., 0.);
  vi->SetCurrentCamera(Camera); 
}

void InttEventDisplay::Loadgeom(){
  // geometry load
  gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry_rphi_black2mm.root");
  TEveGeoTopNode *geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  geom->CanEditMainTransparency();
  geom->SetMainTransparency(70);
  gEve->AddGlobalElement(geom);
}

void InttEventDisplay::RPhiViewer(TGLViewer*vi){

  TEveScene *s = gEve->SpawnNewScene("Rphi Projected Event");
  gEve->GetDefaultViewer()->AddScene(s);
  
  //vi = gEve->GetDefaultGLViewer();
  InttEventDisplay::Setting_rphiGLViewer(vi);
  
  InttEventDisplay::PrintEidNclusters(vi);
  InttEventDisplay::ShowScale(vi);
  
  // projections
  TEveProjectionManager *mng =
    new TEveProjectionManager(TEveProjection::kPT_RPhi);
  s->AddElement(mng);
  gEve->AddToListTree(mng, kTRUE);
  
  gEve->Redraw3D(kFALSE, kTRUE);
  //v->RequestDraw();
}

void InttEventDisplay::RhoZViewer(TGLViewer*vi,TGLViewer::ECameraType Camera){
  TEveScene *s = gEve->SpawnNewScene("Rhoz Projected Event");
  gEve->GetDefaultViewer()->AddScene(s);
  //v = gEve->GetDefaultGLViewer();
  InttEventDisplay::Setting_rhozGLViewer(vi,Camera);
  InttEventDisplay::PrintEidNclusters(vi);
  InttEventDisplay::ShowScale(vi);
  
  // projections
  TEveProjectionManager *mng =
    new TEveProjectionManager(TEveProjection::kPT_RhoZ);
  s->AddElement(mng);
  gEve->AddToListTree(mng, kTRUE);
  
  gEve->Redraw3D(kFALSE, kFALSE);
  //v->RequestDraw();
}

 void InttEventDisplay::ThreeDViewer(TGLViewer*vi){
   InttEventDisplay::Setting_3dGLViewer(vi);
   InttEventDisplay::PrintEidNclusters(vi);
   InttEventDisplay::ShowScale(vi);
   gEve->Redraw3D(kFALSE, kFALSE);
 }

void InttEventDisplay::display(){
  TEveManager::Terminate();
  TEveManager::Create();

  InttEventDisplay::Drawall();
  InttEventDisplay::Loadgeom();

  InttEventDisplay::Make_2DViewerTab();
  InttEventDisplay::Make_3DViewerTab();

  gEve->GetViewers()->SwitchColorSet();
  gEve->GetBrowser()->GetTabRight()->SetTab(1);
}

void InttEventDisplay::Make_2DViewerTab(){
  TEveWindowSlot *slot = 0;
  TEveViewer *rphiv = 0;
  TEveViewer *rhozv = 0;
  TEveWindowPack *pack1;

  slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  pack1 = slot->MakePack();
  pack1->SetShowTitleBar(kFALSE);
  pack1->SetElementName("2D Viewer");
  pack1->SetHorizontal();
  
  // rphiviewer
  slot = pack1->NewSlot();
  rphiv = new TEveViewer("rphi viewer");
  rphiv->SpawnGLEmbeddedViewer(gEve->GetEditor());
  InttEventDisplay::RPhiViewer(rphiv->GetGLViewer());
  slot->ReplaceWindow(rphiv);
  rphiv->SetElementName("rphi viewer");
  gEve->GetViewers()->AddElement(rphiv);
  rphiv->AddScene(gEve->GetGlobalScene());
  rphiv->AddScene(gEve->GetEventScene());

  //rhozviewer
  slot = pack1->NewSlot();
  rhozv = new TEveViewer("rhoz viewer");
  rhozv->SpawnGLViewer(gEve->GetEditor());
  InttEventDisplay::RhoZViewer(rhozv->GetGLViewer());
  slot->ReplaceWindow(rhozv);
  rhozv->SetElementName("rhoz viewer");
  gEve->GetViewers()->AddElement(rhozv);
  rhozv->AddScene(gEve->GetGlobalScene());
  rhozv->AddScene(gEve->GetEventScene());
}

 void InttEventDisplay::Make_3DViewerTab(){
   TEveWindowSlot *slot = 0;
   TEveViewer * evev = 0;
   TEveWindowPack *pack1;
   
   slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
   pack1 = slot->MakePack();
   slot=pack1->NewSlot();
   pack1->SetShowTitleBar(kFALSE);
   pack1->SetElementName("3D Viewer");

   evev = new TEveViewer("3D Viewer");
   evev->SpawnGLEmbeddedViewer(gEve->GetEditor());
   InttEventDisplay::ThreeDViewer(evev->GetGLViewer());
   slot->ReplaceWindow(evev);
   evev->SetElementName("3D viewer");
   gEve->GetViewers()->AddElement(evev);

   evev->AddScene(gEve->GetGlobalScene());
   evev->AddScene(gEve->GetEventScene());
 }
