#include "InttEventDisplay.h"

/// Cluster/Calorimeter includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>

/// Jet includes
//#include <g4jets/Jet.h>
//#include <g4jets/JetMap.h>

/// Tracking includes
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
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
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <trackreco/PHActsSiliconSeeding.h>

/// Truth evaluation includes
#include <g4eval/JetEvalStack.h>
#include <g4eval/SvtxEvalStack.h>

/// HEPMC truth includeshitsetkey
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

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

/// C++ includes
#include <cassert>
#include <sstream>
#include <string>
#include <typeinfo>

// Geometory includes
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <intt/CylinderGeomIntt.h>

//////////////
#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveViewer.h>
#include <TGLViewer.h>
#include <TGeoManager.h>
#include <TEveGeoNode.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TMarker.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <trackbase/InttDefs.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TGLOrthoCamera.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <phfield/PHField.h>
#include <phfield/PHField3DCartesian.h>
#include <phfield/PHFieldConfigv1.h>
#include <phfield/PHFieldUtility.h>
#include <phool/PHIODataNode.h>
#include <g4main/PHG4Reco.h>
#include <TGeoOverlap.h>
#include <phgeom/PHGeomUtility.h>
#include <TGLAnnotation.h>
#include <TEveViewer.h>
#include <TEveWindow.h>
#include <TEveBrowser.h>
#include <TGLCameraOverlay.h>
//#include <g4gdml/PHG4GDMLUtility.hh>
//////////////

using namespace std;

/**
 * This class demonstrates the construction and use of an analysis module
 * within the sPHENIX Fun4All framework. It is intended to show how to
 * obtain physics objects from the analysis tree, and then write them out
 * to a ROOT tree and file for further offline analysis.
 */

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

  //const string inttgdmlgeom = PHGeomUtility::GenerateGeometryFileName("gdml");
  //PHG4Reco::Dump_GDML(inttgdmlgeom);

  return 0;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int InttEventDisplay::process_event(PHCompositeNode *topNode)
{

  cout << "process_event : start" << endl<< endl<< endl<< endl<< endl;
  cout << "process_event : start" << endl<< endl<< endl<< endl<< endl;
  cout << "process_event : start" << endl<< endl<< endl<< endl<< endl;
  cout << "process_event : start" << endl<< endl<< endl<< endl<< endl;
  cout << "process_event : start" << endl<< endl<< endl<< endl<< endl;

  //static int ievt=0;
  cout<<"InttEvt::process evt : "<<ievt++<<endl;

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  CylinderGeomIntt *geom = dynamic_cast<CylinderGeomIntt *>(geom_container->GetLayerGeom(3));
  if (geom == NULL)
    cout << "No CylinderGeomIntt" << endl<< endl<< endl<< endl;
  if (Verbosity() > 5)
  {
    cout << "Beginning process_event in AnaTutorial" << endl;
  }
  /// Get the truth information
  if (m_analyzeTruth && 0)
  {
    getHEPMCTruth(topNode);
    getPHG4Truth(topNode);
  }

  /// Get the tracks
  if (m_analyzeTracks && 0)
  {
    getTracks(topNode);
  }

  /// Get the truth and reconstructed jets  if (m_analyzeJets)
  /*
  i(0){
    getTruthJets(topNode);
    cout<<"gettruthjets"<<endl;
    getReconstructedJets(topNode);
    cout<<"gererconstrutedjets"<<endl;
  }
  */

  // Get some Nodes. INTT geometry + hit container
  getNode(topNode);
  
  m_hits = writeInttHits(topNode);
  cout << "number of hits is "<< m_hits.size()<<endl;
  
  m_clusters = writeInttClusters(topNode);
  cout << "number of clusters is " << m_clusters.size() << endl<< endl<< endl<< endl;

  m_tracks = writeInttTracks(topNode);
  cout << "number of tracks is " << m_tracks.size() << endl<< endl<< endl<< endl;
  m_vertex = writeInttVertex(topNode);
  cout << "number of vertex is "<< m_vertex.size()<<endl;
  
  //cout <<"Event number ="<< se->EventCounter() <<endl;
  //m_bfield=writeMagnetField(topNode);
  //cout<<"bfield (x,y,z)= ("<<m_bfield[0]<<","<<m_bfield[1]<<","<<m_bfield[2]<<")"<<endl;

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

  //cout hit positions
  /*
  for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
    auto cluster = itr[0];
    cout <<"--------------"<<endl;
    cout <<"itr = "  <<*itr << endl;
    cout<<"cluster[0]="<<cluster[0]<<"   cluster[1]="<<cluster[1]<<"   cluster[2]="<<cluster[2]<<endl;
    //cout <<"--------------"<<endl<<endl<<endl<<endl<<endl;
  }
  */

  //cout << "number of points is " << m_clusters.size() << endl<< endl<< endl<< endl;

  /*
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
    std::cout << "x: " << hiter->second->get_avg_x() << std::endl;
    glob(0)=hiter->second->get_avg_x();
    std::cout << "y: " << hiter->second->get_avg_y() << std::endl;
    glob(1)=hiter->second->get_avg_y();
    std::cout << "z: " << hiter->second->get_avg_z() << std::endl;
    glob(2)=hiter->second->get_avg_z();
    m_hits.push_back(glob);
  }
  std::cout << "InttG4HitRead::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  */

  //const string inttgdmlgeom = PHGeomUtility::GenerateGeometryFileName("gdml");
  //Dump_GDML(inttgdmlgeom);

  //gGeoManager->Export("Full_geom.gdml");

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

  return 0;
}

/**
 * This method gets all of the HEPMC truth particles from the node tree
 * and stores them in a ROOT TTree. The HEPMC truth particles are what,
 * for example, directly comes out of PYTHIA and thus gives you all of
 * the associated parton information
 */
void InttEventDisplay::getHEPMCTruth(PHCompositeNode *topNode)
{
  /// Get the node from the node tree
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << endl;
    return;
  }

  /// Could have some print statements for debugging with verbosity
  if (Verbosity() > 1)
  {
    cout << "Getting HEPMC truth particles " << endl;
  }

  /// You can iterate over the number of events in a hepmc event
  /// for pile up events where you have multiple hard scatterings per bunch crossing
  for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
       eventIter != hepmceventmap->end();
       ++eventIter)
  {
    /// Get the event
    PHHepMCGenEvent *hepmcevent = eventIter->second;

    if (hepmcevent)
    {
      /// Get the event characteristics, inherited from HepMC classes
      HepMC::GenEvent *truthevent = hepmcevent->getEvent();
      if (!truthevent)
      {
        cout << PHWHERE
             << "no evt pointer under phhepmvgeneventmap found "
             << endl;
        return;
      }

      /// Get the parton info
      HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();

      /// Get the parton info as determined from HEPMC
      m_partid1 = pdfinfo->id1();
      m_partid2 = pdfinfo->id2();
      m_x1 = pdfinfo->x1();
      m_x2 = pdfinfo->x2();

      /// Are there multiple partonic intercations in a p+p event
      m_mpi = truthevent->mpi();

      /// Get the PYTHIA signal process id identifying the 2-to-2 hard process
      m_process_id = truthevent->signal_process_id();

      if (Verbosity() > 2)
      {
        cout << " Iterating over an event" << endl;
      }
      /// Loop over all the truth particles and get their information
      for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
           iter != truthevent->particles_end();
           ++iter)
      {
        /// Get each pythia particle characteristics
        m_truthenergy = (*iter)->momentum().e();
        m_truthpid = (*iter)->pdg_id();

        m_trutheta = (*iter)->momentum().pseudoRapidity();
        m_truthphi = (*iter)->momentum().phi();
        m_truthpx = (*iter)->momentum().px();
        m_truthpy = (*iter)->momentum().py();
        m_truthpz = (*iter)->momentum().pz();
        m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

        /// Fill the truth tree
        m_hepmctree->Fill();
        m_numparticlesinevent++;
      }
    }
  }
}

/**
 * This function collects the truth PHG4 stable particles that
 * are produced from PYTHIA, or some other event generator. These
 * are the stable particles, e.g. there are not any (for example)
 * partons here.
 */
void InttEventDisplay::getPHG4Truth(PHCompositeNode *topNode)
{
  /// G4 truth particle node
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return;
  }

  /// Get the primary particle range
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

    m_truthphi = atan(m_truthpy / m_truthpx);

    m_trutheta = atanh(m_truthpz / m_truthenergy);
    /// Check for nans
    if (m_trutheta != m_trutheta)
      m_trutheta = -99;
    m_truthpid = truth->get_pid();

    /// Fill the g4 truth tree
    m_truthtree->Fill();
  }
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
 * Method that gets the truth jets and stores them in a tree
 */
/*
void InttEventDisplay::getTruthJets(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    cout << "get the truth jets" << endl;
  }

  /// Get the truth jet node
  JetMap *truth_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Truth_r04");

  /// Get reco jets associated to truth jets to study e.g. jet efficiencies
  JetMap *reco_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Tower_r04");
  if (!m_jetEvalStack)
  {
    m_jetEvalStack = new JetEvalStack(topNode, "AntiKt_Tower_r04",
                                      "AntiKt_Truth_r04");
  }
  m_jetEvalStack->next_event(topNode);
  JetTruthEval *trutheval = m_jetEvalStack->get_truth_eval();

  if (!truth_jets)
  {
    cout << PHWHERE
         << "Truth jet node is missing, can't collect truth jets"
         << endl;
    return;
  }

  /// Iterate over the truth jets
  for (JetMap::Iter iter = truth_jets->begin();
       iter != truth_jets->end();
       ++iter)
  {
    Jet *truthJet = iter->second;

    m_truthjetpt = truthJet->get_pt();

    std::set<PHG4Particle *> truthjetcomp =
        trutheval->all_truth_particles(truthJet);
    int ntruthconstituents = 0;
    // loop over the constituents of the truth jet
    for (std::set<PHG4Particle *>::iterator iter2 = truthjetcomp.begin();
         iter2 != truthjetcomp.end();
         ++iter2)
    {
      // get the particle of the truthjet
      PHG4Particle *truthpart = *iter2;
      if (!truthpart)
      {
        cout << "no truth particles in the jet??" << endl;
        break;
      }

      ntruthconstituents++;
    }

    if (ntruthconstituents < 3)
      continue;
    /// Only collect truthjets above the _minjetpt cut
    if (m_truthjetpt < m_minjetpt)
      continue;

    m_truthjeteta = truthJet->get_eta();
    m_truthjetpx = truthJet->get_px();
    m_truthjetpy = truthJet->get_py();
    m_truthjetpz = truthJet->get_pz();
    m_truthjetphi = truthJet->get_phi();
    m_truthjetp = truthJet->get_p();
    m_truthjetenergy = truthJet->get_e();

    m_recojetpt = 0;
    m_recojetid = 0;
    m_recojetpx = 0;
    m_recojetpy = 0;
    m_recojetpz = 0;
    m_recojetphi = 0;
    m_recojetp = 0;
    m_recojetenergy = 0;
    m_dR = -99;
    float closestjet = 9999;
    /// Iterate over the reconstructed jets
    for (JetMap::Iter recoIter = reco_jets->begin();
         recoIter != reco_jets->end();
         ++recoIter)
    {
      const Jet *recoJet = recoIter->second;
      m_recojetpt = recoJet->get_pt();
      if (m_recojetpt < m_minjetpt - m_minjetpt * 0.4)
        continue;

      m_recojeteta = recoJet->get_eta();
      m_recojetphi = recoJet->get_phi();

      if (Verbosity() > 1)
      {
        cout << "matching by distance jet" << endl;
      }

      float dphi = m_recojetphi - m_truthjetphi;
      if (dphi > TMath::Pi())
        dphi -= TMath::Pi() * 2.;
      if (dphi < -1 * TMath::Pi())
        dphi += TMath::Pi() * 2.;

      float deta = m_recojeteta - m_truthjeteta;
      /// Determine the distance in eta phi space between the reconstructed
      /// and truth jets
      m_dR = sqrt(pow(dphi, 2.) + pow(deta, 2.));

      /// If this truth jet is closer than the previous truth jet, it is
      /// closer and thus should be considered the truth jet
      if (m_dR < truth_jets->get_par() && m_dR < closestjet)
      {
        // Get reco jet characteristics
        m_recojetid = recoJet->get_id();
        m_recojetpx = recoJet->get_px();
        m_recojetpy = recoJet->get_py();
        m_recojetpz = recoJet->get_pz();
        m_recojetphi = recoJet->get_phi();
        m_recojetp = recoJet->get_p();
        m_recojetenergy = recoJet->get_e();
      }
    }

    /// Fill the truthjet tree
    m_truthjettree->Fill();
  }
}
*/
/**
 * Get the reconstructed jets and store them in a tree
 */
/*
void InttEventDisplay::getReconstructedJets(PHCompositeNode *topNode)
{
  /// Get the reconstructed tower jets
  JetMap *reco_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Tower_r04");
  /// Get the truth jets
  JetMap *truth_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Truth_r04");

  if (!m_jetEvalStack)
  {
    m_jetEvalStack = new JetEvalStack(topNode, "AntiKt_Tower_r04",
                                      "AntiKt_Truth_r04");
  }
  m_jetEvalStack->next_event(topNode);
  JetRecoEval *recoeval = m_jetEvalStack->get_reco_eval();
  if (!reco_jets)
  {
    cout << PHWHERE
         << "Reconstructed jet node is missing, can't collect reconstructed jets"
         << endl;
    return;
  }

  if (Verbosity() > 1)
  {
    cout << "Get all Reco Jets" << endl;
  }

  /// Iterate over the reconstructed jets
  for (JetMap::Iter recoIter = reco_jets->begin();
       recoIter != reco_jets->end();
       ++recoIter)
  {
    Jet *recoJet = recoIter->second;
    m_recojetpt = recoJet->get_pt();
    if (m_recojetpt < m_minjetpt)
      continue;

    m_recojeteta = recoJet->get_eta();

    // Get reco jet characteristics
    m_recojetid = recoJet->get_id();
    m_recojetpx = recoJet->get_px();
    m_recojetpy = recoJet->get_py();
    m_recojetpz = recoJet->get_pz();
    m_recojetphi = recoJet->get_phi();
    m_recojetp = recoJet->get_p();
    m_recojetenergy = recoJet->get_e();

    if (Verbosity() > 1)
    {
      cout << "matching by distance jet" << endl;
    }

    /// Set the matched truth jet characteristics to 0
    m_truthjetid = 0;
    m_truthjetp = 0;
    m_truthjetphi = 0;
    m_truthjeteta = 0;
    m_truthjetpt = 0;
    m_truthjetenergy = 0;
    m_truthjetpx = 0;
    m_truthjetpy = 0;
    m_truthjetpz = 0;

    Jet *truthjet = recoeval->max_truth_jet_by_energy(recoJet);
    if (truthjet)
    {
      m_truthjetid = truthjet->get_id();
      m_truthjetp = truthjet->get_p();
      m_truthjetpx = truthjet->get_px();
      m_truthjetpy = truthjet->get_py();
      m_truthjetpz = truthjet->get_pz();
      m_truthjeteta = truthjet->get_eta();
      m_truthjetphi = truthjet->get_phi();
      m_truthjetenergy = truthjet->get_e();
      m_truthjetpt = sqrt(m_truthjetpx * m_truthjetpx + m_truthjetpy * m_truthjetpy);
    }

    /// Check to make sure the truth jet node is available
    else if (truth_jets)
    {
      /// Match the reconstructed jet to the closest truth jet in delta R space
      /// Iterate over the truth jets
      float closestjet = 9999;
      for (JetMap::Iter truthIter = truth_jets->begin();
           truthIter != truth_jets->end();
           ++truthIter)
      {
        const Jet *truthJet = truthIter->second;

        float thisjetpt = truthJet->get_pt();
        if (thisjetpt < m_minjetpt)
          continue;

        float thisjeteta = truthJet->get_eta();
        float thisjetphi = truthJet->get_phi();

        float dphi = m_recojetphi - thisjetphi;
        if (dphi > TMath::Pi())
          dphi -= TMath::Pi() * 2.;
        if (dphi < -1. * TMath::Pi())
          dphi += TMath::Pi() * 2.;

        float deta = m_recojeteta - thisjeteta;
        /// Determine the distance in eta phi space between the reconstructed
        /// and truth jets
        m_dR = sqrt(pow(dphi, 2.) + pow(deta, 2.));

        /// If this truth jet is closer than the previous truth jet, it is
        /// closer and thus should be considered the truth jet

        if (m_dR < reco_jets->get_par() && m_dR < closestjet)
        {
          m_truthjetid = -9999;
          m_truthjetp = truthJet->get_p();
          m_truthjetphi = truthJet->get_phi();
          m_truthjeteta = truthJet->get_eta();
          m_truthjetpt = truthJet->get_pt();
          m_truthjetenergy = truthJet->get_e();
          m_truthjetpx = truthJet->get_px();
          m_truthjetpy = truthJet->get_py();
          m_truthjetpz = truthJet->get_pz();
          closestjet = m_dR;
        }
      }
    }
    m_recojettree->Fill();
  }
}
*/
/**
 * This method gets clusters from the EMCal and stores them in a tree. It
 * also demonstrates how to get trigger emulator information. Clusters from
 * other containers can be obtained in a similar way (e.g. clusters from
 * the IHCal, etc.)
 */
void InttEventDisplay::getEMCalClusters(PHCompositeNode *topNode)
{
  /// Get the raw cluster container
  /// Note: other cluster containers exist as well. Check out the node tree when
  /// you run a simulation
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");

  if (!clusters)
  {
    cout << PHWHERE
         << "EMCal cluster node is missing, can't collect EMCal clusters"
         << endl;
    return;
  }

  /// Get the global vertex to determine the appropriate pseudorapidity of the clusters
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    cout << "AnaTutorial::getEmcalClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    assert(vertexmap); // force quit

    return;
  }

  if (vertexmap->empty())
  {
    cout << "AnaTutorial::getEmcalClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    return;
  }

  GlobalVertex *vtx = vertexmap->begin()->second;
  if (vtx == nullptr)
    return;

  /// Trigger emulator
  CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");

  /// Can obtain some trigger information if desired
  if (trigger)
  {
    m_E_4x4 = trigger->get_best_EMCal_4x4_E();
  }
  RawClusterContainer::ConstRange begin_end = clusters->getClusters();
  RawClusterContainer::ConstIterator clusIter;

  /// Loop over the EMCal clusters
  for (clusIter = begin_end.first;
       clusIter != begin_end.second;
       ++clusIter)
  {
    /// Get this cluster
    const RawCluster *cluster = clusIter->second;

    /// Get cluster characteristics
    /// This helper class determines the photon characteristics
    /// depending on the vertex position
    /// This is important for e.g. eta determination and E_T determination
    CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
    m_clusenergy = E_vec_cluster.mag();
    m_cluseta = E_vec_cluster.pseudoRapidity();
    m_clustheta = E_vec_cluster.getTheta();
    m_cluspt = E_vec_cluster.perp();
    m_clusphi = E_vec_cluster.getPhi();

    if (m_cluspt < m_mincluspt)
      continue;

    m_cluspx = m_cluspt * cos(m_clusphi);
    m_cluspy = m_cluspt * sin(m_clusphi);
    m_cluspz = sqrt(m_clusenergy * m_clusenergy - m_cluspx * m_cluspx - m_cluspy * m_cluspy);

    // fill the cluster tree with all emcal clusters
    m_clustertree->Fill();
  }
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
  // get the geometry node
  /*
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container)
    {
      cout<<PHWHERE
    <<"CylinderGeomContainer node is missing."
    <<endl;
      return;
    }
  else cout<<"PH4CylinderGeomContainer is successfully loaded"<<endl;
  */
  /*
  m_geom = findNode::getClass<CylinderGeomIntt>(topNode, "CYLINDERGEOM_INTT");
  if (!m_geom)
    {
      cout<<PHWHERE
    <<"CylinderGeomIntt node is missing."
    <<endl;
      return;
    }
  */

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

std::vector<Acts::Vector3>  InttEventDisplay ::writeInttHits(PHCompositeNode * topNode)
{
  /* 
  std::vector<TrkrDefs::hitkey> matchedHits;
   std::vector<Acts::Vector3> hits;
  
  for (unsigned int inttlayer = 0; inttlayer < m_nInttLayers; inttlayer++)
      {
  for(const auto& hitkey : m_hitMap->getHitSets(TrkrDefs::TrkrId::inttId, inttlayer+3))
    {
      //cout <<"hitsetkey="<<hitsetkey<<endl;
      auto range = m_hitMap->getHits(hitkey);

      for (auto hitIter=range.first; hitIter!= range.second; ++hitIter)
        {
    const auto hitskey = hitIter->first;
    cout <<"hitkey="<<hitskey<<endl;
    const auto hit = hitIter->second;

    //cout<<"X = "<<hit->getLocalX()<<", Y = "<<hit->getLocalY()<<endl;

    //const auto globalPos = m_tGeometry->getGlobalPosition(cluskey, cluster);
    //cout<<"gX = "<<globalPos(0)<<", gY = "<<globalPos(1)<<" , gZ = "<<globalPos(2)<<endl;

    //clusters.push_back(globalPos);
	}
    }
      }

  */

  /*
  m_hitMap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  auto range = m_hitMap->getHitSets(TrkrDefs::TrkrId::inttId);
  //auto hits = m_hitMap->getHits();


  for (auto hitIter=range.first; hitIter!= range.second; ++hitIter)
    {
      const auto hitkey = hitIter->first;
      cout <<"hitkey="<<hitkey<<endl;
      auto hit = getHit(hitkey);

      const auto m_hits = hitIter->second;

      for(auto itr = m_hits.begin(); itr != m_hits.end();++itr ){
  auto hits = itr[0];
  cout <<"--------------"<<endl;
  cout<<"x="<<hits[0]<<"   y="<<hits[1]<<"   z="<<hits[2]<<endl;
  cout <<"--------------"<<endl<<endl<<endl<<endl<<endl;
      }

      //const auto globalPos = m_tGeometry->getGlobalPosition(cluskey, cluster);
      //cout<<"gX = "<<globalPos(0)<<", gY = "<<globalPos(1)<<" , gZ = "<<globalPos(2)<<endl;

      //clusters.push_back(globalPos);
    }
  */
  /*
  for(auto itr = m_hits.begin(); itr != m_hits.end();++itr ){
    auto hits = itr[0];
    cout <<"--------------"<<endl;
    //cout <<"itr = "  <<*itr << endl;
    cout<<"x="<<hits[0]<<"   y="<<hits[1]<<"   z="<<hits[2]<<endl;
    cout <<"--------------"<<endl<<endl<<endl<<endl<<endl;
  }
  */
  // return m_hits;

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
    //std::cout << "x: " << hiter->second->get_avg_x() << std::endl;
    glob(0)=hiter->second->get_avg_x();
    //std::cout << "y: " << hiter->second->get_avg_y() << std::endl;
    glob(1)=hiter->second->get_avg_y();
    //std::cout << "z: " << hiter->second->get_avg_z() << std::endl;
    glob(2)=hiter->second->get_avg_z();
    hits.push_back(glob);
  }
  std::cout << "InttG4HitRead::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  
  /*
  TrkrHitSetContainer *trkrhit = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkrhit)
  {
    std::cout << "Could not locate g4 hit node TRKR_HITSET " << std::endl;
    exit(1);
  }
  TrkrHitSetContainer::ConstRange hit_begin_end = trkrhit->getHitSets();
  for (TrkrHitSetContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    Acts::Vector3 glob;
    //std::cout << "x: " << hiter->second->get_avg_x() << std::endl;
    glob(0)=hiter->second->get_avg_x();
    //std::cout << "y: " << hiter->second->get_avg_y() << std::endl;
    glob(1)=hiter->second->get_avg_y();
    //std::cout << "z: " << hiter->second->get_avg_z() << std::endl;
    glob(2)=hiter->second->get_avg_z();
    hits.push_back(glob);
  }
  std::cout << "TRKRHITSET::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  */
  return hits;
}

std::vector<Acts::Vector3> InttEventDisplay ::writeInttClusters(PHCompositeNode * /*topNode*/)
{
  std::vector<TrkrDefs::cluskey> matchedClusters;
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
      //cout << "px =" << m_tr_px << endl;
      m_tr_py = track->get_py();
      //cout << "py =" << m_tr_py << endl;
      globt(1) = m_tr_py;
      m_tr_pz = track->get_pz();
      //cout << "pz =" << m_tr_pz << endl;
      globt(2) = m_tr_pz;
      tracks.push_back(globt);

      m_tr_p = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py + m_tr_pz * m_tr_pz);
      //cout << "m_tr_p =" << m_tr_p << endl;   
    }
    
  //Tracking = false;
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

      //m_tr_p = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py + m_tr_pz * m_tr_pz);
      //cout << "m_tr_p =" << m_tr_p << endl;
      
    }
  return vertex;
}


//vector<double> InttEventDisplay::writeMagnetField(PHCompositeNode *topNode)
//{
  //PHG4Reco *g4Reco = new PHG4Reco();
  /*
  default_config =  PHFieldUtility::DefaultFieldConfig();
  fieldmap = PHFieldUtility::GetFieldMapNode(default_config, topNode, Verbosity());
  //fieldmap = PHFieldUtility::GetFieldMapNode(default_config, topNode ,0);
  //default_config->identify();
  */
  /*
  if (fieldmap)
  {
    std::cout << "Found or created fieldmap" << std::endl;
  }
  else
  {
    std::cout << "Fieldmap not found or created" << std::endl;
  }
  */  

  //double field[3];
  /*
  //field[0]=fieldmap->get_field_mag_x();
  field[0]=default_config->get_field_mag_x();
  cout<<"mag field x ="<<field[0]<<endl;
  //field[1]=fieldmap->get_field_mag_y();
  field[1]=default_config->get_field_mag_y();
  cout<<"mag field y ="<<field[1]<<endl;
  //field[2]=fieldmap->get_field_mag_z();
  field[2]=default_config->get_field_mag_z();
  cout<<"mag field z ="<<field[2]<<endl;
  */
  /*
  //fieldmap->GetFieldValue(Point, field);
  g4Reco->PHG4MagneticField::GetFieldValue(Point, field);
  std::cout << "bx: " << field[0] / CLHEP::tesla
            << " by: " << field[1] / CLHEP::tesla
            << " bz: " << field[2] / CLHEP::tesla
            << std::endl;

  std::vector<double> bfield( std::begin(field), std::end(field) );
  */

  //return bfield;
//}

void InttEventDisplay::Display_3D()
{
  TEveManager::Terminate();
  TEveManager::Create();
  
  InttEventDisplay::drawall();

  // geometry load
  gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry.root");
  //InttEventDisplay::extractinttgeom();
  TEveGeoTopNode *geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  geom->CanEditMainTransparency();
  geom->SetMainTransparency(50);
  gEve->AddGlobalElement(geom);

  // x,y,z axis show
  TEveViewer *ev = gEve->GetDefaultViewer();
  v = ev->GetGLViewer();
  v->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

  PrintEidNclusters();

  // Camera control
  Double_t camcenter[3] = {0., 0., 0.};
  gSystem->ProcessEvents();
  v->CurrentCamera().RotateRad(0, -3.14 / 2);
  v->SetPerspectiveCamera(TGLViewer::kCameraPerspXOZ, 2, 30, &camcenter[0], 0, 0);

  v->RequestDraw();
  gEve->Redraw3D(kFALSE, kFALSE);
}

void InttEventDisplay::Display_rphi()
{
  TEveManager::Terminate();
  TEveManager::Create();
 
  InttEventDisplay::drawall();

  // geometry load
  gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry_rphi_black2mm.root");
  //InttEventDisplay::extractinttgeom();
  TEveGeoTopNode *geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  geom->CanEditMainTransparency();
  geom->SetMainTransparency(0);
  gEve->AddGlobalElement(geom);

  // camera
  TEveScene *s = gEve->SpawnNewScene("Projected Event");
  //gEve->GetDefaultViewer()->AddScene(s);
  v = gEve->GetDefaultGLViewer();
  Double_t camcenter[3] = {0., 0., 0.};
  v->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 5, 100, &camcenter[0], 0., 0.);
  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  v->SetClearColor(kWhite);

  InttEventDisplay::PrintEidNclusters();

  // projections
  TEveProjectionManager *mng =
      new TEveProjectionManager(TEveProjection::kPT_RPhi);
  s->AddElement(mng);
  TEveProjectionAxes *axes = new TEveProjectionAxes(mng);
  axes->SetTitle("TEveProjections demo");
  s->AddElement(axes);
  gEve->AddToListTree(axes, kTRUE);

  //scale
  TGLCameraOverlay* glovlay = v->GetCameraOverlay();
  glovlay->SetShowOrthographic(1);
  
  gEve->GetDefaultViewer()->AddScene(s);
  gEve->AddToListTree(mng, kTRUE);

  gEve->Redraw3D(kFALSE, kTRUE);
  v->RequestDraw();
}

void InttEventDisplay::Display_rhoz()
{
  TEveManager::Terminate();
  TEveManager::Create();
 
  InttEventDisplay::drawall();

  // geometry load
  gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry_rphi_black2mm.root");
  TEveGeoTopNode *geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  geom->CanEditMainTransparency();
  geom->SetMainTransparency(0);
  gEve->AddGlobalElement(geom);

  // camera
  TEveScene *s = gEve->SpawnNewScene("Projected Event");
  gEve->GetDefaultViewer()->AddScene(s);
  v = gEve->GetDefaultGLViewer();
  Double_t camcenter[3] = {0., 0., 0.};
  v->SetOrthoCamera(TGLViewer::kCameraOrthoZOY, 5, 100, &camcenter[0], 0., 0.);
  v->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
  v->SetClearColor(kWhite);

  InttEventDisplay::PrintEidNclusters();
  
  // projections
  TEveProjectionManager *mng =
      new TEveProjectionManager(TEveProjection::kPT_RhoZ);
  s->AddElement(mng);
  TEveProjectionAxes *axes = new TEveProjectionAxes(mng);
  axes->SetTitle("TEveProjections demo");
  s->AddElement(axes);
  gEve->AddToListTree(axes, kTRUE);
  gEve->AddToListTree(mng, kTRUE);

  //scale
  TGLCameraOverlay* glovlay =v->GetCameraOverlay();
  glovlay->SetShowOrthographic(1);

  gEve->Redraw3D(kFALSE, kFALSE);
  v->RequestDraw();
}

void InttEventDisplay::drawTracks()
{
  int numtrack = m_tracks.size();
  cout << "number of tracks = " << numtrack << endl;
  m_list = new TEveTrackList();

  if(numtrack==0){
    cout<<"No track!"<<endl;
    m_list = nullptr;
    return;
  }

  TEveTrackPropagator *prop = m_list->GetPropagator();
  TEveTrack *track[numtrack];
  TEveRecTrackD *rc = new TEveRecTrackD();
  int counter = 0;

  //prop->SetMagFieldObj(new TEveMagFieldConst(m_bfield[0],m_bfield[1], m_bfield[2]));
  //cout<<"(bfield[0],bfield[1], bfield[2])="<<m_bfield[0]<<","<<m_bfield[1]<<","<<","<<m_bfield[2]<<endl;
  prop->SetMagFieldObj(new TEveMagFieldConst(0.,0.,0.));

  m_list->SetElementName(Form("%s, constB", m_list->GetElementName()));
  m_list->SetLineColor(kMagenta);
  auto vertex = m_vertex[0];
  cout<<"vertex x =" << vertex[0]<<"vertex y =" << vertex[1]<<"vertex z =" << vertex[2]<<endl;
  for (auto itr = m_tracks.begin(); itr != m_tracks.end(); ++itr)
  {
    auto moment = itr[0];
    rc->fV.Set(vertex[0],vertex[1],vertex[2]);
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

void InttEventDisplay::drawClusters()
{
  int npoints = m_clusters.size();
  cout << "npoints = " << npoints << endl;
  m_ps = new TEvePointSet(npoints);

  if(npoints==0){
    cout<<"No Cluster!"<<endl;
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
    //cout<<"xyz = "<<cluster[0]<<","<< cluster[1]<<","<< cluster[2]<<endl;
    m_ps->SetPointId(new TNamed(Form("Point %d", counter), ""));
    counter++;
  }
  
  m_ps->SetMarkerColor(2);
  m_ps->SetMarkerSize(1.2);
  m_ps->SetMarkerStyle(4);
}

void InttEventDisplay::drawVertex()
{
  
  int nvertex = m_vertex.size();
  cout << "nvertex = " << nvertex << endl;
  //m_psv = new TEvePointSet(nvertex);

  if(nvertex==0){
    cout<<"No vertex!"<<endl;
    m_psv=nullptr;
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

void InttEventDisplay::drawHits()
{
  int nhits = m_hits.size();
  cout << "nhits = " << nhits << endl;
  m_psh = new TEvePointSet(nhits);

  if(nhits==0){
    cout<<"No hit!"<<endl;
    m_psh=nullptr;
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

void InttEventDisplay::getMatrices(string name, TGeoVolume *world, TGeoNode *node)
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
            //cout Volume name & TGeoMatrix & add node to intt volume
            //cout << "k =" << k << endl;
            //std::cout << "Volume name: " << volume->GetName() << std::endl;
            //std::cout << "(x, y, z): (" << matrix->GetTranslation()[0] << ", " << matrix->GetTranslation()[1] << ", " << matrix->GetTranslation()[2] << ")" << std::endl;
            //std::cout << "rotation:" << std::endl;
            //matrix->Print();
            //std::cout << std::endl;
            world->AddNode(volume, k, matrix);
            k++;
        }
    }
    for (Int_t i = 0; i < node->GetNdaughters(); ++i)
    {
        TGeoNode *daughterNode = node->GetDaughter(i);
        getMatrices(name, world, daughterNode);
    }
}

void InttEventDisplay::extractinttgeom()
{
    // get GeoMatrix from nodes
    TGeoNode *topNode = gGeoManager->GetTopNode();
    TGeoVolume *intt = gGeoManager->MakeBox("BOX",gGeoManager->GetTopVolume()->GetMedium(), 1, 1, 1);
    for (int inttlayer = 0; inttlayer < 4; inttlayer++)
    {
        for (int itype = 0; itype < 2; itype++)
        {
            string name = "ladder_" + to_string(inttlayer) + "_" + to_string(itype);
            getMatrices(name, intt, topNode);
        }
    }

    gGeoManager->SetTopVolume(intt);
    intt->SetVisibility(kFALSE);
    gGeoManager->CloseGeometry();    
}

void InttEventDisplay::PrintEidNclusters(){
  string legend = "eid :"+to_string(ievt-1)+", Nclusters: "+to_string(m_clusters.size());
  an = new TGLAnnotation(v,legend.c_str(),0.05,0.95);
  an->SetTextSize(0.05);
  an->SetTextColor(kWhite);
}

void InttEventDisplay::drawall(){
  
  InttEventDisplay::drawHits();
  if(m_psh!=nullptr){
    gEve->AddElement(m_psh);
  }

  InttEventDisplay::drawClusters();
  if(m_ps!=nullptr){
    gEve->AddElement(m_ps);
  }
  
  InttEventDisplay::drawTracks();
  if(m_list!=nullptr){
    gEve->AddElement(m_list);
  }
  
  InttEventDisplay::drawVertex();
  if(m_psv!=nullptr){
    gEve->AddElement(m_psv);
  }
}

void InttEventDisplay::Display_2D()
{
  /*
  TEveManager::Terminate();
  TEveManager::Create();
 
  InttEventDisplay::drawall();

  TEveScene *xyscene = gEve->SpawnNewScene("xy view");
  gEve->GetDefaultViewer()->AddScene(xyscene);
  //TGLViewer *xyviewer = gEve->GetDefaultGLViewer();
  Double_t camcenter[3] = {0., 0., 0.};
  //xyviewer->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 5, 100, &camcenter[0], 0., 0.);
  //xyviewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);

  TEveProjectionManager * xyproman =
      new TEveProjectionManager(TEveProjection::kPT_RPhi);
  //xys->AddElement(xyman);

  {
         TEveProjectionAxes* a = new TEveProjectionAxes(xyproman);
         a->SetMainColor(kWhite);
         a->SetTitle("xyview");
         a->SetTitleSize(0.05);
         a->SetTitleFont(102);
         a->SetLabelSize(0.025);
         a->SetLabelFont(102);
         xyscene->AddElement(a);
  }
  gEve->AddToListTree(xyproman, kFALSE);

  TEveWindowSlot *slot = 0;
  slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  TEveViewer * xyviewer = gEve->SpawnNewViewer("xyview","");
  xyviewer = GetGLViewer()->SetOrthoCamera(TGLViewer::kCameraOrthoXOY, 5, 100, &camcenter[0], 0., 0.);
  xyviewer -> GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  xyviewer -> AddScene(xyscene);
  */
}
