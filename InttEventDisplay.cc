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
#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

/// Tracking includes
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterv4.h>            
#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include<Acts/Surfaces/Surface.hpp>
#include<Acts/Surfaces/CylinderSurface.hpp>
#include<trackreco/PHActsSiliconSeeding.h>

/// Truth evaluation includes
#include <g4eval/JetEvalStack.h>
#include <g4eval/SvtxEvalStack.h>

/// HEPMC truth includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored"-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
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
  : SubsysReco(name)
  , m_outfilename(filename)
  , m_hm(nullptr)
  , m_minjetpt(5.0)
  , m_mincluspt(0.25)
  , m_analyzeTracks(true)
  , m_analyzeClusters(true)
  , m_analyzeJets(true)
  , m_analyzeTruth(false)
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
  if(m_c1!=NULL) delete m_c1;
}

/**
 * Initialize the module and prepare looping over events
 */
int InttEventDisplay::Init(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 5)
  {
    cout << "Beginning Init in InttEventDisplay" << endl;
  }
 
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");

  m_phi_h = new TH1D("phi_h", ";Counts;#phi [rad]", 50, -6, 6);
  m_eta_phi_h = new TH2F("phi_eta_h", ";#eta;#phi [rad]", 10, -1, 1, 50, -6, 6);

  //cout<<"print print print"<<endl;
  //topNode->print();
  //cout<<"pripripirpint"<<endl;
  cout<<"finish intteventdisplay init"<<endl;
  return 0;
}

/**
 * Main workhorse function where each event is looped over and 
 * data from each event is collected from the node tree for analysis
 */
int InttEventDisplay::process_event(PHCompositeNode *topNode)
{

  cout << "process_event : start"<<endl<<endl<<endl<<endl<<endl;
  cout << "process_event : start"<<endl<<endl<<endl<<endl<<endl;
  cout << "process_event : start"<<endl<<endl<<endl<<endl<<endl;
  cout << "process_event : start"<<endl<<endl<<endl<<endl<<endl;
  cout << "process_event : start"<<endl<<endl<<endl<<endl<<endl;
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  CylinderGeomIntt* geom = dynamic_cast<CylinderGeomIntt*>(geom_container->GetLayerGeom(3));
  if (geom==NULL) cout<<"No CylinderGeomIntt"<<endl<<endl<<endl<<endl;
  float pitch = geom->get_strip_y_spacing();
  float length = geom->get_strip_z_spacing();
  float radius = geom->get_radius();
  //int layer = geom->get_layer();
  // float tilt = geom->get_strip_tilt();

  cout<<"pitch : "<<pitch<<"   length:"<<length<<"    radius:"<<radius<<endl;
  // cout<<"layer : "<<layer<<endl;
  //cout<<"strip : "<<strip<<endl;


  if (Verbosity() > 5)
  {
    cout << "Beginning process_event in AnaTutorial" << endl;
  }
  /// Get the truth information
  if (m_analyzeTruth && 0)
  {
    getHEPMCTruth(topNode);
    cout<<"getHEPMCTruth"<<endl;
    getPHG4Truth(topNode);
    cout<<"getPHG4Truth"<<endl;
  }

  /// Get the tracks
  if (m_analyzeTracks && 0)
  {
    getTracks(topNode);
    cout<<"gettracks"<<endl;
    
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

  //Get some Nodes. INTT geometry + hit container
  cout<<"bofore getnode"<<endl;
  getNode(topNode);
  cout<<"getnode"<<endl;

  //std::vector<Acts::Vector3> 
  m_clusters = writeInttClusters(topNode);
  cout<<"number of clusters is "<<m_clusters.size()<<endl<<endl<<endl<<endl;
  m_tracks = writeInttTracks(topNode);
  cout<<"number of tracks is "<<m_tracks.size()<<endl<<endl<<endl<<endl;

  for(auto itr = m_tracks.begin(); itr != m_tracks.end();++itr ){
    auto track = itr[0];
    cout <<"--------------"<<endl;	     
    cout <<"itr = "  <<*itr << endl;
    cout<<"px="<<track[0]<<"   py="<<track[1]<<"   pz="<<track[2]<<endl;
    cout <<"--------------"<<endl<<endl<<endl<<endl<<endl;
  }

  
  /*
  for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
    auto cluster = itr[0];
    cout <<"--------------"<<endl;	     
    cout <<"itr = "  <<*itr << endl;
    cout<<"cluster[0]="<<cluster[0]<<"   cluster[1]="<<cluster[1]<<"   cluster[2]="<<cluster[2]<<endl;
    //cout <<"--------------"<<endl<<endl<<endl<<endl<<endl;
  }
  */
  

  cout<<"number of points is "<<m_clusters.size()<<endl<<endl<<endl<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int InttEventDisplay::End(PHCompositeNode */*topNode*/)
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

  return 0;}

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
  if(!m_svtxEvalStack)
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

    count1 = count1 +1;
    m_tracktree->Fill();
  }
  cout<<"get track loop ="<<count1<<endl;
}

/**
 * Method that gets the truth jets and stores them in a tree
 */
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
  if(!m_jetEvalStack)
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
    //loop over the constituents of the truth jet
    for (std::set<PHG4Particle *>::iterator iter2 = truthjetcomp.begin();
         iter2 != truthjetcomp.end();
         ++iter2)
    {
      //get the particle of the truthjet
      PHG4Particle *truthpart = *iter2;
      if (!truthpart)
      {
        cout << "no truth particles in the jet??" << endl;
        break;
      }

      ntruthconstituents++;
    }

    if(ntruthconstituents < 3)
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
	if (m_recojetpt < m_minjetpt-m_minjetpt*0.4)
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

/**
 * Get the reconstructed jets and store them in a tree
 */
void InttEventDisplay::getReconstructedJets(PHCompositeNode *topNode)
{
  /// Get the reconstructed tower jets
  JetMap *reco_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Tower_r04");
  /// Get the truth jets
  JetMap *truth_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Truth_r04");

  if(!m_jetEvalStack)
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
    if(truthjet){
      m_truthjetid = truthjet->get_id();
      m_truthjetp = truthjet->get_p();
      m_truthjetpx = truthjet->get_px();
      m_truthjetpy = truthjet->get_py();
      m_truthjetpz = truthjet->get_pz();
      m_truthjeteta = truthjet->get_eta();
      m_truthjetphi = truthjet->get_phi();
      m_truthjetenergy = truthjet->get_e();
      m_truthjetpt = sqrt(m_truthjetpx * m_truthjetpx 
			  + m_truthjetpy * m_truthjetpy);
    }

    /// Check to make sure the truth jet node is available
    else if(truth_jets)
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
        if (dphi >TMath::Pi())
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
    assert(vertexmap);  // force quit

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
  if(trigger)
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

    //fill the cluster tree with all emcal clusters
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
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << "No ActsGeometry on node tree. Bailing."
		<< std::endl;
      return;
    }
  else cout<<"ActGeometry on node tree is successfully loaded"<<endl;


  if(m_useTruthClusters){
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							    "TRKR_CLUSTER_TRUTH");
    cout<<"TRKR_CLUSTER_TRUTH"<<endl;
  }
  else{
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,
							     "TRKR_CLUSTER");
    cout<<"TRKR_CLUSTER"<<endl;
  }

  if(!m_clusterMap)
    {
      cout<<PHWHERE
	  <<"TrkrClusterContainer node is missing."
	  <<endl;
      //return ABORTEVENT;
      return;
    }
  else cout<<"TrkrClusterContainr is successfully loaded"<<endl;

  return;
}


std::vector<Acts::Vector3> InttEventDisplay :: writeInttClusters(PHCompositeNode */*topNode*/)
{
    std::vector<TrkrDefs::cluskey> matchedClusters;
    std::vector<Acts::Vector3> clusters;
    for (unsigned int inttlayer = 0; inttlayer < m_nInttLayers; inttlayer++)
      {
	for( const auto& hitsetkey : m_clusterMap->getHitSetKeys(TrkrDefs::TrkrId::inttId, inttlayer+3))
	  {
	    //cout <<"hitsetkey="<<hitsetkey<<endl;
	    auto range = m_clusterMap->getClusters(hitsetkey);
	    

	    for (auto clusIter=range.first; clusIter!= range.second; ++clusIter)
	      {
		const auto cluskey = clusIter->first;
		//cout <<"cluskey="<<cluskey<<endl;
		const auto cluster = clusIter->second;

		//cout<<"X = "<<cluster->getLocalX()<<", Y = "<<cluster->getLocalY()<<endl;

		const auto globalPos = m_tGeometry->getGlobalPosition(cluskey, cluster);
		//cout<<"gX = "<<globalPos(0)<<", gY = "<<globalPos(1)<<" , gZ = "<<globalPos(2)<<endl;

		clusters.push_back(globalPos);
              }
	  }
      }
    

    //get segments locations
    /*
    uint8_t type;
    uint8_t ladder_index_max;
    int time_bucket = 0; 

    for(uint8_t lyr=3;lyr<7;lyr++){
      
      //number of ladders
      switch(lyr){
      case 3:
	ladder_index_max=12;
	cout<<"ladder_index_max ="<<ladder_index_max<<endl;
	break;
      case 4:
	ladder_index_max=12;
	cout<<"ladder_index_max ="<<ladder_index_max<<endl;
	break;
      case 5:
	ladder_index_max=16;
	cout<<"ladder_index_max ="<<ladder_index_max<<endl;
	break;
      case 6:
	ladder_index_max=16;
	cout<<"ladder_index_max ="<<ladder_index_max<<endl;
	break;
      }
    

      //typeA position
      type = 0;
      for(uint8_t ladder_index=0;ladder_index<ladder_index_max;ladder_index++){
	auto genhitkeyintt = InttDefs::genHitSetKey(lyr,type,ladder_index,time_bucket);
	auto surfintt = m_tGeometry->maps().getSiliconSurface(genhitkeyintt);
	double ladderLocation[3] = {0.,0.,0.};
	m_geom->find_segment_center(surfintt,m_tGeometry,ladderLocation);
	//cout << "ladderLocationA=(" << ladderLocation[0] << "," << ladderLocation[1] << "," << ladderLocation[2] << ")" << endl;

	ladderALocationX.push_back(ladderLocation[0]);
	ladderALocationY.push_back(ladderLocation[1]);
	ladderALocationZ.push_back(ladderLocation[2]);
	//cout<<"(lyr,ladder_index)=("<<unsigned(lyr)<<","<<unsigned(ladder_index)<<")"<<endl;
      }

      type = 2;
      for(uint8_t ladder_index=0;ladder_index<ladder_index_max;ladder_index++){
	auto genhitkeyintt = InttDefs::genHitSetKey(lyr,type,ladder_index,time_bucket);
	auto surfintt = m_tGeometry->maps().getSiliconSurface(genhitkeyintt);
	double ladderLocation[3] = {0.,0.,0.};
	m_geom->find_segment_center(surfintt,m_tGeometry,ladderLocation);
	//cout << "ladderLocationA=(" << ladderLocation[0] << "," << ladderLocation[1] << "," << ladderLocation[2] << ")" << endl;

	ladderALocationX.push_back(ladderLocation[0]);
	ladderALocationY.push_back(ladderLocation[1]);
	ladderALocationZ.push_back(ladderLocation[2]);
	//cout<<"(lyr,ladder_index)=("<<unsigned(lyr)<<","<<unsigned(ladder_index)<<")"<<endl;
      }

      //typeB position
      type = 1;
      for(uint8_t ladder_index=0;ladder_index<ladder_index_max;ladder_index++){
	auto genhitkeyintt = InttDefs::genHitSetKey(lyr,type,ladder_index,time_bucket);
	auto surfintt = m_tGeometry->maps().getSiliconSurface(genhitkeyintt);
	double ladderLocation[3] = {0.,0.,0.};
	m_geom->find_segment_center(surfintt,m_tGeometry,ladderLocation);
	//cout << "ladderLocationB=(" << ladderLocation[0] << "," << ladderLocation[1] << "," << ladderLocation[2] << ")" << endl;

	ladderBLocationX.push_back(ladderLocation[0]);
	ladderBLocationY.push_back(ladderLocation[1]);
	ladderBLocationZ.push_back(ladderLocation[2]);
	//cout<<"(lyr,ladder_index)=("<<unsigned(lyr)<<","<<unsigned(ladder_index)<<")"<<endl;
      }


      type = 3;
      for(uint8_t ladder_index=0;ladder_index<ladder_index_max;ladder_index++){
	auto genhitkeyintt = InttDefs::genHitSetKey(lyr,type,ladder_index,time_bucket);
	auto surfintt = m_tGeometry->maps().getSiliconSurface(genhitkeyintt);
	double ladderLocation[3] = {0.,0.,0.};
	m_geom->find_segment_center(surfintt,m_tGeometry,ladderLocation);
	//cout << "ladderLocationB=(" << ladderLocation[0] << "," << ladderLocation[1] << "," << ladderLocation[2] << ")" << endl;

	ladderBLocationX.push_back(ladderLocation[0]);
	ladderBLocationY.push_back(ladderLocation[1]);
	ladderBLocationZ.push_back(ladderLocation[2]);
	//cout<<"(lyr,ladder_index)=("<<unsigned(lyr)<<","<<unsigned(ladder_index)<<")"<<endl;
      }
      
    

      }
    */

      return clusters;
}

std::vector<Acts::Vector3> InttEventDisplay :: writeInttTracks(PHCompositeNode *topNode){
  std::vector<Acts::Vector3> tracks;
    /// SVTX tracks node
    SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    /*
    if (!trackmap)
    {
        cout << PHWHERE
             << "SvtxTrackMap node is missing, can't collect tracks"
             << endl;
    }
    */

   /// EvalStack for truth track matching
    if (!m_svtxEvalStack)
    {
        m_svtxEvalStack = new SvtxEvalStack(topNode);
        m_svtxEvalStack->set_verbosity(Verbosity());
    }

    m_svtxEvalStack->next_event(topNode);

    /// Get the track evaluator
    // SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();

    /// Get the range for primary tracks
    // PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    Acts::Vector3 glob;
    if (Verbosity() > 1)
    {
        cout << "Get the SVTX tracks" << endl;
    }
    for (unsigned int inttlayer = 0; inttlayer < m_nInttLayers; inttlayer++)
    {
      //for (const auto &hitsetkey : m_clusterMap->getHitSetKeys(TrkrDefs::TrkrId::inttId, inttlayer + 3)){
            // cout <<"hitsetkey="<<hitsetkey<<endl;
            //auto range = m_clusterMap->getClusters(hitsetkey);

            //for (auto clusIter = range.first; clusIter != range.second; ++clusIter){
                for (SvtxTrackMap::Iter iter = trackmap->begin();iter != trackmap->end();++iter)
                {
                    SvtxTrack *track = iter->second;

                    /// Get the reconstructed track info
                    m_tr_px = track->get_px();
                    glob(0) = m_tr_px;
                    m_tr_py = track->get_py();
                    glob(1) = m_tr_py;
                    m_tr_pz = track->get_pz();
                    glob(2) = m_tr_pz;
                    tracks.push_back(glob);
                }
		//}
    
    }
    Tracking = true;
    return tracks;
    }

void InttEventDisplay :: make_ladderlocationfile(){
  TFile * location = TFile::Open("/sphenix/u/mfujiwara/Documents/segmentlocation.root","RECREATE");

  location->WriteObject(&ladderALocationX,"segmentALocationX");
  location->WriteObject(&ladderALocationY,"segmentALocationY");
  location->WriteObject(&ladderALocationZ,"segmentALocationZ");

  location->WriteObject(&ladderBLocationX,"segmentBLocationX");
  location->WriteObject(&ladderBLocationY,"segmentBLocationY");
  location->WriteObject(&ladderBLocationZ,"segmentBLocationZ");

  location->Close();

  cout<<"finish process"<<endl;
}

void InttEventDisplay :: DrawHit_rphi(){
  //TFile::SetCacheFileDir(".");
  TEveManager::Terminate();
   TEveManager::Create();

   //open geom file
   gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry_rphi.root");
   TEveGeoTopNode *geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
   gEve->AddGlobalElement(geom);

   /*
   // camera
   TEveScene* s = gEve->SpawnNewScene("Projected Event");
   gEve->GetDefaultViewer()->AddScene(s);
   TGLViewer* v = gEve->GetDefaultGLViewer();
   v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   Double_t camcenter[3]={0.,0.,0.};
   v->SetOrthoCamera(TGLViewer::kCameraOrthoXOY,15,50,&camcenter[0],0.,0.);
   v->RequestDraw();
   //TGLOrthoCamera& cam = (TGLOrthoCamera&) v->CurrentCamera();
   //cam.SetZoomMinMax(5.0,10.);
   
   // projections
   TEveProjectionManager* mng =
      new TEveProjectionManager(TEveProjection::kPT_RPhi);
   s->AddElement(mng);
   TEveProjectionAxes* axes = new TEveProjectionAxes(mng);
   axes->SetTitle("TEveProjections demo");
   s->AddElement(axes);
   gEve->AddToListTree(axes, kTRUE);
   gEve->AddToListTree(mng, kTRUE);
   */

   //hit point
   int npoints = m_clusters.size();
    cout<<"npoints = " <<npoints<<endl;
    TEvePointSet* ps = new TEvePointSet(npoints);
    ps->SetOwnIds(kTRUE);

    cout<<"new TEvePointSet"<<endl;

    int counter  = 0;
    for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
      auto cluster = itr[0];
      ps->SetNextPoint(cluster[0], cluster[1], cluster[2]);
      ps->SetPointId(new TNamed(Form("Point %d", counter), ""));
      counter ++ ;
      //cout<<"itr = "<<*itr<<endl;
      //cout<<"cluster[0] = "<<cluster[0]<<"   cluster[1] = "<<cluster[1]<<"     cluster[2] = "<<cluster[2]<<endl;
    }
 
    ps->SetMarkerColor(2);
    ps->SetMarkerSize(1.2);
    ps->SetMarkerStyle(4);
    gEve->AddElement(ps);

    geom->CanEditMainColor();

    //gEve->Redraw3D(kTRUE);

    if(Tracking == true){
    //set Magfield and draw Tracks
    int numtrack = m_tracks.size();
    cout<<"number of tracks = " <<numtrack<<endl;
    
    TEveTrackList *list = new TEveTrackList();
    TEveTrackPropagator *prop = list->GetPropagator();
    TEveTrack *track[numtrack];
    TEveRecTrackD *rc = new TEveRecTrackD();
    counter = 0;

    prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 0.));

    list->SetElementName(Form("%s, constB", list->GetElementName()));
    list->SetLineColor(kMagenta);
    for (auto itr = m_tracks.begin(); itr != m_tracks.end();++itr)
    {   
        auto moment = itr[0];
        rc->fV.Set(0.0, 0.0, 0.0);
        rc->fP.Set(moment[0],moment[1],moment[2]);
        rc->fSign = 1;
        track[counter] = new TEveTrack(rc, prop);
	list->AddElement(track[counter]);
	track[counter]->SetLineColor(list->GetLineColor());
	track[counter]->MakeTrack();
        counter ++;
    }
    
    gEve->AddElement(list);
    //gEve->Redraw3D(kTRUE);
    //v->RequestDraw();

    cout<<"number of hit points = " <<npoints<<endl;
    cout<<"number of tracks ="<<numtrack<<endl;
    }

    // camera
   TEveScene* s = gEve->SpawnNewScene("Projected Event");
   gEve->GetDefaultViewer()->AddScene(s);
   TGLViewer* v = gEve->GetDefaultGLViewer();
   //v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   Double_t camcenter[3]={0.,0.,0.};
   v->SetOrthoCamera(TGLViewer::kCameraOrthoXOY,5,100,&camcenter[0],0.,0.);
   //v->RequestDraw();
   v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   //TGLOrthoCamera& cam = (TGLOrthoCamera&) v->CurrentCamera();
   //cam.SetZoomMinMax(5.0,10.);
   
   // projections
   TEveProjectionManager* mng =
      new TEveProjectionManager(TEveProjection::kPT_RPhi);
   //mng->SetCenter(0,0,0);
   s->AddElement(mng);
   TEveProjectionAxes* axes = new TEveProjectionAxes(mng);
   axes->SetTitle("TEveProjections demo");
   s->AddElement(axes);
   gEve->AddToListTree(axes, kTRUE);
   gEve->AddToListTree(mng, kTRUE);
   
   gEve->Redraw3D(kFALSE,kFALSE);
   v->RequestDraw();  
}

void InttEventDisplay :: DrawHits()
{
  TEveManager::Terminate();
    TEveManager::Create();

    //cout<<"TEve created"<<endl;

    int npoints = m_clusters.size();
    cout<<"npoints = " <<npoints<<endl;
    TEvePointSet* ps = new TEvePointSet(npoints);
    ps->SetOwnIds(kTRUE);

    //cout<<"new TEvePointSet"<<endl;

    int counter  = 0;
    for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
      auto cluster = itr[0];
      ps->SetNextPoint(cluster[0], cluster[1], cluster[2]);
      ps->SetPointId(new TNamed(Form("Point %d", counter), ""));
      counter ++ ;
      //cout<<"itr = "<<*itr<<endl;
      //cout<<"cluster[0] = "<<cluster[0]<<"   cluster[1] = "<<cluster[1]<<"     cluster[2] = "<<cluster[2]<<endl;
    }
 
    ps->SetMarkerColor(2);
    ps->SetMarkerSize(1.2);
    ps->SetMarkerStyle(4);

    //typeA center location
    int Apoints = ladderALocationX.size();
    //cout <<"Apoints ="<< Apoints <<endl;
    TEvePointSet* ps2 = new TEvePointSet(Apoints);
    ps2->SetOwnIds(kTRUE);

    for(int i=0;i<Apoints;i++){
      ps2->SetNextPoint(ladderALocationX[i],ladderALocationY[i],ladderALocationZ[i]);
      ps2->SetPointId(new TNamed(Form("Point %d", i), ""));
    }
    
    ps2->SetMarkerColor(3);
    ps2->SetMarkerSize(1.0);
    ps2->SetMarkerStyle(4);
 
    //typeB center location
    int Bpoints = ladderBLocationX.size();
    //cout <<"Bpoints ="<< Bpoints <<endl;
    TEvePointSet* ps3 = new TEvePointSet(Bpoints);
    ps3->SetOwnIds(kTRUE);

    for(int i=0;i<Bpoints;i++){
      ps3->SetNextPoint(ladderBLocationX[i],ladderBLocationY[i],ladderBLocationZ[i]);
      ps3->SetPointId(new TNamed(Form("Point %d", i), ""));
    }
    
    ps3->SetMarkerColor(5);
    ps3->SetMarkerSize(1.0);
    ps3->SetMarkerStyle(4);

    /*
    if (parent)
    {
       parent->AddElement(ps);
    }
    */
    //.X else
    
       gEve->AddElement(ps);
       //gEve->AddElement(ps2);
       // gEve->AddElement(ps3);

       //set Magfield and draw Tracks
       if(Tracking == true){
       int numtrack = m_tracks.size();
       cout<<"number of tracks = " <<numtrack<<endl;
    
       TEveTrackList *list = new TEveTrackList();
       TEveTrackPropagator *prop = list->GetPropagator();
       TEveTrack *track[numtrack];
       TEveRecTrackD *rc = new TEveRecTrackD();
       counter = 0;

       prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 0.));

       list->SetElementName(Form("%s, constB", list->GetElementName()));
       list->SetLineColor(kMagenta);
       for (auto itr = m_tracks.begin(); itr != m_tracks.end();++itr)
	 {   
	   auto moment = itr[0];
	   rc->fV.Set(0.0, 0.0, 0.0);
	   rc->fP.Set(moment[0],moment[1],moment[2]);
	   rc->fSign = 1;
	   track[counter] = new TEveTrack(rc, prop);
	   list->AddElement(track[counter]);
	   track[counter]->SetLineColor(list->GetLineColor());
	   track[counter]->MakeTrack();
	   counter ++;
	 }
    
       cout<<"counter ="<<counter<<endl;
       gEve->AddElement(list);
       //gEve->Redraw3D(kTRUE);
       }
       
       //geometry load
       gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry.root");
       TEveGeoTopNode* geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
       geom->CanEditMainTransparency();
       geom->SetMainTransparency(50);
       gEve->AddGlobalElement(geom);
       

       //x,y,z axis show
       TEveViewer *ev = gEve->GetDefaultViewer();
       TGLViewer *gv = ev->GetGLViewer();
       gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0); 
       //gEve->Redraw3D(kTRUE);

       //Camera control
       Double_t camcenter[3]={0.,0.,0.};
       gSystem->ProcessEvents();
       gv->CurrentCamera().RotateRad(0,-3.14/2);
       gv->SetPerspectiveCamera(TGLViewer::kCameraPerspXOZ,2,30,&camcenter[0],0,0);
       //TGLPerspectiveCamera& gpc =(TGLPerspectiveCamera&) gv->CurrentCamera();
       //gpc.Zoom(10,kTRUE,kTRUE);

       gv->RequestDraw();
       gEve->Redraw3D(kFALSE,kFALSE);
       cout<<"number of hit points = " <<npoints<<endl;
    //return ps;
}

void InttEventDisplay::Draw_rhoz_display(){
  TEveManager::Terminate();
  TEveManager::Create();
  
  InttEventDisplay::drawhits();
  gEve->AddElement(m_ps);
  InttEventDisplay::DrawTracks();
  gEve->AddElement(m_list);

  //geometry load
  gGeoManager = gEve->GetGeometry("/sphenix/u/mfujiwara/Documents/inttgeometry.root");
  TEveGeoTopNode* geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  geom->CanEditMainTransparency();
  geom->SetMainTransparency(20);
  gEve->AddGlobalElement(geom);

  // camera
   TEveScene* s = gEve->SpawnNewScene("Projected Event");
   gEve->GetDefaultViewer()->AddScene(s);
   TGLViewer* v = gEve->GetDefaultGLViewer();
   Double_t camcenter[3]={0.,0.,0.};
   v->SetOrthoCamera(TGLViewer::kCameraOrthoZOY,5,100,&camcenter[0],0.,0.);
   v->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   
   // projections
   TEveProjectionManager* mng =
      new TEveProjectionManager(TEveProjection::kPT_RhoZ);
   s->AddElement(mng);
   TEveProjectionAxes* axes = new TEveProjectionAxes(mng);
   axes->SetTitle("TEveProjections demo");
   s->AddElement(axes);
   gEve->AddToListTree(axes, kTRUE);
   gEve->AddToListTree(mng, kTRUE);
   
   gEve->Redraw3D(kFALSE,kFALSE);
   v->RequestDraw();
   
   
}

void InttEventDisplay::DrawTracks(){
  int numtrack = m_tracks.size();
       cout<<"number of tracks = " <<numtrack<<endl;
    
       m_list = new TEveTrackList();
       TEveTrackPropagator *prop = m_list->GetPropagator();
       TEveTrack *track[numtrack];
       TEveRecTrackD *rc = new TEveRecTrackD();
       int counter = 0;

       prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 0.));

       m_list->SetElementName(Form("%s, constB", m_list->GetElementName()));
       m_list->SetLineColor(kMagenta);
       for (auto itr = m_tracks.begin(); itr != m_tracks.end();++itr)
	 {   
	   auto moment = itr[0];
	   rc->fV.Set(0.0, 0.0, 0.0);
	   rc->fP.Set(moment[0],moment[1],moment[2]);
	   rc->fSign = 1;
	   track[counter] = new TEveTrack(rc, prop);
	   m_list->AddElement(track[counter]);
	   track[counter]->SetLineColor(m_list->GetLineColor());
	   track[counter]->MakeTrack();
	   counter ++;
	 }
    
       cout<<"counter ="<<counter<<endl;
       //gEve->AddElement(list);

       //return list;
}

void InttEventDisplay::drawhits(){
int npoints = m_clusters.size();
    cout<<"npoints = " <<npoints<<endl;
    m_ps = new TEvePointSet(npoints);
    m_ps->SetOwnIds(kTRUE);

    //cout<<"new TEvePointSet"<<endl;

    int counter  = 0;
    for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
      auto cluster = itr[0];
      m_ps->SetNextPoint(cluster[0], cluster[1], cluster[2]);
      m_ps->SetPointId(new TNamed(Form("Point %d", counter), ""));
      counter ++ ;
      //cout<<"itr = "<<*itr<<endl;
      //cout<<"cluster[0] = "<<cluster[0]<<"   cluster[1] = "<<cluster[1]<<"     cluster[2] = "<<cluster[2]<<endl;
    }
 
    m_ps->SetMarkerColor(2);
    m_ps->SetMarkerSize(1.2);
    m_ps->SetMarkerStyle(4);

    //return ps;
}



void InttEventDisplay::drawCanvas()
{
  if(m_c1==NULL){
    m_c1 = new TCanvas("c1", "c1", 500, 500);
    m_c1->Range(-100, -100, 100, 100);
  }
}

void InttEventDisplay::drawHit()
{
  if(m_c1==NULL) return;

  clear();



    int npoints = m_clusters.size();
    cout<<"npoints = " <<npoints<<endl;
    TEvePointSet* ps = new TEvePointSet(npoints);
    ps->SetOwnIds(kTRUE);

    cout<<"new TEvePointSet"<<endl;

    std::vector<double> gX;
    std::vector<double> gY;
    for(auto itr = m_clusters.begin(); itr != m_clusters.end();++itr ){
      auto cluster = itr[0];
      gX.push_back(cluster[0]);
      gY.push_back(cluster[1]);
    }

    TMarker *pos = new TMarker(gX[0], gY[0], gX.size());
    pos->SetMarkerColor(2);
    pos->Draw("same");

  vPos.push_back(pos);
}

void InttEventDisplay::clear()
{
  for(int i=0; i<(int)vPos.size(); i++){
    delete vPos[i];
  }
  vPos.clear();
}
