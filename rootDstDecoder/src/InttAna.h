#ifndef INTTANA_H__
#define INTTANA_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <fstream>

/// Class declarations for use in the analysis module
class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TNtuple;

/// Definition of this analysis module class
class InttAna : public SubsysReco
{
 public:
  /// Constructor
  InttAna(const std::string &name = "InttAna",
          const std::string &fname = "AnaTutorial.root");

  // Destructor
  virtual ~InttAna();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);


  void readRawHit(PHCompositeNode *);

 private:
  std::string fname_;
  TFile* anafile_;
  TH1*   h_dca2d_zero; 
  TH2*   h2_dca2d_zero; 
  TH2*   h2_dca2d_len; 
  TNtuple*  h_ntp_clus; 
  TNtuple*  h_ntp_cluspair; 

};
#endif
