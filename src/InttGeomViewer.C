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

using namespace std;

void InttGeomViewer(string fname = "inttgeometry_rphi_black2mm.root")
{ 
  TEveManager::Create();

  gGeoManager = gEve->GetGeometry(fname.c_str());
  TEveGeoTopNode *geom = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  geom->CanEditMainTransparency();
  geom->CanEditMainColor();
  geom->SetMainColor(kRed);
  geom->SetMainTransparency(50);
  gEve->AddGlobalElement(geom);

  TEveViewer *ev = gEve->GetDefaultViewer();
  TGLViewer *gv = ev->GetGLViewer();
  gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  gEve->Redraw3D(kTRUE);
  //set background color
  gv->SetClearColor(kWhite);

  gSystem->ProcessEvents();                    
//gv->CurrentCamera().RotateRad(0, -3.14 / 2); 
  gv->RequestDraw();
}
