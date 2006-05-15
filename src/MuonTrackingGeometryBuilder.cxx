///////////////////////////////////////////////////////////////////
// MuonTrackingGeometryBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonTrackingGeometryBuilder.h"
// Trk
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeHelper.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeDisplayer.h"
#include "TrkDetDescrUtils/BinUtility1DR.h"
#include "TrkDetDescrUtils/BinUtility1DZ.h"
#include "TrkDetDescrUtils/BinnedArray.h"
//#include "TrkDetDescrUtils/BinningType.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
//#include "TrkGeometry/SimplifiedMaterialProperties.h"
#include "TrkMagFieldTools/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/TrackingGeometry.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/GlueVolumesDescriptor.h"
#include<fstream>

// BField
#include "BFieldAth/MagFieldAthena.h"

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

// STD
#include <map>

// Gaudi
#include "GaudiKernel/MsgStream.h"


// constructor
Muon::MuonTrackingGeometryBuilder::MuonTrackingGeometryBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgTool(t,n,p),
  m_magFieldTool(0),
  m_magFieldToolName("Trk::MagneticFieldTool"),
  m_magFieldToolInstanceName("ATLAS_TrackingMagFieldTool"),
  //m_layerArrayCreator(0),
  //m_layerArrayCreatorName("Trk::LayerArrayCreator"),
  //m_layerArrayCreatorInstanceName("LayerArrayCreator"),
  m_trackingVolumeArrayCreator(0),
  m_trackingVolumeArrayCreatorName("Trk::TrackingVolumeArrayCreator"),
  m_trackingVolumeArrayCreatorInstanceName("TrackingVolumeArrayCreator"),
  m_trackingVolumeHelper(0),
  m_trackingVolumeHelperName("Trk::TrackingVolumeHelper"),
  m_trackingVolumeHelperInstanceName("TrackingVolumeHelper"),
  m_trackingVolumeDisplayer(0),
  m_trackingVolumeDisplayerName("Trk::TrackingVolumeDisplayer"),
  m_trackingVolumeDisplayerInstanceName("TrackingVolumeDisplayer"),
  m_muonSimple(false),
  m_enclosingVolumes(true),
  m_innerBarrelRadius(4300.),
  m_outerBarrelRadius(13000.),
  m_barrelZ(6740.),
  m_innerEndcapZ(12900.),
  m_outerEndcapZ(22030.),
  m_beamPipeRadius(70.),
  m_outerEndcapRadius(12500.)
{
  declareInterface<Trk::IGeometryBuilder>(this);

  declareProperty("SimpleMuonGeometry",               m_muonSimple);  
  declareProperty("MagneticFieldTool",                m_magFieldToolName);  
  declareProperty("MagneticFieldTool",                m_magFieldToolName);  
  declareProperty("MagneticFieldToolInstance",        m_magFieldToolInstanceName);
  //declareProperty("LayerArrayCreator",                m_layerArrayCreatorName);
  //declareProperty("LayerArrayCreatorInstance",        m_layerArrayCreatorInstanceName); 
  // the overall dimensions
  declareProperty("InnerBarrelRadius",              m_innerBarrelRadius);
  declareProperty("OuterBarrelRadius",              m_outerBarrelRadius);
  declareProperty("BarrelZ",                        m_barrelZ);
  declareProperty("InnerEndcapZ",                   m_innerEndcapZ);
  declareProperty("OuterEndcapZ",                   m_outerEndcapZ);
}

// destructor
Muon::MuonTrackingGeometryBuilder::~MuonTrackingGeometryBuilder()
{}

// Athena standard methods
// initialize
StatusCode Muon::MuonTrackingGeometryBuilder::initialize()
{
    
    MsgStream log(msgSvc(), name());

    StatusCode s = AlgTool::initialize();
    
    s = toolSvc()->retrieveTool(m_magFieldToolName, m_magFieldToolInstanceName, m_magFieldTool);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_magFieldToolName << " from ToolSvc. MagneticField will be 0. " << endreq;
    }
    // Retrieve the tracking volume helper tool
    s = toolSvc()->retrieveTool(m_trackingVolumeHelperName, m_trackingVolumeHelperInstanceName, m_trackingVolumeHelper);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_trackingVolumeHelperName << " from ToolSvc. ";
      log <<" Creation of Gap Volumes will fail." << endreq;
    }
    // Retrieve the tracking volume helper tool
    s = toolSvc()->retrieveTool(m_trackingVolumeDisplayerName, m_trackingVolumeDisplayerInstanceName, m_trackingVolumeDisplayer);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_trackingVolumeDisplayerName << " from ToolSvc. ";
      log <<" Creation of Gap Volumes will fail." << endreq;
    }

    /*
    s = toolSvc()->retrieveTool(m_layerArrayCreatorName, m_layerArrayCreatorInstanceName, m_layerArrayCreator);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_layerArrayCreatorName << " from ToolSvc. ";
      log <<" Creation of LayerArrays might fail." << endreq;
    }
    */

    s = toolSvc()->retrieveTool(m_trackingVolumeArrayCreatorName, m_trackingVolumeArrayCreatorInstanceName, m_trackingVolumeArrayCreator);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_trackingVolumeArrayCreatorName << " from ToolSvc. ";
      log <<" Creation of LayerArrays might fail." << endreq;
    }
        
    log << MSG::INFO  << name() <<" initialize() successful" << endreq;    
    
  return StatusCode::SUCCESS;
}


const Trk::TrackingGeometry* Muon::MuonTrackingGeometryBuilder::trackingGeometry(const Trk::TrackingVolume* tvol) const
{

    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building tracking geometry" << endreq;    

// 0) Preparation //////////////////////////////////////////////////////////////////////////////////////
  // if no muon materials are declared, take default ones
  if (m_muonMaterialProperties.size() < 3){
    // set 0. / 0. / 0. / 0.
    m_muonMaterialProperties = std::vector<double>();
    m_muonMaterialProperties.push_back(0.);
    m_muonMaterialProperties.push_back(0.);
    m_muonMaterialProperties.push_back(0.);
  }

   m_muonMaterial = Trk::MaterialProperties(m_muonMaterialProperties[0],
                                            m_muonMaterialProperties[1],
                                            m_muonMaterialProperties[2]);

   Trk::MagneticFieldProperties muonMagneticFieldProperties(m_magFieldTool, Trk::RealisticField);    


   Trk::VolumeBounds* globalBounds = new Trk::CylinderVolumeBounds(m_outerBarrelRadius,
                                                                   m_outerEndcapZ);
                                                              
   if (m_muonSimple) {             
       Trk::TrackingVolume* topVolume = new Trk::TrackingVolume(new HepTransform3D,
                                                                globalBounds,
                                                                m_muonMaterial,
                                                                m_muonMagneticField,
                                                                0,0,
                                                                "GlobalVolume");
       m_standaloneTrackingVolume = topVolume;
       return new Trk::TrackingGeometry(topVolume);
   }     

    log << MSG::INFO  << name() <<"building barrel+innerEndcap+outerEndcap" << endreq;    

////////////////////////////////////////////////////////////////////////////////////////////////////////
// MuonSpectrometer contains:
//       - Barrel
//       - Endcaps inner/outer
     const Trk::TrackingVolume* muonBarrel = 0;
     Trk::CylinderVolumeBounds* barrelBounds = 0;
     const Trk::TrackingVolume* negativeMuonInnerEndcap = 0;
     Trk::CylinderVolumeBounds* negativeInnerEndcapBounds = 0;
     const Trk::TrackingVolume* negativeMuonOuterEndcap = 0;
     Trk::CylinderVolumeBounds* negativeOuterEndcapBounds = 0;
     const Trk::TrackingVolume* positiveMuonInnerEndcap = 0;
     Trk::CylinderVolumeBounds* positiveInnerEndcapBounds = 0;
     const Trk::TrackingVolume* positiveMuonOuterEndcap = 0;
     Trk::CylinderVolumeBounds* positiveOuterEndcapBounds = 0;
// volumes needed to close the geometry
     const Trk::TrackingVolume* negBeamPipe = 0;
     const Trk::TrackingVolume* posBeamPipe = 0;
     Trk::CylinderVolumeBounds* negBeamPipeBounds = 0;
     Trk::CylinderVolumeBounds* posBeamPipeBounds = 0;
     Trk::CylinderVolumeBounds* enclosedBounds = 0;
     const Trk::TrackingVolume* negCavern = 0;
     const Trk::TrackingVolume* posCavern = 0;
     Trk::CylinderVolumeBounds* negCavernBounds = 0;
     Trk::CylinderVolumeBounds* posCavernBounds = 0;
// enclosed volumes, if any
     std::vector<const Trk::TrackingVolume*> enclosedNegativeFaceVolumes;
     std::vector<const Trk::TrackingVolume*> enclosedCentralFaceVolumes;
     std::vector<const Trk::TrackingVolume*> enclosedPositiveFaceVolumes;
// if input, redefine dimensions
    if (tvol){
      // get dimensions
      const Trk::CylinderVolumeBounds* enclosedDetectorBounds 
            = dynamic_cast<const Trk::CylinderVolumeBounds*>(&(tvol->volumeBounds()));

      double enclosedDetectorHalflength = enclosedDetectorBounds->halflengthZ();
      double enclosedDetectorOuterRadius = enclosedDetectorBounds->outerRadius();
      // WHAT ABOUT BEAMPIPE ????????????
      log << MSG::INFO  << name() <<" dimensions of enclosed detectors (halfZ,outerR):"
          <<  enclosedDetectorHalflength<<","<< enclosedDetectorOuterRadius << endreq;    
      // check if input makes sense - gives warning if cuts into muon envelope
      if ( enclosedDetectorHalflength <= m_barrelZ ) {
        m_barrelZ =  enclosedDetectorHalflength ;
      } else {    
        log << MSG::WARNING  << name() <<" enclosed volume too long, cuts into muon envelope ?" << endreq;    
      }
      if ( enclosedDetectorOuterRadius <= m_innerBarrelRadius ) {
        m_innerBarrelRadius = enclosedDetectorOuterRadius ;
      } else {    
        log << MSG::WARNING  << name() <<" enclosed volume too wide, cuts into muon envelope ?" << endreq;    
      }
      // get the glue volumes
      const Trk::GlueVolumesDescriptor& enclosedDetGlueVolumes = tvol->glueVolumesDescriptor();
      // at negative face     
      enclosedNegativeFaceVolumes = enclosedDetGlueVolumes.glueVolumes(Trk::negativeFaceXY);
      // at cylinder cover     
      enclosedCentralFaceVolumes = enclosedDetGlueVolumes.glueVolumes(Trk::cylinderCover);
      // at positive face     
      enclosedPositiveFaceVolumes = enclosedDetGlueVolumes.glueVolumes(Trk::positiveFaceXY);
    }

// define all volumes
// muon barrel
   barrelBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
                                                m_outerBarrelRadius,
                                                m_barrelZ);
   muonBarrel = new Trk::TrackingVolume(0,
                                    barrelBounds,
                                    m_muonMaterial,
                                    m_muonMagneticField,
				    0,0,
                                    "Muon::Detectors::Barrel");
// inner Endcap
   double innerEndcapZHalfSize = 0.5*(m_innerEndcapZ - m_barrelZ);
   negativeInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerBarrelRadius,
                                                             innerEndcapZHalfSize);
   Hep3Vector negInnerEndcapPosition(0.,0.,-m_barrelZ-innerEndcapZHalfSize);
   //HepTransform3D* negInnerEndcapTransf = new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition);
   negativeMuonInnerEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition),
						 negativeInnerEndcapBounds,
                                                 m_muonMaterial,
                                                 m_muonMagneticField,
						 0,0,
						 "Muon::Detectors::NegativeInnerEndcap");
   positiveInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerBarrelRadius,
                                                             innerEndcapZHalfSize);
   Hep3Vector posInnerEndcapPosition(0.,0.,m_barrelZ+innerEndcapZHalfSize);
   positiveMuonInnerEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posInnerEndcapPosition),
						 positiveInnerEndcapBounds,
                                                 m_muonMaterial,
                                                 m_muonMagneticField,
						 0,0,
						 "Muon::Detectors::PositiveInnerEndcap");
// outer Endcap
   double outerEndcapZHalfSize = 0.5*(m_outerEndcapZ - m_innerEndcapZ);
   negativeOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerEndcapRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector negOuterEndcapPosition(0.,0.,-m_innerEndcapZ-outerEndcapZHalfSize);
   negativeMuonOuterEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition),
						     negativeOuterEndcapBounds,
                                                     m_muonMaterial,
                                                     m_muonMagneticField,
                                                     0,0,
						     "Muon::Detectors::NegativeOuterEndcap");
   positiveOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerEndcapRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector posOuterEndcapPosition(0.,0., m_innerEndcapZ+outerEndcapZHalfSize);
   positiveMuonOuterEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),
						     positiveOuterEndcapBounds, 
                                                     m_muonMaterial,
                                                     m_muonMagneticField,
                                                     0,0,
						     "Muon::Detectors::PositiveOuterEndcap");
// beamPipe
   negBeamPipeBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  outerEndcapZHalfSize+innerEndcapZHalfSize);
   posBeamPipeBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  outerEndcapZHalfSize+innerEndcapZHalfSize);
   Hep3Vector posBeamPipePosition(0.,0., m_outerEndcapZ-innerEndcapZHalfSize-outerEndcapZHalfSize);
   Hep3Vector negBeamPipePosition(0.,0.,-m_outerEndcapZ+innerEndcapZHalfSize+outerEndcapZHalfSize);
   negBeamPipe = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negBeamPipePosition),
                                      negBeamPipeBounds,
                                      m_muonMaterial,
                                      m_muonMagneticField,
                                      0,0,
                                      "Muons::Gaps::NegativeBeamPipe");
   posBeamPipe = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posBeamPipePosition),
                                      posBeamPipeBounds,
                                      m_muonMaterial,
                                      m_muonMagneticField,
                                      0,0,
                                      "Muons::Gaps::PositiveBeamPipe");
// enclosed (ID+Calo)
   if (!tvol) {
      enclosedBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                               m_innerBarrelRadius,
                                               m_barrelZ);
      tvol = new Trk::TrackingVolume(0,
                                   enclosedBounds,
                                   m_muonMaterial,
                                   m_muonMagneticField,
                                   0,0,
                                   "Muon::Detectors::EnclosedCentralDetectors");
   }
// cavern   
   negCavernBounds = new Trk::CylinderVolumeBounds(m_outerEndcapRadius,
                                                   m_outerBarrelRadius,
                                                   outerEndcapZHalfSize);
   posCavernBounds = new Trk::CylinderVolumeBounds(m_outerEndcapRadius,
                                                   m_outerBarrelRadius,
                                                   outerEndcapZHalfSize);

   negCavern = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition),
                                      negCavernBounds,
                                      m_muonMaterial,
                                      m_muonMagneticField,
				      0,0,
                                      "Muon::Gaps::NegativeCavern");
   posCavern = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),
                                      posCavernBounds,
                                      m_muonMaterial,
                                      m_muonMagneticField,
				       0,0,
                                      "Muon::Gaps::PositiveCavern");

    log << MSG::INFO  << name() <<" volumes defined " << endreq;    
//
// glue volumes at navigation level 
// radially , create corresponding volume arrays
    if ( enclosedCentralFaceVolumes.size()== 0 ) {
       m_trackingVolumeHelper->glueTrackingVolumes(*muonBarrel, Trk::tubeInnerCover,
                                                   *tvol, Trk::tubeOuterCover);
       enclosedNegativeFaceVolumes.push_back(tvol);
       enclosedNegativeFaceVolumes.push_back(muonBarrel);
       enclosedPositiveFaceVolumes.push_back(tvol);
       enclosedPositiveFaceVolumes.push_back(muonBarrel);
    } else {
       m_trackingVolumeHelper->glueTrackingVolumes(*muonBarrel, Trk::tubeInnerCover,
				                   enclosedCentralFaceVolumes, Trk::tubeOuterCover);
       enclosedNegativeFaceVolumes.push_back(muonBarrel);
       enclosedPositiveFaceVolumes.push_back(muonBarrel);
    }
   m_trackingVolumeHelper->glueTrackingVolumes(*negCavern, Trk::tubeInnerCover,
                                      *negativeMuonOuterEndcap, Trk::tubeOuterCover);
   m_trackingVolumeHelper->glueTrackingVolumes(*posCavern, Trk::tubeInnerCover,
                                      *positiveMuonOuterEndcap, Trk::tubeOuterCover);
    log << MSG::INFO  << name() <<" first radial gluing ok " << endreq;    
/// create corresponding volumeArrays  
// barrel
   std::vector<const Trk::TrackingVolume*> barrelNavVolumes;
   barrelNavVolumes.push_back(tvol);
   barrelNavVolumes.push_back(muonBarrel);
   Trk::BinnedArray<Trk::TrackingVolume>* barrelNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(barrelNavVolumes, true);
   Trk::BinnedArray<Trk::TrackingVolume>* barrelNegNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(enclosedNegativeFaceVolumes, true);
   Trk::BinnedArray<Trk::TrackingVolume>* barrelPosNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(enclosedPositiveFaceVolumes, true);
// negative outer endcap
   std::vector<const Trk::TrackingVolume*> negOuterEndcapNavVolumes;
   negOuterEndcapNavVolumes.push_back(negativeMuonOuterEndcap);
   negOuterEndcapNavVolumes.push_back(negCavern);
   Trk::BinnedArray<Trk::TrackingVolume>* negOuterEndcapNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(negOuterEndcapNavVolumes, true);
// positive outer endcap
   std::vector<const Trk::TrackingVolume*> posOuterEndcapNavVolumes;
   posOuterEndcapNavVolumes.push_back(positiveMuonOuterEndcap);
   posOuterEndcapNavVolumes.push_back(posCavern);
   Trk::BinnedArray<Trk::TrackingVolume>* posOuterEndcapNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(posOuterEndcapNavVolumes, true);
// glue volumes/volume arrays in Z
// endcaps
   m_trackingVolumeHelper->setOutsideTrackingVolume(*negativeMuonOuterEndcap,Trk::positiveFaceXY,negativeMuonInnerEndcap);
   m_trackingVolumeHelper->setOutsideTrackingVolume(*negCavern,Trk::positiveFaceXY,negativeMuonInnerEndcap);
   m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*negativeMuonInnerEndcap,Trk::negativeFaceXY,negOuterEndcapNavArray);
   m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*positiveMuonInnerEndcap,Trk::positiveFaceXY,posOuterEndcapNavArray);
   m_trackingVolumeHelper->setOutsideTrackingVolume(*positiveMuonOuterEndcap,Trk::negativeFaceXY,positiveMuonInnerEndcap);
   m_trackingVolumeHelper->setOutsideTrackingVolume(*posCavern,Trk::negativeFaceXY,positiveMuonInnerEndcap);
// barrel glued to endcap
   for (unsigned int i=0; i<enclosedNegativeFaceVolumes.size();i++){ 
     m_trackingVolumeHelper->setOutsideTrackingVolume(*(enclosedNegativeFaceVolumes[i]),Trk::negativeFaceXY,negativeMuonInnerEndcap);
   }
   for (unsigned int i=0; i<enclosedPositiveFaceVolumes.size();i++){ 
     m_trackingVolumeHelper->setOutsideTrackingVolume(*(enclosedPositiveFaceVolumes[i]),Trk::positiveFaceXY,positiveMuonInnerEndcap);
   }
   m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*negativeMuonInnerEndcap,Trk::positiveFaceXY,barrelNegNavArray);
   m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*positiveMuonInnerEndcap,Trk::negativeFaceXY,barrelPosNavArray);
// glue beam pipe - > at navigation level in Z
   m_trackingVolumeHelper->setOutsideTrackingVolume(*negBeamPipe,Trk::positiveFaceXY,enclosedNegativeFaceVolumes[0]);
   m_trackingVolumeHelper->setOutsideTrackingVolume(*posBeamPipe,Trk::negativeFaceXY,enclosedPositiveFaceVolumes[0]);
// glue beam pipe - > at navigation level radially
   std::vector<const Trk::TrackingVolume*> negBPVol;
   negBPVol.push_back(negativeMuonOuterEndcap);
   negBPVol.push_back(negativeMuonInnerEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* negBPGlueArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(negBPVol, true);
   m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*negBeamPipe, Trk::cylinderCover,
                                                 negBPGlueArray);
   m_trackingVolumeHelper->setInsideTrackingVolume(*negativeMuonOuterEndcap, Trk::tubeInnerCover,negBeamPipe);
   m_trackingVolumeHelper->setInsideTrackingVolume(*negativeMuonInnerEndcap, Trk::tubeInnerCover,negBeamPipe);
   std::vector<const Trk::TrackingVolume*> posBPVol;
   posBPVol.push_back(positiveMuonInnerEndcap);
   posBPVol.push_back(positiveMuonOuterEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* posBPGlueArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(posBPVol, true);
   m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*posBeamPipe, Trk::cylinderCover,
                                                 posBPGlueArray);
   m_trackingVolumeHelper->setInsideTrackingVolume(*positiveMuonOuterEndcap, Trk::tubeInnerCover,posBeamPipe);
   m_trackingVolumeHelper->setInsideTrackingVolume(*positiveMuonInnerEndcap, Trk::tubeInnerCover,posBeamPipe);
   log << MSG::INFO  << name() <<"volumes glued to beam pipe " << endreq;    
// define enveloping volumes
   barrelBounds = new Trk::CylinderVolumeBounds(0.,                       // central beam pipe belongs to enclosed 
                                                m_outerBarrelRadius,
                                                m_barrelZ);
   Trk::TrackingVolume* barrel  =  new Trk::TrackingVolume(0,
                                   barrelBounds,
                                   m_muonMaterial,
                                   m_muonMagneticField,
                                   0,barrelNavArray,
                                   "All::Container::Barrel");

   Trk::CylinderVolumeBounds* negOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                        m_outerBarrelRadius,
                                                        outerEndcapZHalfSize);
   Trk::TrackingVolume* negOuterEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition),
                                           negOuterEndcapBounds,
                                           m_muonMaterial,
                                           m_muonMagneticField,
                                           0,negOuterEndcapNavArray,
                                           "Muon::Container::NegativeOuterEndcap");
   Trk::CylinderVolumeBounds* posOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                        m_outerBarrelRadius,
                                                        outerEndcapZHalfSize);
   Trk::TrackingVolume* posOuterEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),
                                           posOuterEndcapBounds,
                                           m_muonMaterial,
                                           m_muonMagneticField,
                                           0,posOuterEndcapNavArray, 
                                           "Muon::Container::PositiveOuterEndcap");
// create volume array and envelope
   std::vector<const Trk::TrackingVolume*> endcapNegNavVolumes;
   endcapNegNavVolumes.push_back(negOuterEndcap);
   endcapNegNavVolumes.push_back(negativeMuonInnerEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* endcapNegNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(endcapNegNavVolumes, true);
   Trk::CylinderVolumeBounds* negNavEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  m_outerBarrelRadius,
                                                  negBeamPipeBounds->halflengthZ());
   Trk::TrackingVolume* negNavEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negBeamPipePosition),
                                                             negNavEndcapBounds,
                                                             m_muonMaterial,
                                                             m_muonMagneticField,
                                                             0,endcapNegNavArray,
                                                             "Muon::Container::NegativeEndcap");
// create volume array and envelope
   std::vector<const Trk::TrackingVolume*> endcapPosNavVolumes;
   endcapPosNavVolumes.push_back(positiveMuonInnerEndcap);
   endcapPosNavVolumes.push_back(posOuterEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* endcapPosNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(endcapPosNavVolumes, true);
   Trk::CylinderVolumeBounds* posNavEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  m_outerBarrelRadius,
                                                  posBeamPipeBounds->halflengthZ());
   Trk::TrackingVolume* posNavEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posBeamPipePosition),
                                                             posNavEndcapBounds,
                                                             m_muonMaterial,
                                                             m_muonMagneticField,
                                                             0,endcapPosNavArray,
                                                             "Muon::Container::PositiveEndcap");
// add beam pipe to complete negative endcap 
   std::vector<const Trk::TrackingVolume*> endcapNegVolumes;
   endcapNegVolumes.push_back(negBeamPipe);
   endcapNegVolumes.push_back(negNavEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* endcapNegArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(endcapNegVolumes, true);
   Trk::CylinderVolumeBounds* negEndcapBounds = new Trk::CylinderVolumeBounds( m_outerBarrelRadius,
                                                                               negBeamPipeBounds->halflengthZ());
   Trk::TrackingVolume* negEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negBeamPipePosition),
                                                             negEndcapBounds,
                                                             m_muonMaterial,
                                                             m_muonMagneticField,
                                                             0,endcapNegArray,
                                                             "All::Container::NegativeEndcap");
// add beam pipe to complete positive endcap
   std::vector<const Trk::TrackingVolume*> endcapPosVolumes;
   endcapPosVolumes.push_back(posBeamPipe);
   endcapPosVolumes.push_back(posNavEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* endcapPosArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(endcapPosVolumes, true);
   Trk::CylinderVolumeBounds* posEndcapBounds = new Trk::CylinderVolumeBounds( m_outerBarrelRadius,
                                                                               posBeamPipeBounds->halflengthZ());
   Trk::TrackingVolume* posEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posBeamPipePosition),
                                                             posEndcapBounds,
                                                             m_muonMaterial,
                                                             m_muonMagneticField,
                                                             0,endcapPosArray,
                                                             "All::Container::PositiveEndcap");
//  
   std::vector<const Trk::TrackingVolume*> detectorVolumes;
   detectorVolumes.push_back(negEndcap);
   detectorVolumes.push_back(barrel);
   detectorVolumes.push_back(posEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* detectorArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(detectorVolumes, true);
   Trk::CylinderVolumeBounds* detectorBounds = new Trk::CylinderVolumeBounds( m_outerBarrelRadius,
                                                                              m_outerEndcapZ);
   Trk::TrackingVolume* detector  =  new Trk::TrackingVolume(0,
                                                             detectorBounds,
                                                             m_muonMaterial,
                                                             m_muonMagneticField,
                                                             0,detectorArray,
                                                             "All::Container::CompleteDetector");
    log << MSG::INFO  << name() <<" detector envelope created " << endreq;    
//
   Trk::TrackingGeometry* trackingGeometry = new Trk::TrackingGeometry(detector);

   log << MSG::INFO  << name() <<" returning tracking geometry " << endreq;    
   return trackingGeometry;  
}

// finalize
StatusCode Muon::MuonTrackingGeometryBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    return StatusCode::SUCCESS;
}
