///////////////////////////////////////////////////////////////////
// MuonTrackingGeometryBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonTrackingGeometryBuilder.h"
// Trk
//#include "TrkDetDescrInterfaces/ITrackingVolumeBuilder.h"
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkDetDescrUtils/BinUtility1DR.h"
#include "TrkDetDescrUtils/BinUtility1DZ.h"
#include "TrkDetDescrUtils/BinnedArray.h"
//#include "TrkDetDescrUtils/BinningType.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkGeometry/CylinderLayer.h"
#include "TrkGeometry/DiscLayer.h"
#include "TrkGeometry/PlaneLayer.h"
//#include "TrkGeometry/SimplifiedMaterialProperties.h"
#include "TrkMagFieldTools/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/TrackingGeometry.h"
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
  m_layerArrayCreator(0),
  m_layerArrayCreatorName("Trk::LayerArrayCreator"),
  m_layerArrayCreatorInstanceName("LayerArrayCreator"),
  m_trackingVolumeArrayCreator(0),
  m_trackingVolumeArrayCreatorName("Trk::TrackingVolumeArrayCreator"),
  m_trackingVolumeArrayCreatorInstanceName("TrackingVolumeArrayCreator"),
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
  declareProperty("LayerArrayCreator",                m_layerArrayCreatorName);
  declareProperty("LayerArrayCreatorInstance",        m_layerArrayCreatorInstanceName); 
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

    s = toolSvc()->retrieveTool(m_layerArrayCreatorName, m_layerArrayCreatorInstanceName, m_layerArrayCreator);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_layerArrayCreatorName << " from ToolSvc. ";
      log <<" Creation of LayerArrays might fail." << endreq;
    }

    s = toolSvc()->retrieveTool(m_trackingVolumeArrayCreatorName, m_trackingVolumeArrayCreatorInstanceName, m_trackingVolumeArrayCreator);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_trackingVolumeArrayCreatorName << " from ToolSvc. ";
      log <<" Creation of LayerArrays might fail." << endreq;
    }
        
    log << MSG::INFO  << name() <<" initialize() successful" << endreq;    
    
  return StatusCode::SUCCESS;
}


const Trk::TrackingGeometry* Muon::MuonTrackingGeometryBuilder::trackingGeometry(const Trk::TrackingVolume* tvol, std::vector<const Trk::TrackingVolume*>* gluevol) const
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

    log << MSG::INFO  << name() <<" barrel+innerEndcap+outerEndcap" << endreq;    
// if input, redefine dimensions

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
     const Trk::TrackingVolume* beamPipe = 0;
     Trk::CylinderVolumeBounds* beamPipeBounds = 0;
     const Trk::TrackingVolume* other = 0;
     Trk::CylinderVolumeBounds* otherBounds = 0;
     const Trk::TrackingVolume* negCavern = 0;
     const Trk::TrackingVolume* posCavern = 0;
     Trk::CylinderVolumeBounds* negCavernBounds = 0;
     Trk::CylinderVolumeBounds* posCavernBounds = 0;
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
                                    "MuonBarrel");
// inner Endcap
   double innerEndcapZSize = m_innerEndcapZ - m_barrelZ;
   negativeInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerBarrelRadius,
                                                             0.5*innerEndcapZSize);
   Hep3Vector negInnerEndcapPosition(0.,0.,-m_barrelZ-0.5*innerEndcapZSize);
   //HepTransform3D* negInnerEndcapTransf = new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition);
   negativeMuonInnerEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition),
						 negativeInnerEndcapBounds,
                                                 m_muonMaterial,
                                                 m_muonMagneticField,
						 0,0,
						 "MuonNegativeInnerEndcap");
   positiveInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerBarrelRadius,
                                                             0.5*innerEndcapZSize);
   Hep3Vector posInnerEndcapPosition(0.,0.,m_barrelZ+0.5*innerEndcapZSize);
   positiveMuonInnerEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posInnerEndcapPosition),
						 positiveInnerEndcapBounds,
                                                 m_muonMaterial,
                                                 m_muonMagneticField,
						 0,0,
						 "MuonPositiveInnerEndcap");
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
						     "MuonNegativeOuterEndcap");
   positiveOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerEndcapRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector posOuterEndcapPosition(0.,0., m_innerEndcapZ+outerEndcapZHalfSize);
   positiveMuonOuterEndcap = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),
						     positiveOuterEndcapBounds, 
                                                     m_muonMaterial,
                                                     m_muonMagneticField,
                                                     0,0,
						     "MuonPositiveOuterEndcap");
// beamPipe
   beamPipeBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  m_outerEndcapZ);
   beamPipe = new Trk::TrackingVolume(0,
                                      beamPipeBounds,
                                      m_muonMaterial,
                                      m_muonMagneticField,
                                      0,0,
                                      "BeamPipe");
// other (ID+Calo)
   otherBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                               m_innerBarrelRadius,
                                               m_barrelZ);
   other = new Trk::TrackingVolume(0,
                                   otherBounds,
                                   m_muonMaterial,
                                   m_muonMagneticField,
                                   0,0,
                                   "Other");
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
                                      "NegativeCavern");
   posCavern = new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),
                                      posCavernBounds,
                                      m_muonMaterial,
                                      m_muonMagneticField,
				       0,0,
                                      "PositiveCavern");

    log << MSG::INFO  << name() <<" volumes defined " << endreq;    
// glue volumes  
// radially
   Trk::IGeometryBuilder::glueVolumes(*muonBarrel, Trk::tubeInnerCover,
                                      *other, Trk::tubeOuterCover);
   Trk::IGeometryBuilder::glueVolumes(*negCavern, Trk::tubeInnerCover,
                                      *negativeMuonOuterEndcap, Trk::tubeOuterCover);
   Trk::IGeometryBuilder::glueVolumes(*posCavern, Trk::tubeInnerCover,
                                      *positiveMuonOuterEndcap, Trk::tubeOuterCover);
    log << MSG::INFO  << name() <<" first radial gluing ok " << endreq;    
/// create corresponding volumeArrays  
// barrel
   std::vector<const Trk::TrackingVolume*> barrelNavVolumes;
   barrelNavVolumes.push_back(other);
   barrelNavVolumes.push_back(muonBarrel);
   Trk::BinnedArray<Trk::TrackingVolume>* barrelNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(barrelNavVolumes, true);
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
   Trk::IGeometryBuilder::setOutsideVolume(*negativeMuonOuterEndcap,Trk::positiveFaceXY,negativeMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolume(*negCavern,Trk::positiveFaceXY,negativeMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolumeArray(*negativeMuonInnerEndcap,Trk::negativeFaceXY,negOuterEndcapNavArray);
   Trk::IGeometryBuilder::setOutsideVolumeArray(*positiveMuonInnerEndcap,Trk::positiveFaceXY,posOuterEndcapNavArray);
   Trk::IGeometryBuilder::setOutsideVolume(*positiveMuonOuterEndcap,Trk::negativeFaceXY,positiveMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolume(*posCavern,Trk::negativeFaceXY,positiveMuonInnerEndcap);
// barrel glued to endcap
   Trk::IGeometryBuilder::setOutsideVolume(*other,Trk::negativeFaceXY,negativeMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolume(*other,Trk::positiveFaceXY,positiveMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolume(*muonBarrel,Trk::negativeFaceXY,negativeMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolume(*muonBarrel,Trk::positiveFaceXY,positiveMuonInnerEndcap);
   Trk::IGeometryBuilder::setOutsideVolumeArray(*negativeMuonInnerEndcap,Trk::positiveFaceXY,barrelNavArray);
   Trk::IGeometryBuilder::setOutsideVolumeArray(*positiveMuonInnerEndcap,Trk::negativeFaceXY,barrelNavArray);
// define enveloping volumes
   barrelBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                m_outerBarrelRadius,
                                                m_barrelZ);
   Trk::TrackingVolume* barrel  =  new Trk::TrackingVolume(0,
                                   barrelBounds,
                                   m_muonMaterial,
                                   m_muonMagneticField,
                                   0,barrelNavArray,
                                   "Barrel");

   Trk::CylinderVolumeBounds* negOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                        m_outerBarrelRadius,
                                                        outerEndcapZHalfSize);
   Trk::TrackingVolume* negOuterEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition),
                                           negOuterEndcapBounds,
                                           m_muonMaterial,
                                           m_muonMagneticField,
                                           0,negOuterEndcapNavArray,
                                           "NegativeOuterEndcap");
   Trk::CylinderVolumeBounds* posOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                        m_outerBarrelRadius,
                                                        outerEndcapZHalfSize);
   Trk::TrackingVolume* posOuterEndcap  =  new Trk::TrackingVolume(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),
                                           posOuterEndcapBounds,
                                           m_muonMaterial,
                                           m_muonMagneticField,
                                           0,posOuterEndcapNavArray, 
                                           "PositiveOuterEndcap");
    log << MSG::INFO  << name() <<" enveloping volumes ok " << endreq;    
// create volume array and envelope
   std::vector<const Trk::TrackingVolume*> detectorNavVolumes;
   detectorNavVolumes.push_back(negOuterEndcap);
   detectorNavVolumes.push_back(negativeMuonInnerEndcap);
   detectorNavVolumes.push_back(barrel);
   detectorNavVolumes.push_back(positiveMuonInnerEndcap);
   detectorNavVolumes.push_back(posOuterEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* detectorNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(detectorNavVolumes, true);
   Trk::CylinderVolumeBounds* detectorBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  m_outerBarrelRadius,
                                                  m_outerEndcapZ);
   Trk::TrackingVolume* detector  =  new Trk::TrackingVolume(0,
                                                             detectorBounds,
                                                             m_muonMaterial,
                                                             m_muonMagneticField,
                                                             0,detectorNavArray,
                                                             "Detector");
    log << MSG::INFO  << name() <<" detector envelope created " << endreq;    
// finally, glue detectors with beam pipe
   std::vector<const Trk::TrackingVolume*> detectorGlueVolumes;
   detectorGlueVolumes.push_back(negativeMuonOuterEndcap);
   detectorGlueVolumes.push_back(negativeMuonInnerEndcap);
   detectorGlueVolumes.push_back(other);
   detectorGlueVolumes.push_back(positiveMuonInnerEndcap);
   detectorGlueVolumes.push_back(positiveMuonOuterEndcap);
   Trk::BinnedArray<Trk::TrackingVolume>* detectorGlueArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInZ(detectorGlueVolumes, true);
   Trk::IGeometryBuilder::setOutsideVolumeArray(*beamPipe, Trk::tubeOuterCover,
                                                 detectorGlueArray);
   Trk::IGeometryBuilder::setInsideVolume(*negativeMuonOuterEndcap, Trk::tubeInnerCover,beamPipe);
   Trk::IGeometryBuilder::setInsideVolume(*negativeMuonInnerEndcap, Trk::tubeInnerCover,beamPipe);
   Trk::IGeometryBuilder::setInsideVolume(*other, Trk::tubeInnerCover,beamPipe);
   Trk::IGeometryBuilder::setInsideVolume(*positiveMuonInnerEndcap, Trk::tubeInnerCover,beamPipe);
   Trk::IGeometryBuilder::setInsideVolume(*positiveMuonOuterEndcap, Trk::tubeInnerCover,beamPipe);
   log << MSG::INFO  << name() <<"volumes glued to beam pipe " << endreq;    
//
   std::vector<const Trk::TrackingVolume*>  globalNavVolumes;
   globalNavVolumes.push_back(beamPipe);
   globalNavVolumes.push_back(detector);
   Trk::BinnedArray<Trk::TrackingVolume>* globalNavArray =
   m_trackingVolumeArrayCreator->cylinderVolumesArrayInR(globalNavVolumes, true);
// close geometry ( global boundaries already defined )
   Trk::TrackingVolume* topVolume = new Trk::TrackingVolume(0,
                                                            globalBounds,
                                                            m_muonMaterial,
                                                            m_muonMagneticField,
                                                            0,globalNavArray,
                                                            "GlobalVolume");

   Trk::TrackingGeometry* trackingGeometry = new Trk::TrackingGeometry(topVolume);

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
