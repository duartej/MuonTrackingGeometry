///////////////////////////////////////////////////////////////////
// MuonTrackingGeometryBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonTrackingGeometryBuilder.h"
#include "MuonGeoModel/GlobalUtilities.h" 
// Trk
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeHelper.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeDisplayer.h"
#include "TrkDetDescrUtils/BinUtility1DR.h"
#include "TrkDetDescrUtils/BinUtility1DZ.h"
#include "TrkDetDescrUtils/BinUtility2DPhiZ.h"
#include "TrkDetDescrUtils/BinnedArray.h"
#include "TrkDetDescrUtils/BinnedArray2D.h"
//#include "TrkDetDescrUtils/BinningType.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkDetDescrTools/TrackingVolumeArrayCreator.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/DoubleTrapezoidVolumeBounds.h"
#include "TrkVolumes/BevelledCylinderVolumeBounds.h"
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
  m_stationBuilder(0),
  m_stationBuilderName("Muon::MuonStationBuilder"),
  m_stationBuilderInstanceName("MuonStationBuilder"),
  m_inertBuilder(0),
  m_inertBuilderName("Muon::MuonInertMaterialBuilder"),
  m_inertBuilderInstanceName("MuonInertMaterialBuilder"),
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
  m_muonActive(false),
  m_muonInert(false),
  m_innerBarrelRadius(4300.),
  m_outerBarrelRadius(13000.),
  m_barrelZ(6740.),
  m_innerEndcapZ(12900.),
  //m_outerEndcapZ(22030.),
  m_outerEndcapZ(25000.),
  m_beamPipeRadius(70.),
  m_outerEndcapRadius(12500.),
  m_barrelEtaPartition(1),
  m_innerEndcapEtaPartition(1),
  m_outerEndcapEtaPartition(1),
  m_phiPartition(1)
{
  declareInterface<Trk::IGeometryBuilder>(this);

  declareProperty("SimpleMuonGeometry",               m_muonSimple);  
  declareProperty("BuildActiveMaterial",              m_muonActive);  
  declareProperty("BuildInertMaterial",               m_muonInert);
  declareProperty("MagneticFieldTool",                m_magFieldToolName);  
  declareProperty("MagneticFieldTool",                m_magFieldToolName);  
  declareProperty("MagneticFieldToolInstance",        m_magFieldToolInstanceName);
  // the overall dimensions
  declareProperty("InnerBarrelRadius",              m_innerBarrelRadius);
  declareProperty("OuterBarrelRadius",              m_outerBarrelRadius);
  declareProperty("BarrelZ",                        m_barrelZ);
  declareProperty("InnerEndcapZ",                   m_innerEndcapZ);
  declareProperty("OuterEndcapZ",                   m_outerEndcapZ);
  declareProperty("EtaBarrelPartitions",            m_barrelEtaPartition);
  declareProperty("EtaInnerEndcapPartitions",       m_innerEndcapEtaPartition);
  declareProperty("EtaOuterEndcapPartitions",       m_outerEndcapEtaPartition);
  declareProperty("PhiPartitions",                  m_phiPartition);
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

    if (m_muonActive) { 
      s = toolSvc()->retrieveTool(m_stationBuilderName, m_stationBuilderInstanceName, m_stationBuilder);
      if (s.isFailure())
      {
        log << MSG::ERROR << "Could not retrieve " << m_stationBuilderName << " from ToolSvc. ";
        log <<" Creation of stations might fail." << endreq;
      }
    }

    if (m_muonInert) {
      s = toolSvc()->retrieveTool(m_inertBuilderName, m_inertBuilderInstanceName, m_inertBuilder);
      if (s.isFailure())
      {
        log << MSG::ERROR << "Could not retrieve " << m_inertBuilderName << " from ToolSvc. ";
        log <<" Creation of inert material objects might fail." << endreq;
      }
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


const Trk::TrackingGeometry* Muon::MuonTrackingGeometryBuilder::trackingGeometry(const Trk::TrackingVolume* tvol) const
{

    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building tracking geometry" << endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////
    // process muon material objects
    m_stations =0;
    m_inertObjs=0;
    if (m_muonActive && m_stationBuilder) m_stations = m_stationBuilder->buildDetachedTrackingVolumes();
    if (m_muonInert && m_inertBuilder) m_inertObjs = m_inertBuilder->buildDetachedTrackingVolumes();

    // find object's span with tolerance for the alignment ( 50 cm in z, 5 degrees in phi - starting guess)
    // const std::vector<const Span*>* stationsSpan = findVolumeSpan(stations,50.*cm,5*deg);
    m_stationSpan = findVolumeSpan(m_stations, 0.,0.);
    m_inertSpan = findVolumeSpan(m_inertObjs,0.,0.);

    /*
    for (unsigned int i=0; i<m_stations->size();i++) {
      std::cout << "span of chamber:" << (*m_stations)[i]->name() << ":" << 
	(*((*m_stationSpan)[i]))[0] << "," <<	(*((*m_stationSpan)[i]))[1] <<"," <<
	(*((*m_stationSpan)[i]))[2] << "," <<	(*((*m_stationSpan)[i]))[3] << std::endl;   
    }
    */
// 0) Preparation //////////////////////////////////////////////////////////////////////////////////////
  // if no muon materials are declared, take default ones
  if (m_muonMaterialProperties.size() < 3){
    // set 0. / 0. / 0. / 0.
    m_muonMaterialProperties = std::vector<double>();
    m_muonMaterialProperties.push_back(0.);
    m_muonMaterialProperties.push_back(10e10);
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

// define basic volumes
// muon barrel
   barrelBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
                                                m_outerBarrelRadius,
                                                m_barrelZ);
   const Trk::Volume barrelVol(new HepTransform3D(),barrelBounds);
// process volume
   muonBarrel = processVolume( &barrelVol,m_barrelEtaPartition,m_phiPartition,"Muon::Detectors::Barrel"); 
// inner Endcap
   double innerEndcapZHalfSize = 0.5*(m_innerEndcapZ - m_barrelZ);
   negativeInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerBarrelRadius,
                                                             innerEndcapZHalfSize);
   Hep3Vector negInnerEndcapPosition(0.,0.,-m_barrelZ-innerEndcapZHalfSize);
   HepTransform3D* negInnerEndcapTransf = new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition);
   const Trk::Volume negIECVol(negInnerEndcapTransf,negativeInnerEndcapBounds);
   negativeMuonInnerEndcap = processVolume( &negIECVol,m_innerEndcapEtaPartition,m_phiPartition,
					    "Muon::Detectors::NegativeInnerEndcap" ); 
//
   positiveInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerBarrelRadius,
                                                             innerEndcapZHalfSize);
   Hep3Vector posInnerEndcapPosition(0.,0.,m_barrelZ+innerEndcapZHalfSize);
   HepTransform3D* posInnerEndcapTransf = new HepTransform3D(Trk::s_idRotation,posInnerEndcapPosition);
   const Trk::Volume posIECVol(posInnerEndcapTransf,positiveInnerEndcapBounds);
   positiveMuonInnerEndcap = processVolume( &posIECVol,m_innerEndcapEtaPartition,m_phiPartition,
					    "Muon::Detectors::PositiveInnerEndcap" ); 
// outer Endcap
   double outerEndcapZHalfSize = 0.5*(m_outerEndcapZ - m_innerEndcapZ);
   negativeOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerEndcapRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector negOuterEndcapPosition(0.,0.,-m_innerEndcapZ-outerEndcapZHalfSize);
   HepTransform3D* negOuterEndcapTransf = new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition);
   const Trk::Volume negOECVol(negOuterEndcapTransf,negativeOuterEndcapBounds);
   negativeMuonOuterEndcap = processVolume( &negOECVol,m_outerEndcapEtaPartition,m_phiPartition,
					    "Muon::Detectors::NegativeOuterEndcap" ); 
//
   positiveOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                             m_outerEndcapRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector posOuterEndcapPosition(0.,0., m_innerEndcapZ+outerEndcapZHalfSize);
   HepTransform3D* posOuterEndcapTransf = new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition);
   const Trk::Volume posOECVol(posOuterEndcapTransf,positiveOuterEndcapBounds);
   positiveMuonOuterEndcap = processVolume( &posOECVol,m_outerEndcapEtaPartition,m_phiPartition,
					    "Muon::Detectors::PositiveOuterEndcap" ); 
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
      enclosedBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
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
                                                   *tvol, Trk::cylinderCover);
       enclosedNegativeFaceVolumes.push_back(tvol);
       enclosedNegativeFaceVolumes.push_back(muonBarrel);
       enclosedPositiveFaceVolumes.push_back(tvol);
       enclosedPositiveFaceVolumes.push_back(muonBarrel);
    } else {
       m_trackingVolumeHelper->glueTrackingVolumes(*muonBarrel, Trk::tubeInnerCover,
				                   enclosedCentralFaceVolumes, Trk::cylinderCover);
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

const std::vector<const Muon::Span*>* Muon::MuonTrackingGeometryBuilder::findVolumeSpan(const std::vector<const Trk::DetachedTrackingVolume*>*& objs, double zTol, double phiTol) const
{
  MsgStream log(msgSvc(), name());

  if (!objs) return 0;
  std::vector<const Span*> volSpans;
  for (unsigned int iobj=0; iobj<objs->size(); iobj++) {
    // get object boundaries, transform
    // Trk::VolumeBounds  bounds = (*objs)[iobj]->trackingVolume()->volumeBounds();    
    HepTransform3D  transform = (*objs)[iobj]->trackingVolume()->transform();
    // volume shape
    const Trk::CuboidVolumeBounds* box = dynamic_cast<const Trk::CuboidVolumeBounds*>
      (&((*objs)[iobj]->trackingVolume()->volumeBounds()));
    const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*>
      (&((*objs)[iobj]->trackingVolume()->volumeBounds()));
    const Trk::DoubleTrapezoidVolumeBounds* dtrd = dynamic_cast<const Trk::DoubleTrapezoidVolumeBounds*>
      (&((*objs)[iobj]->trackingVolume()->volumeBounds()));
    const Trk::BevelledCylinderVolumeBounds* bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*>
      (&((*objs)[iobj]->trackingVolume()->volumeBounds()));
    const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*>
      (&((*objs)[iobj]->trackingVolume()->volumeBounds()));
    // loop over edges ...
    double minZ = m_outerEndcapZ ; double maxZ = - m_outerEndcapZ;
    double minPhi = 2*M_PI; double maxPhi = 0.;
    std::vector<Trk::GlobalPosition> edges;
    Muon::Span span;  
    if (box) {
      edges.push_back( Trk::GlobalPosition(box->halflengthX(),box->halflengthY(),box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-box->halflengthX(),box->halflengthY(),box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(box->halflengthX(),-box->halflengthY(),box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-box->halflengthX(),-box->halflengthY(),box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(box->halflengthX(),box->halflengthY(),-box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-box->halflengthX(),box->halflengthY(),-box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(box->halflengthX(),-box->halflengthY(),-box->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-box->halflengthX(),-box->halflengthY(),-box->halflengthZ()) );
    }
    if (trd) {
      edges.push_back( Trk::GlobalPosition(trd->maxHalflengthX(),trd->halflengthY(),trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-trd->maxHalflengthX(),trd->halflengthY(),trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(trd->minHalflengthX(),-trd->halflengthY(),trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-trd->minHalflengthX(),-trd->halflengthY(),trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(trd->maxHalflengthX(),trd->halflengthY(),-trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-trd->maxHalflengthX(),trd->halflengthY(),-trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(trd->minHalflengthX(),-trd->halflengthY(),-trd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-trd->minHalflengthX(),-trd->halflengthY(),-trd->halflengthZ()) );
    }
    if (dtrd) {
      edges.push_back( Trk::GlobalPosition( dtrd->maxHalflengthX(),2*dtrd->halflengthY2(),dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-dtrd->maxHalflengthX(),2*dtrd->halflengthY2(),dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition( dtrd->medHalflengthX(),0.,dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-dtrd->medHalflengthX(),0.,dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition( dtrd->minHalflengthX(),-2*dtrd->halflengthY1(),dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-dtrd->minHalflengthX(),-2*dtrd->halflengthY1(),dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition( dtrd->maxHalflengthX(),2*dtrd->halflengthY2(),-dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-dtrd->maxHalflengthX(),2*dtrd->halflengthY2(),-dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition( dtrd->medHalflengthX(),0.,-dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-dtrd->medHalflengthX(),0.,-dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition( dtrd->minHalflengthX(),-2*dtrd->halflengthY1(),-dtrd->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(-dtrd->minHalflengthX(),-2*dtrd->halflengthY1(),-dtrd->halflengthZ()) );
    }
    if (bcyl) {
      edges.push_back( Trk::GlobalPosition(0.,0., bcyl->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(0.,0.,-bcyl->halflengthZ()) );
    }
    if (cyl) {
      edges.push_back( Trk::GlobalPosition(0.,0., cyl->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(0.,0.,-cyl->halflengthZ()) );
    }
    // apply transform and get span
    for (unsigned int ie=0; ie < edges.size() ; ie++) {
      Trk::GlobalPosition gp = transform*edges[ie];
      double phi = gp.phi()+M_PI; 
      if ( gp[2]<minZ ) minZ = gp[2];
      if ( gp[2]>maxZ ) maxZ = gp[2];
      if ( phi < minPhi ) minPhi = phi; 
      if ( phi > maxPhi ) maxPhi = phi; 
    }
    if ( maxPhi-minPhi > M_PI ) {
      double medPhi = minPhi;
      minPhi = maxPhi;
      maxPhi = medPhi;
    }
    if ( box || trd || dtrd ) {
      span.push_back( minZ - zTol );  
      span.push_back( maxZ + zTol );  
      span.push_back( minPhi - phiTol );  
      span.push_back( maxPhi + phiTol );  
    } else if (bcyl) {
      span.push_back( minZ - bcyl->outerRadius()-zTol );
      span.push_back( maxZ + bcyl->outerRadius()+zTol );
      span.push_back( minPhi - phiTol );
      span.push_back( maxPhi + phiTol );
    } else if (cyl) {
      span.push_back( minZ - cyl->outerRadius()-zTol );
      span.push_back( maxZ + cyl->outerRadius()+zTol );
      span.push_back( minPhi - phiTol );
      span.push_back( maxPhi + phiTol );
    } else {
      log << MSG::ERROR  << name() <<" volume shape not recognized: "<< (*objs)[iobj]->name() << endreq;
      for (int i=0; i<4; i++) span.push_back(0.);
    }
    const Muon::Span* newSpan=new Muon::Span(span);
    volSpans.push_back(newSpan);
  }

  if (volSpans.size()>0) {
    const std::vector<const Muon::Span*>* spans = new std::vector<const Muon::Span*>(volSpans);
    return spans;
  }
  return 0;
}

const Trk::TrackingVolume* Muon::MuonTrackingGeometryBuilder::processVolume(const Trk::Volume* vol,unsigned int etaN ,unsigned int phiN, std::string volumeName) const
{
  MsgStream log(msgSvc(), name());

  const Trk::TrackingVolume* tVol = 0;

  // partitions ? include protection against wrong setup
  if (etaN < 1 || phiN < 1) {
    log << MSG::ERROR << name() << "wrong partition setup" << endreq;
    etaN = 1;
    phiN = 1;
  }
  if ( etaN * phiN > 1 ) {  // partition
    const Trk::CylinderVolumeBounds* cyl=dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
    if (!cyl) {
      log << MSG::ERROR << " process volume: volume cylinder boundaries not retrieved, return 0 " << endreq;
      return 0; 
    }
    // subvolume boundaries
    Trk::CylinderVolumeBounds* subBds;
    if ( phiN > 1 && etaN > 1 ) { 
      subBds = new Trk::CylinderVolumeBounds(cyl->innerRadius(),
					     cyl->outerRadius(),
                                             M_PI/phiN, 
					     (cyl->halflengthZ())/etaN);
    } else if ( phiN > 1 ) {
      subBds = new Trk::CylinderVolumeBounds(cyl->innerRadius(),
					     cyl->outerRadius(),
                                             M_PI/phiN, 
					     cyl->halflengthZ());
    } else {
      subBds = new Trk::CylinderVolumeBounds(cyl->innerRadius(),
					     cyl->outerRadius(),
					     (cyl->halflengthZ())/etaN);
    } 
    // create subvolumes & BinnedArray
    std::vector<Trk::TrackingVolumeOrderPosition> subVolumes;
    for (unsigned int eta = 0; eta < etaN; eta++) {
      for (unsigned int phi = 0; phi < phiN; phi++) {
        // define subvolume
        double posZ = vol->center()[2]+subBds->halflengthZ() * (2*eta+1-etaN) ;
        HepTransform3D* transf = new HepTransform3D( HepRotateZ3D(2*M_PI/phiN*phi)*HepTranslateZ3D(posZ));
        const Trk::Volume subVol(transf, subBds);     
        // enclosed muon objects ? 
        std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( &subVol );
	std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) ; 
        const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( subVol,
								   m_muonMaterial,
								   m_muonMagneticField,
								   detVols,
								   volName );
        // reference position 
	Trk::GlobalPosition gp(subBds->outerRadius(),0.,0.);
        subVolumes.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
                                                             new Trk::GlobalPosition((*transf)*gp)));
      }
    }
    Trk::BinUtility2DPhiZ* volBinUtil=new Trk::BinUtility2DPhiZ(phiN,etaN,subBds->outerRadius(),cyl->halflengthZ(),-M_PI, new HepTransform3D());
    Trk::BinnedArray2D<Trk::TrackingVolume>* subVols=new Trk::BinnedArray2D<Trk::TrackingVolume>(subVolumes,volBinUtil);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    0,subVols,
				    volumeName);
  } else {
    // enclosed muon objects ? 
    std::vector<const Trk::DetachedTrackingVolume*>* muonObjs = getDetachedObjects( vol );

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    muonObjs,
				    volumeName);
  }

  return tVol;
} 

std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonTrackingGeometryBuilder::getDetachedObjects(const Trk::Volume* vol ) const
{

  std::vector<const Trk::DetachedTrackingVolume*>* detTVs = 0;

  if (!m_stations && !m_inertObjs) return detTVs;
  
  // get min/max Z/Phi from volume (allways a cylinder volume )
  Trk::GlobalPosition center = vol->center();  
  const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  if (!cyl) return 0; 

  double zMin = center[2] - cyl->halflengthZ(); 
  double zMax = center[2] + cyl->halflengthZ();
  double pMin = 0.;
  double pMax = +2*M_PI;
  bool phiLim = false;
  if (cyl->halfPhiSector() < M_PI) { 
    pMin = center.phi() - cyl->halfPhiSector() + M_PI; 
    pMax = center.phi() + cyl->halfPhiSector() + M_PI;
    phiLim = true;
  } 
   
  std::vector<const Trk::DetachedTrackingVolume*> detached;
  // active
  if (m_stationSpan) {
    for (unsigned int i=0; i<m_stationSpan->size() ; i++) {
      const Muon::Span s = *((*m_stationSpan)[i]);
      if ( s[0] <= zMax && s[1] >= zMin && m_stations->size()>i) {
	if (phiLim) {
	  if (pMin>=0 && pMax<=2*M_PI) {
	    if ( s[2]<=s[3] && s[2] <= pMax && s[3] >= pMin ) detached.push_back((*m_stations)[i]);
	    if ( s[2]>s[3] && (s[2] <= pMax || s[3] >= pMin) ) detached.push_back((*m_stations)[i]);
	  } else if (pMin < 0) {
	    if ( s[2]<=s[3] && (s[2] <= pMax || s[3] >= pMin+2*M_PI) ) detached.push_back((*m_stations)[i]);
	    if ( s[2]>s[3]  ) detached.push_back((*m_stations)[i]);
	  } else if (pMax > 2*M_PI) {
	    if ( s[2]<=s[3] && (s[2] <= pMax-2*M_PI || s[3] >= pMin) ) detached.push_back((*m_stations)[i]);
	    if ( s[2]>s[3]  ) detached.push_back((*m_stations)[i]);
	  }
	} else {
	  detached.push_back((*m_stations)[i]);
	}
      } 
    }
  }
  // passive
  if (m_inertSpan) {
    for (unsigned int i=0; i<m_inertSpan->size() ; i++) {
      const Muon::Span s = (*(*m_inertSpan)[i]);
      if ( s[0] <= zMax && s[1] >= zMin && m_inertObjs->size()>i) {
	if (phiLim) {
	  if (pMin>=0 && pMax<=2*M_PI) {
	    if ( s[2]<=s[3] && s[2] <= pMax && s[3] >= pMin ) detached.push_back((*m_inertObjs)[i]);
	    if ( s[2]>s[3] && (s[2] <= pMax || s[3] >= pMin) ) detached.push_back((*m_inertObjs)[i]);
	  } else if (pMin < 0) {
	    if ( s[2]<=s[3] && (s[2] <= pMax || s[3] >= pMin+2*M_PI) ) detached.push_back((*m_inertObjs)[i]);
	    if ( s[2]>s[3]  ) detached.push_back((*m_inertObjs)[i]);
	  } else if (pMax > 2*M_PI) {
	    if ( s[2]<=s[3] && (s[2] <= pMax-2*M_PI || s[3] >= pMin) ) detached.push_back((*m_inertObjs)[i]);
	    if ( s[2]>s[3]  ) detached.push_back((*m_inertObjs)[i]);
	  }
	} else {
	  detached.push_back((*m_inertObjs)[i]);
	}
      } 
    }
  }
  if (detached.size()>0) detTVs = new std::vector<const Trk::DetachedTrackingVolume*>(detached); 
  return detTVs;
}
