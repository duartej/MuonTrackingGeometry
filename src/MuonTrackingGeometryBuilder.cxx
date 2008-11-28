//////////////////////////////////////////////////////////////////
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
#include "TrkDetDescrUtils/BinUtility1DZZ.h"
#include "TrkDetDescrUtils/BinUtility1DF.h"
#include "TrkDetDescrUtils/BinUtility1DH.h"
#include "TrkDetDescrUtils/BinUtility2DPhiZ.h"
#include "TrkDetDescrUtils/BinUtility2DZPhi.h"
#include "TrkDetDescrUtils/BinUtility2DZF.h"
#include "TrkDetDescrUtils/BinUtility3DZFH.h"
#include "TrkDetDescrUtils/BinnedArray.h"
#include "TrkDetDescrUtils/BinnedArray1D.h"
#include "TrkDetDescrUtils/BinnedArray2D.h"
#include "TrkDetDescrUtils/BinnedArray3D.h"
//#include "TrkDetDescrUtils/BinningType.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkDetDescrTools/TrackingVolumeArrayCreator.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/DoubleTrapezoidVolumeBounds.h"
#include "TrkVolumes/BevelledCylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/SubtractedVolumeBounds.h"
#include "TrkVolumes/CombinedVolumeBounds.h"
#include "TrkVolumes/SimplePolygonBrepVolumeBounds.h"
#include "TrkVolumes/PrismVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
//#include "TrkGeometry/SimplifiedMaterialProperties.h"
#include "TrkMagFieldInterfaces/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/TrackingGeometry.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/GlueVolumesDescriptor.h"
#include "TrkGeometry/BendingCorrector.h"
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
  m_magFieldTool("Trk::MagneticFieldTool/AtlasMagneticFieldTool"),
  m_stationBuilder("Muon::MuonStationBuilder/MuonStationBuilder"),
  m_inertBuilder("Muon::MuonInertMaterialBuilder/MuonInertMaterialBuilder"),
  m_trackingVolumeArrayCreator("Trk::TrackingVolumeArrayCreator/TrackingVolumeArrayCreator"),
  m_trackingVolumeHelper("Trk::TrackingVolumeHelper/TrackingVolumeHelper"),
  m_trackingVolumeDisplayer("Trk::TrackingVolumeDisplayer/TrackingVolumeDisplayer"),
  m_tvSvc("Trk::TrackingVolumesSvc/TrackingVolumesSvc",n),
  m_muonSimple(false),
  m_loadMSentry(false),
  m_muonActive(true),
  m_muonInert(true),
  m_innerBarrelRadius(4251.),
  m_outerBarrelRadius(13000.),
  m_barrelZ(6750.),
  m_innerEndcapZ(12900.),
  m_outerEndcapZ(25000.),
  m_beamPipeRadius(70.),
  m_innerShieldRadius(908.),
  m_outerShieldRadius(1476.),
  m_diskShieldZ(6910.),
  m_barrelEtaPartition(9),
  m_innerEndcapEtaPartition(3),
  m_outerEndcapEtaPartition(3),
  m_phiPartition(16),
  m_adjustStatic(true),
  m_static3d(true),
  m_blendInertMaterial(true),
  m_alignTolerance(0.),
  m_chronoStatSvc( "ChronoStatSvc", n )
{
  m_stationSpan = 0;
  m_inertSpan = 0;
  m_stations = 0;
  m_inertObjs = 0;

  declareInterface<Trk::IGeometryBuilder>(this);

  declareProperty("SimpleMuonGeometry",               m_muonSimple);  
  declareProperty("LoadMSEntry",                      m_loadMSentry);  
  declareProperty("BuildActiveMaterial",              m_muonActive);  
  declareProperty("BuildInertMaterial",               m_muonInert);
  declareProperty("MagneticFieldTool",                m_magFieldTool);  
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
  declareProperty("AdjustStatic",                   m_adjustStatic);
  declareProperty("StaticPartition3D",              m_static3d);
  declareProperty("BlendInertMaterial",             m_blendInertMaterial);
  declareProperty("AlignmentPositionTolerance",     m_alignTolerance);
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
    if (s.isFailure()) log << MSG::INFO << " failing to initialize ?" << endreq;


    // Retrieve the magnetic field tool   ----------------------------------------------------    
    if (m_magFieldTool.retrieve().isFailure())
    {
      log << MSG::FATAL << "Failed to retrieve tool " << m_magFieldTool << endreq;
      return StatusCode::FAILURE;
    } else
      log << MSG::INFO << "Retrieved tool " << m_magFieldTool << endreq;


    // Retrieve the tracking volume helper   -------------------------------------------------    
    if (m_trackingVolumeHelper.retrieve().isFailure())
    {
      log << MSG::FATAL << "Failed to retrieve tool " << m_trackingVolumeHelper << endreq;
      return StatusCode::FAILURE;
    } else
      log << MSG::INFO << "Retrieved tool " << m_trackingVolumeHelper << endreq;

    // Retrieve the tracking volume array creator   -------------------------------------------    
    if (m_trackingVolumeArrayCreator.retrieve().isFailure())
    {
      log << MSG::FATAL << "Failed to retrieve tool " << m_trackingVolumeArrayCreator << endreq;
      return StatusCode::FAILURE;
    } else
      log << MSG::INFO << "Retrieved tool " << m_trackingVolumeArrayCreator << endreq;

    
    // Retrieve the station builder (if configured) -------------------------------------------    
    if (m_muonActive) { 

      if (m_stationBuilder.retrieve().isFailure())
      {
          log << MSG::ERROR << "Failed to retrieve tool " << m_stationBuilder << endreq;
          log <<" Creation of stations might fail." << endreq;
      } else
          log << MSG::INFO << "Retrieved tool " << m_trackingVolumeArrayCreator << endreq;
    }

    // Retrieve the inert material builder builder (if configured) -------------------------------------------    

    if (m_muonInert || m_blendInertMaterial) {

      if (m_inertBuilder.retrieve().isFailure())
      {
          log << MSG::ERROR << "Failed to retrieve tool " << m_inertBuilder << endreq;
          log <<"Creation of inert material objects might fail." << endreq;
      } else
          log << MSG::INFO << "Retrieved tool " << m_trackingVolumeArrayCreator << endreq;
    }

    if (m_loadMSentry ) {
      // Retrieve the tracking volumes Svc (if configured) -------------------------------------------    
      if (m_tvSvc.retrieve().isFailure()) {
	log << MSG::WARNING << "Failed to load " << m_tvSvc << " switch to default volume size" << endreq;
	m_loadMSentry = false;
      }
    }

    if (m_chronoStatSvc.retrieve().isFailure()) {
      log << MSG::WARNING << "Could not retrieve Tool " << m_chronoStatSvc << ". Exiting."<<endreq;
    }
    
        
    log << MSG::INFO  << name() <<" initialize() successful" << endreq;    
    
  return StatusCode::SUCCESS;
}


const Trk::TrackingGeometry* Muon::MuonTrackingGeometryBuilder::trackingGeometry(const Trk::TrackingVolume* tvol) const
{

  MsgStream log( msgSvc(), name() );
  
  log << MSG::INFO  << name() <<" building tracking geometry" << endreq;    
  m_chronoStatSvc->chronoStart("MS::build-up");
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // check setup
  if (m_muonInert && m_blendInertMaterial) {
    if (!m_adjustStatic || !m_static3d) {
      log << MSG::INFO  << name() <<" diluted inert material hardcoded for 3D volume frame, adjusting setup" << endreq;
      m_adjustStatic = true;
      m_static3d = true;
    }
  }
  // process muon material objects
  if (m_muonActive && m_stationBuilder && !m_stations) m_stations = m_stationBuilder->buildDetachedTrackingVolumes();
  if (m_muonInert && m_inertBuilder && !m_inertObjs) m_inertObjs = m_inertBuilder->buildDetachedTrackingVolumes();
  //if (m_inertObjs && m_blendInertMaterial) {
  //  //m_inertBlend.resize(m_inertObjs->size());
  //  //getVolumeFractions();    
  //  getDilutingFactors();    
  //}
  
  // find object's span with tolerance for the alignment 
  if (!m_stationSpan) m_stationSpan = findVolumesSpan(m_stations, 100.*m_alignTolerance, m_alignTolerance*deg);
  if (!m_inertSpan)   m_inertSpan = findVolumesSpan(m_inertObjs,0.,0.);
 
  /*
  if (m_inertSpan) {
    for (unsigned int ii=0;ii<m_inertSpan->size();ii++) {
      const Muon::Span* s=(*m_inertSpan)[ii];
      log << MSG::DEBUG << (*m_inertObjs)[ii]->name()<<":"<<
	(*s)[0]<<","<<(*s)[1]<<","<<(*s)[2]<<","<<(*s)[3]<<","<<(*s)[4]<<","<<(*s)[5]<<endreq;
    }
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
  
  Trk::MagneticFieldProperties muonMagneticFieldProperties(&(*m_magFieldTool), Trk::RealisticField);    
  m_muonMagneticField = muonMagneticFieldProperties;
  
  // dummy substructures
  const Trk::LayerArray* dummyLayers = 0;
  const Trk::TrackingVolumeArray* dummyVolumes = 0;
  Trk::BendingCorrector* bcorr = 0;
  
  if (m_muonSimple) {             
    Trk::VolumeBounds* globalBounds = new Trk::CylinderVolumeBounds(m_outerBarrelRadius,
								    m_outerEndcapZ);
    Trk::TrackingVolume* topVolume = new Trk::TrackingVolume(new HepTransform3D,
							     globalBounds,
							     m_muonMaterial,
							     m_muonMagneticField,
							     dummyLayers,dummyVolumes,
							     "GlobalVolume",bcorr);
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
  const Trk::TrackingVolume* enclosed = 0;
  Trk::CylinderVolumeBounds* enclosedBounds = 0;
  const Trk::TrackingVolume* negDiskShield = 0;
  Trk::CylinderVolumeBounds* negDiskShieldBounds = 0;
  const Trk::TrackingVolume* posDiskShield = 0;
  Trk::CylinderVolumeBounds* posDiskShieldBounds = 0;
  const Trk::TrackingVolume* negInnerShield = 0;
  const Trk::TrackingVolume* posInnerShield = 0;
  Trk::CylinderVolumeBounds* negInnerShieldBounds = 0;
  Trk::CylinderVolumeBounds* posInnerShieldBounds = 0;
  const Trk::TrackingVolume* negOuterShield = 0;
  const Trk::TrackingVolume* posOuterShield = 0;
  Trk::CylinderVolumeBounds* negOuterShieldBounds = 0;
  Trk::CylinderVolumeBounds* posOuterShieldBounds = 0;
// if input, redefine dimensions to fit expected MS entry 
  if (tvol){
    bool msEntryDefined = false;
    if ( tvol->volumeName() == "MuonSpectrometerEntrance" ) msEntryDefined = true; 
    // get dimensions
    const Trk::CylinderVolumeBounds* enclosedDetectorBounds 
      = dynamic_cast<const Trk::CylinderVolumeBounds*>(&(tvol->volumeBounds()));
    
    double enclosedDetectorHalfZ = enclosedDetectorBounds->halflengthZ();
    double enclosedDetectorOuterRadius = enclosedDetectorBounds->outerRadius();
    // get subvolumes at navigation level and check THEIR dimensions      
    const Trk::GlueVolumesDescriptor& enclosedDetGlueVolumes = tvol->glueVolumesDescriptor();
    std::vector<const Trk::TrackingVolume*> enclosedCentralFaceVolumes = enclosedDetGlueVolumes.glueVolumes(Trk::cylinderCover);
    std::vector<const Trk::TrackingVolume*> enclosedNegativeFaceVolumes = enclosedDetGlueVolumes.glueVolumes(Trk::negativeFaceXY);
    std::vector<const Trk::TrackingVolume*> enclosedPositiveFaceVolumes = enclosedDetGlueVolumes.glueVolumes(Trk::positiveFaceXY);
    if (enclosedCentralFaceVolumes.size()) {
      const Trk::CylinderVolumeBounds* cylR = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(enclosedCentralFaceVolumes[0]->volumeBounds()));
      if (cylR && cylR->outerRadius() != enclosedDetectorOuterRadius ) {
        enclosedDetectorOuterRadius = cylR->outerRadius();
        log << MSG::WARNING << name() << " enclosed volume envelope outer radius does not correspond to radius of glue volumes : adjusted " << endreq;
      }
    }
    if (enclosedNegativeFaceVolumes.size() && enclosedPositiveFaceVolumes.size()) {
      double negZ = -enclosedDetectorHalfZ;
      double posZ =  enclosedDetectorHalfZ;
      const Trk::CylinderVolumeBounds* cylN = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(enclosedNegativeFaceVolumes[0]->volumeBounds()));
      if (cylN) negZ = enclosedNegativeFaceVolumes[0]->center()[2]-cylN->halflengthZ();
      const Trk::CylinderVolumeBounds* cylP = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(enclosedPositiveFaceVolumes[0]->volumeBounds()));
      if (cylP) posZ = enclosedPositiveFaceVolumes[0]->center()[2]+cylP->halflengthZ();
      if ( fabs(negZ+enclosedDetectorHalfZ)>0.001 ||  fabs(posZ-enclosedDetectorHalfZ)>0.001 ) {
	log << MSG::WARNING << name() << " enclosed volume envelope z dimension does not correspond to that of glue volumes " << endreq;
        if ( fabs(negZ+posZ)<0.001) {
          enclosedDetectorHalfZ = posZ;
          log << MSG::WARNING << name() << " z adjusted " << endreq;
	} else {
	  log << MSG::ERROR << name() << "assymetric Z dimensions - cannot recover "<< negZ <<","<< posZ << endreq;
          return 0;
	}
     }
    }
    //
    

    // 
    log << MSG::INFO  << name() <<" dimensions of enclosed detectors (halfZ,outerR):"
	<<  enclosedDetectorHalfZ<<","<< enclosedDetectorOuterRadius << endreq;    
    // check if input makes sense - gives warning if cuts into muon envelope
    // adjust radius
    const Trk::TrackingVolume* barrelR = 0;
    if ( enclosedDetectorOuterRadius > m_innerBarrelRadius ) {
      log << MSG::WARNING  << name() <<" enclosed volume too wide, cuts into muon envelope, abandon :R:" <<enclosedDetectorOuterRadius << endreq;
      return 0;   
    } else {    
      if ( enclosedDetectorOuterRadius < m_innerBarrelRadius ) {
	Trk::CylinderVolumeBounds* barrelRBounds = new Trk::CylinderVolumeBounds(enclosedDetectorOuterRadius,
										 m_innerBarrelRadius,
										 enclosedDetectorHalfZ);
	const Trk::TrackingVolume* barrelRBuffer = new Trk::TrackingVolume(new HepTransform3D(),
									   barrelRBounds,
									   m_muonMaterial,
									   m_muonMagneticField,
									   dummyLayers,dummyVolumes,
									   "BarrelRBuffer", bcorr);
	
	log << MSG::INFO << name() << "glue barrel R buffer + ID/Calo volumes" << endreq;
	barrelR = m_trackingVolumeHelper->glueTrackingVolumeArrays(*barrelRBuffer,Trk::tubeInnerCover,
								   *tvol, Trk::cylinderCover,
								   "All::Gaps::BarrelR");    
      } else  barrelR = tvol;
      
    }
    // adjust z
    if ( enclosedDetectorHalfZ > m_barrelZ ) {
      log << MSG::WARNING  << name() <<" enclosed volume too long, cuts into muon envelope, abandon :Z:"<< enclosedDetectorHalfZ << endreq;    
      return 0;
    } else {    
      if ( enclosedDetectorHalfZ < m_barrelZ ) {
	Trk::CylinderVolumeBounds* barrelZPBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
										 0.5*(m_barrelZ - enclosedDetectorHalfZ) );
	Trk::CylinderVolumeBounds* barrelZMBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
										 0.5*(m_barrelZ - enclosedDetectorHalfZ) );
	double zbShift =  0.5*(m_barrelZ + enclosedDetectorHalfZ);
	const Trk::TrackingVolume* barrelZPBuffer = new Trk::TrackingVolume(new HepTransform3D(HepTranslateZ3D(zbShift)),
									    barrelZPBounds,
									    m_muonMaterial,
									    m_muonMagneticField,
									    dummyLayers,dummyVolumes,
									    "BarrelRZPosBuffer", bcorr);
	const Trk::TrackingVolume* barrelZMBuffer = new Trk::TrackingVolume(new HepTransform3D(HepTranslateZ3D(-zbShift)),
									    barrelZMBounds,
									    m_muonMaterial,
									    m_muonMagneticField,
									    dummyLayers,dummyVolumes,
									    "BarrelRZNegBuffer", bcorr);
	
	log << MSG::INFO << name() << "glue barrel R  + barrel Z buffer" << endreq;
	const Trk::TrackingVolume* barrelZP = m_trackingVolumeHelper->glueTrackingVolumeArrays(*barrelR, Trk::positiveFaceXY,
											       *barrelZPBuffer,Trk::negativeFaceXY, 
											       "All::Gaps::BarrelZP");    
        // set name
	std::string nameEncl = msEntryDefined ? "All::Gaps::Barrel" : "MuonSpectrometerEntrance" ;
	enclosed = m_trackingVolumeHelper->glueTrackingVolumeArrays(*barrelZP, Trk::negativeFaceXY,
								    *barrelZMBuffer,Trk::positiveFaceXY, 
								    nameEncl);    
 
      } else enclosed = barrelR;
      
    }
    
  } else {     // no input, create the enclosed volume
    if (m_loadMSentry && m_tvSvc ) {
      const Trk::CylinderVolumeBounds* cylMS = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(m_tvSvc->volume(ITrackingVolumesSvc::MuonSpectrometerEntryLayer).volumeBounds()));
      if ( cylMS ) {
        double rMS = cylMS->outerRadius();
        double zMS = cylMS->halflengthZ();
        if ( rMS <= m_innerBarrelRadius && zMS <= m_barrelZ ) {
	  enclosedBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
							 m_barrelZ);
	} else {
	  log << MSG::INFO << name() << " input MSEntrance size (R,Z:"<< rMS <<","<<zMS<<") clashes with MS material, switch to default values (R,Z:" << m_innerBarrelRadius <<","<< m_barrelZ << ")" << endreq;
	}
      }
    }
    
    if (!enclosedBounds) enclosedBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
									m_barrelZ);
    enclosed = new Trk::TrackingVolume(new HepTransform3D(),
				       enclosedBounds,
				       m_muonMaterial,
				       m_muonMagneticField,
				       dummyLayers,dummyVolumes,
				       "MuonSpectrometerEntrance", bcorr);
  }
  
  // create central volume ("enclosed" + disk shields ) - this is to allow safe gluing with 3D MS binning
  
  negDiskShieldBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
						      0.5*(m_diskShieldZ - m_barrelZ) );
  Hep3Vector negDiskShieldPosition(0.,0.,-0.5*(m_diskShieldZ+m_barrelZ));
  Trk::Volume negDiskVol(new HepTransform3D(Trk::s_idRotation,negDiskShieldPosition),negDiskShieldBounds);
  negDiskShield = processVolume(&negDiskVol,1,1,"Muons::Detectors::NegativeDiskShield");

  posDiskShieldBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
						      0.5*(m_diskShieldZ - m_barrelZ) );
  Hep3Vector posDiskShieldPosition(0.,0., 0.5*(m_diskShieldZ+m_barrelZ));
  Trk::Volume posDiskVol(new HepTransform3D(Trk::s_idRotation,posDiskShieldPosition),posDiskShieldBounds);
  posDiskShield = processVolume(&posDiskVol,1,1,"Muons::Detectors::PositiveDiskShield");
	   
  log << MSG::INFO << name() << "glue enclosed  + disk shields" << endreq;
  const Trk::TrackingVolume* centralP = m_trackingVolumeHelper->glueTrackingVolumeArrays(*enclosed, Trk::positiveFaceXY,
											 *posDiskShield,Trk::negativeFaceXY, 
											 "Muon::Container::CentralP");    
  const Trk::TrackingVolume* central  = m_trackingVolumeHelper->glueTrackingVolumeArrays(*centralP, Trk::negativeFaceXY,
											 *negDiskShield,Trk::positiveFaceXY, 
											 "Muon::Container::Central");    
  
// define basic volumes
  if (m_adjustStatic) { getZParts(); getPhiParts(); getHParts();}

// muon barrel
  barrelBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
					       m_outerBarrelRadius,
					       m_diskShieldZ);
  Trk::Volume barrelVol(new HepTransform3D(),barrelBounds);
// process volume
// barrel
  if (m_adjustStatic && m_static3d) muonBarrel = processVolume( &barrelVol,0,"Muon::Detectors::Barrel"); 
  else if (m_adjustStatic) muonBarrel = processVolume( &barrelVol,-1,"Muon::Detectors::Barrel"); 
  else muonBarrel = processVolume( &barrelVol,m_barrelEtaPartition,m_phiPartition,"Muon::Detectors::Barrel"); 
// inner Endcap
   double innerEndcapZHalfSize = 0.5*(m_innerEndcapZ - m_diskShieldZ);
   negativeInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_innerShieldRadius,
                                                             m_outerBarrelRadius,
                                                             innerEndcapZHalfSize);
   Hep3Vector negInnerEndcapPosition(0.,0.,-m_diskShieldZ-innerEndcapZHalfSize);
   HepTransform3D* negInnerEndcapTransf = new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition);
   Trk::Volume negIECVol(negInnerEndcapTransf,negativeInnerEndcapBounds);
   if (m_adjustStatic && m_static3d) negativeMuonInnerEndcap = processVolume( &negIECVol,1, "Muon::Detectors::NegativeInnerEndcap" ); 
   else if (m_adjustStatic) negativeMuonInnerEndcap = processVolume( &negIECVol,-1, "Muon::Detectors::NegativeInnerEndcap" ); 
   else negativeMuonInnerEndcap = processVolume( &negIECVol,m_innerEndcapEtaPartition,m_phiPartition,
					    "Muon::Detectors::NegativeInnerEndcap" ); 
//
   positiveInnerEndcapBounds = new Trk::CylinderVolumeBounds(m_innerShieldRadius,
                                                             m_outerBarrelRadius,
                                                             innerEndcapZHalfSize);
   Hep3Vector posInnerEndcapPosition(0.,0.,m_diskShieldZ+innerEndcapZHalfSize);
   HepTransform3D* posInnerEndcapTransf = new HepTransform3D(Trk::s_idRotation,posInnerEndcapPosition);
   Trk::Volume posIECVol(posInnerEndcapTransf,positiveInnerEndcapBounds);
   if (m_adjustStatic && m_static3d) positiveMuonInnerEndcap = processVolume( &posIECVol,1, "Muon::Detectors::PositiveInnerEndcap" ); 
   else if (m_adjustStatic) positiveMuonInnerEndcap = processVolume( &posIECVol,-1, "Muon::Detectors::PositiveInnerEndcap" ); 
   else positiveMuonInnerEndcap = processVolume( &posIECVol,m_innerEndcapEtaPartition,m_phiPartition,
					    "Muon::Detectors::PositiveInnerEndcap" ); 
// inner shields   
   negInnerShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_innerShieldRadius,
							innerEndcapZHalfSize);
   Trk::Volume negisVol(new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition),negInnerShieldBounds);
   negInnerShield = processVolume(&negisVol,1,1,"Muons::Detectors::NegativeInnerShield");

   posInnerShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_innerShieldRadius,
							innerEndcapZHalfSize);
   Trk::Volume posisVol(new HepTransform3D(Trk::s_idRotation,posInnerEndcapPosition),posInnerShieldBounds);
   posInnerShield = processVolume(&posisVol,1,1,"Muons::Detectors::PositiveInnerShield");

// outer Endcap
   double outerEndcapZHalfSize = 0.5*(m_outerEndcapZ - m_innerEndcapZ);
   negativeOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_outerShieldRadius,
                                                             m_outerBarrelRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector negOuterEndcapPosition(0.,0.,-m_innerEndcapZ-outerEndcapZHalfSize);
   HepTransform3D* negOuterEndcapTransf = new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition);
   Trk::Volume negOECVol(negOuterEndcapTransf,negativeOuterEndcapBounds);
   if (m_adjustStatic && m_static3d) negativeMuonOuterEndcap = processVolume( &negOECVol,2,"Muon::Detectors::NegativeOuterEndcap" ); 
   else if (m_adjustStatic) negativeMuonOuterEndcap = processVolume( &negOECVol,-1,"Muon::Detectors::NegativeOuterEndcap" ); 
   else negativeMuonOuterEndcap = processVolume( &negOECVol,m_outerEndcapEtaPartition,m_phiPartition,
						 "Muon::Detectors::NegativeOuterEndcap" ); 
//
   positiveOuterEndcapBounds = new Trk::CylinderVolumeBounds(m_outerShieldRadius,
                                                             m_outerBarrelRadius,
                                                             outerEndcapZHalfSize);
   Hep3Vector posOuterEndcapPosition(0.,0., m_innerEndcapZ+outerEndcapZHalfSize);
   HepTransform3D* posOuterEndcapTransf = new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition);
   Trk::Volume posOECVol(posOuterEndcapTransf,positiveOuterEndcapBounds);
   if (m_adjustStatic && m_static3d) positiveMuonOuterEndcap = processVolume( &posOECVol,2,"Muon::Detectors::PositiveOuterEndcap" ); 
   else if (m_adjustStatic) positiveMuonOuterEndcap = processVolume( &posOECVol,-1,"Muon::Detectors::PositiveOuterEndcap" ); 
   else positiveMuonOuterEndcap = processVolume( &posOECVol,m_outerEndcapEtaPartition,m_phiPartition,
						 "Muon::Detectors::PositiveOuterEndcap" ); 
// outer shields   
   negOuterShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_outerShieldRadius,
							outerEndcapZHalfSize);
   Trk::Volume negosVol(new HepTransform3D(Trk::s_idRotation,negOuterEndcapPosition),negOuterShieldBounds);
   negOuterShield = processVolume(&negosVol,1,1,"Muons::Detectors::NegativeOuterShield");
   posOuterShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_outerShieldRadius,
							outerEndcapZHalfSize);
   Trk::Volume pososVol(new HepTransform3D(Trk::s_idRotation,posOuterEndcapPosition),posOuterShieldBounds);
   posOuterShield = processVolume(&pososVol,1,1,"Muons::Detectors::PositiveOuterShield");

// beamPipe
   negBeamPipeBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  outerEndcapZHalfSize+innerEndcapZHalfSize);
   posBeamPipeBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
                                                  outerEndcapZHalfSize+innerEndcapZHalfSize);
   Hep3Vector posBeamPipePosition(0.,0., m_outerEndcapZ-innerEndcapZHalfSize-outerEndcapZHalfSize);
   Hep3Vector negBeamPipePosition(0.,0.,-m_outerEndcapZ+innerEndcapZHalfSize+outerEndcapZHalfSize);
   Trk::Volume negbpVol(new HepTransform3D(Trk::s_idRotation,negBeamPipePosition),negBeamPipeBounds);
   negBeamPipe = processVolume(&negbpVol,1,1,"Muons::Gaps::NegativeBeamPipe");
   Trk::Volume posbpVol(new HepTransform3D(Trk::s_idRotation,posBeamPipePosition),posBeamPipeBounds);
   posBeamPipe = processVolume(&posbpVol,1,1,"Muons::Gaps::PositiveBeamPipe");

   log << MSG::INFO  << name() <<" volumes defined " << endreq;    
//
// glue volumes at navigation level, create enveloping volume 
// radially
// central + barrel
   log << MSG::INFO << name() << "glue barrel+enclosed volumes" << endreq;
   const Trk::TrackingVolume* barrel = m_trackingVolumeHelper->glueTrackingVolumeArrays(*muonBarrel,Trk::tubeInnerCover,
											*central, Trk::cylinderCover, 
                                                                                        "All::Container::Barrel");
   //checkVolume(barrel);
// shield+outerEndcap
   log << MSG::INFO << name() << "glue shield+outerEndcap" << endreq;
   const Trk::TrackingVolume* negOuterEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*negativeMuonOuterEndcap, Trk::tubeInnerCover,
												*negOuterShield, Trk::tubeOuterCover,
												"Muon::Container::NegativeOuterEndcap");
   //checkVolume(negOuterEndcap);
   const Trk::TrackingVolume* posOuterEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*positiveMuonOuterEndcap, Trk::tubeInnerCover,
												*posOuterShield,Trk::tubeOuterCover,
												"Muon::Container::PositiveOuterEndcap");
   //checkVolume(posOuterEndcap);
// shield+innerEndcap
   log << MSG::INFO << name() << "glue shield+innerEndcap" << endreq;
   const Trk::TrackingVolume* negInnerEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*negativeMuonInnerEndcap, Trk::tubeInnerCover,
												*negInnerShield, Trk::tubeOuterCover,
												"Muon::Container::NegativeInnerEndcap");
   //checkVolume(negInnerEndcap);
   const Trk::TrackingVolume* posInnerEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*positiveMuonInnerEndcap, Trk::tubeInnerCover,
												*posInnerShield,Trk::tubeOuterCover,
												"Muon::Container::PositiveInnerEndcap");
   //checkVolume(posInnerEndcap);
// inner+outerEndcap
   log << MSG::INFO << name() << "glue inner+outerEndcap" << endreq;
   const Trk::TrackingVolume* negNavEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*negOuterEndcap, Trk::positiveFaceXY,
											      *negInnerEndcap, Trk::negativeFaceXY, 
											      "Muon::Container::NegativeEndcap");  
   const Trk::TrackingVolume* posNavEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*posInnerEndcap, Trk::positiveFaceXY,
											      *posOuterEndcap, Trk::negativeFaceXY,
											      "Muon::Container::PositiveEndcap");  
   //checkVolume(negNavEndcap);
   //checkVolume(posNavEndcap);
// beam pipe + endcaps
   log << MSG::INFO << name() << "glue beamPipe+endcaps" << endreq;
   const Trk::TrackingVolume* negEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*negNavEndcap, Trk::tubeInnerCover,
											   *negBeamPipe, Trk::cylinderCover,
										           "All::Container::NegativeEndcap");  
   const Trk::TrackingVolume* posEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*posNavEndcap, Trk::tubeInnerCover,
											   *posBeamPipe, Trk::cylinderCover,
										           "All::Container::PositiveEndcap");  
   //checkVolume(negEndcap);
   //checkVolume(posEndcap);
// barrel + endcaps
   log << MSG::INFO << name() << "glue barrel+endcaps" << endreq;

   const Trk::TrackingVolume* negDet = m_trackingVolumeHelper->glueTrackingVolumeArrays(*negEndcap, Trk::positiveFaceXY,
                                                                                        *barrel, Trk::negativeFaceXY, 
									                "All::Container::NegDet");  
   const Trk::TrackingVolume* detector = m_trackingVolumeHelper->glueTrackingVolumeArrays(*posEndcap, Trk::negativeFaceXY, 
											  *negDet, Trk::positiveFaceXY,
											  "All::Container::CompleteDetector");
//
   Trk::TrackingGeometry* trackingGeometry = new Trk::TrackingGeometry(detector,Trk::globalSearch);

   if (m_blendInertMaterial) blendMaterial();
   //log << MSG::INFO << name() << " print volume hierarchy" << endreq;
   //trackingGeometry->printVolumeHierarchy(log);

   m_chronoStatSvc->chronoStop("MS::build-up");

   log << MSG::INFO  << name() <<" returning tracking geometry " << endreq;    
   return trackingGeometry;  
}

// finalize
StatusCode Muon::MuonTrackingGeometryBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    if (m_stations) {
      for (size_t i = 0; i < m_stations->size(); i++)
	if ((*m_stations)[i]) delete (*m_stations)[i];
        else log << MSG::DEBUG << name() << " station pointer corrupted ! " << endreq; 
      delete m_stations; m_stations = 0;
    } 
    if (m_inertObjs) {
      for (size_t i = 0; i < m_inertObjs->size(); i++)
	if ((*m_inertObjs)[i])	delete (*m_inertObjs)[i];
        else log << MSG::DEBUG << name() << " inert object pointer corrupted ! " << endreq; 
      delete m_inertObjs; m_inertObjs = 0;
    } 
    if (m_stationSpan) {
      for (size_t i = 0; i < m_stationSpan->size(); i++)
        delete (*m_stationSpan)[i];
      delete m_stationSpan; m_stationSpan = 0;
    }
    if (m_inertSpan) {
      for (size_t i = 0; i < m_inertSpan->size(); i++)
        delete (*m_inertSpan)[i];
      delete m_inertSpan; m_inertSpan = 0;
    }

    for (std::map<const Trk::DetachedTrackingVolume*,std::vector<const Trk::TrackingVolume*>* >::iterator it = m_blendMap.begin();
         it != m_blendMap.end();
         ++it)
    {
      delete it->second;
    }

    m_chronoStatSvc->chronoPrint("MS::build-up");

    log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    return StatusCode::SUCCESS;
}
const Muon::Span* Muon::MuonTrackingGeometryBuilder::findVolumeSpan(const Trk::VolumeBounds* volBounds, HepTransform3D transform, double zTol, double phiTol) const
{
  MsgStream log(msgSvc(), name());

  if (!volBounds) return 0;
  // volume shape
  const Trk::CuboidVolumeBounds* box = dynamic_cast<const Trk::CuboidVolumeBounds*> (volBounds);
  const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (volBounds);
  const Trk::DoubleTrapezoidVolumeBounds* dtrd = dynamic_cast<const Trk::DoubleTrapezoidVolumeBounds*> (volBounds);
  const Trk::BevelledCylinderVolumeBounds* bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (volBounds);
  const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (volBounds);
  const Trk::SubtractedVolumeBounds* sub = dynamic_cast<const Trk::SubtractedVolumeBounds*> (volBounds);
  const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*> (volBounds);
  const Trk::SimplePolygonBrepVolumeBounds* spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*> (volBounds);
  const Trk::PrismVolumeBounds* prism = dynamic_cast<const Trk::PrismVolumeBounds*> (volBounds);

  if (sub) return findVolumeSpan(&(sub->outer()->volumeBounds()),sub->outer()->transform(),zTol,phiTol);

  if (comb) {
    const Muon::Span* s1 = findVolumeSpan(&(comb->first()->volumeBounds()),comb->first()->transform(),zTol,phiTol);     
    const Muon::Span* s2 = findVolumeSpan(&(comb->second()->volumeBounds()),comb->second()->transform(),zTol,phiTol);     
    Muon::Span scomb;
    scomb.reserve(6); 
    scomb.push_back(fmin((*s1)[0],(*s2)[0]));
    scomb.push_back(fmax((*s1)[1],(*s2)[1]));
    scomb.push_back(fmin((*s1)[2],(*s2)[2]));
    scomb.push_back(fmax((*s1)[3],(*s2)[3]));
    scomb.push_back(fmin((*s1)[4],(*s2)[4]));
    scomb.push_back(fmax((*s1)[5],(*s2)[5]));
    return new Muon::Span(scomb);
  }

  // loop over edges ...
  double minZ = m_outerEndcapZ ; double maxZ = - m_outerEndcapZ;
  double minPhi = 2*M_PI; double maxPhi = 0.;
  double minR = m_outerBarrelRadius; double maxR = 0.;
  std::vector<Trk::GlobalPosition> edges;
  edges.reserve(16);
  Muon::Span span;  
  span.reserve(6); 
  
  double cylZcorr = 0.; 
  if (box) {
    edges.push_back( Trk::GlobalPosition(box->halflengthX(),box->halflengthY(),box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(-box->halflengthX(),box->halflengthY(),box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(box->halflengthX(),-box->halflengthY(),box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(-box->halflengthX(),-box->halflengthY(),box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(box->halflengthX(),box->halflengthY(),-box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(-box->halflengthX(),box->halflengthY(),-box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(box->halflengthX(),-box->halflengthY(),-box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(-box->halflengthX(),-box->halflengthY(),-box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(0.,0.,-box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(0.,0., box->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(-box->halflengthX(),0.,0.) );
    edges.push_back( Trk::GlobalPosition( box->halflengthX(),0.,0.) );
    edges.push_back( Trk::GlobalPosition(0.,-box->halflengthY(),0.) );
    edges.push_back( Trk::GlobalPosition(0., box->halflengthY(),0.) );
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
    edges.push_back( Trk::GlobalPosition(0.,0.,-trd->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(0.,0., trd->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(-0.5*(trd->minHalflengthX()+trd->maxHalflengthX()),0.,0.) );
    edges.push_back( Trk::GlobalPosition( 0.5*(trd->minHalflengthX()+trd->maxHalflengthX()),0.,0.) );
    edges.push_back( Trk::GlobalPosition(0.,-trd->halflengthY(),0.) );
    edges.push_back( Trk::GlobalPosition(0., trd->halflengthY(),0.) );
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
    edges.push_back( Trk::GlobalPosition(0.,0.,-dtrd->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(0.,0., dtrd->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(0.,-2*dtrd->halflengthY1(),0.) );
    edges.push_back( Trk::GlobalPosition(0., 2*dtrd->halflengthY2(),0.) );
  }
  if (bcyl) {
    edges.push_back( Trk::GlobalPosition(0.,0., bcyl->halflengthZ()) );
    Trk::GlobalPosition gp = transform.getRotation()*edges.back();
    cylZcorr = bcyl->outerRadius()*gp.perp()/bcyl->halflengthZ(); 
    edges.push_back( Trk::GlobalPosition(0.,0.,-bcyl->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(bcyl->innerRadius(),0.,0.) );
    edges.push_back( Trk::GlobalPosition(bcyl->outerRadius(),0.,0.) );
  }
  if (cyl) {
    edges.push_back( Trk::GlobalPosition(0.,0., cyl->halflengthZ()) );
    Trk::GlobalPosition gp = transform.getRotation()*edges.back();
    cylZcorr = cyl->outerRadius()*gp.perp()/cyl->halflengthZ();
    edges.push_back( Trk::GlobalPosition(0.,0.,-cyl->halflengthZ()) );
    edges.push_back( Trk::GlobalPosition(cyl->innerRadius(),0.,0.) );
    edges.push_back( Trk::GlobalPosition(cyl->outerRadius(),0.,0.) );
  }
  if (spb) {
    const std::vector<std::pair<double, double> > vtcs = spb->xyVertices();
    for (unsigned int i=0;i<vtcs.size();i++) {
      edges.push_back( Trk::GlobalPosition(vtcs[i].first,vtcs[i].second, spb->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(vtcs[i].first,vtcs[i].second, -spb->halflengthZ()) );
    }
  }
  if (prism) {
    const std::vector<std::pair<double, double> > vtcs = prism->xyVertices();
    for (unsigned int i=0;i<vtcs.size();i++) {
      edges.push_back( Trk::GlobalPosition(vtcs[i].first,vtcs[i].second, prism->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(vtcs[i].first,vtcs[i].second, -prism->halflengthZ()) );
    }
  }
  // apply transform and get span
  double minP0=M_PI; double maxP0 = 0.; double minP1=2*M_PI; double maxP1=M_PI;
  for (unsigned int ie=0; ie < edges.size() ; ie++) {
    Trk::GlobalPosition gp = transform*edges[ie];
    double phi = gp.phi()+M_PI; 
    //log << MSG::DEBUG << "edges:"<< ie<<","<<gp<<","<< phi<< endreq;
    double rad = gp.perp();
    if ( gp[2]<minZ ) minZ = gp[2];
    if ( gp[2]>maxZ ) maxZ = gp[2];
    if ( phi < minP0 ) minP0 = phi; 
    if ( phi > maxP0 ) maxP0 = phi; 
    if ( phi < minP1 ) minP1 = phi; 
    if ( phi > maxP1 ) maxP1 = phi; 
    //if ( phi < minPhi ) minPhi = phi; 
    //if ( phi > maxPhi ) maxPhi = phi; 
    if ( rad < minR ) minR = rad; 
    if ( rad > maxR ) maxR = rad; 
    if (cyl || bcyl ) {
      double radius = 0.;
      if (cyl) radius = cyl->outerRadius();
      if (bcyl) radius = bcyl->outerRadius();
      if ( gp.perp() < radius ) {
	minPhi = 0.;
	maxPhi = 2*M_PI;
      } else {
	double dPhi = asin(radius/gp.perp());
	minPhi = fmin(minPhi, phi - dPhi);
	maxPhi = fmax(maxPhi, phi + dPhi);
      }
    }
  }
  if (maxP0>=minP0 && maxP1<minP1) { minPhi = minP0; maxPhi = maxP0; }
  else if ( maxP1>=minP1 && maxP0<minP0) { minPhi = minP1; maxPhi = maxP1; }
  else if ( maxP1 - minP0 < (maxP0 - minP1-2*M_PI) ) { minPhi = minP0; maxPhi = maxP1; }
  else { minPhi = minP1 ; maxPhi = maxP0; }  
  //if ( maxPhi-minPhi > M_PI  &&  maxPhi-minPhi < 2*M_PI ) {
  //  double medPhi = minPhi;
  //  minPhi = maxPhi;
  //  maxPhi = medPhi;
  //}
  if ( box || trd || dtrd || spb ) {
    span.push_back( minZ - zTol );  
    span.push_back( maxZ + zTol );  
    span.push_back( minPhi - phiTol );  
    span.push_back( maxPhi + phiTol );  
    span.push_back( fmax(m_beamPipeRadius+0.001, minR - zTol) );  
    span.push_back( maxR + zTol );  
  } else if (bcyl) {
    span.push_back( minZ - cylZcorr -zTol );
    span.push_back( maxZ + cylZcorr +zTol );
    span.push_back( minPhi - phiTol );
    span.push_back( maxPhi + phiTol );
    span.push_back( fmax(m_beamPipeRadius+0.001, minR - zTol) );  
    span.push_back( maxR + zTol );  
  } else if (cyl) {
    span.push_back( minZ - cylZcorr -zTol );
    span.push_back( maxZ + cylZcorr +zTol );
    span.push_back( minPhi - phiTol );
    span.push_back( maxPhi + phiTol );
    span.push_back( fmax(m_beamPipeRadius+0.001, minR - zTol) );  
    span.push_back( maxR + zTol );  
  } else {
    log << MSG::ERROR  << name() <<" volume shape not recognized: "<< endreq;
    for (int i=0; i<4; i++) span.push_back(0.);
  }
  const Muon::Span* newSpan=new Muon::Span(span);
  return newSpan;
}

const std::vector<const Muon::Span*>* Muon::MuonTrackingGeometryBuilder::findVolumesSpan(const std::vector<const Trk::DetachedTrackingVolume*>*& objs, double zTol, double phiTol) const
{
  MsgStream log(msgSvc(), name());

  if (!objs) return 0;
  std::vector<const Span*> volSpans;
  for (unsigned int iobj=0; iobj<objs->size(); iobj++) {
    HepTransform3D  transform = (*objs)[iobj]->trackingVolume()->transform();
    const Muon::Span* span = findVolumeSpan(&((*objs)[iobj]->trackingVolume()->volumeBounds()), transform, zTol, phiTol);
    //std::cout << "span:"<<(*objs)[iobj]->name()<< ","<<(*span)[0]<<","<< (*span)[1]<<","<<(*span)[2]<<","
    //<< (*span)[3]<<","<< (*span)[4]<<","<< (*span)[5]<<std::endl; 
    volSpans.push_back(span);
  }

  if (volSpans.size()>0) {
    const std::vector<const Muon::Span*>* spans = new std::vector<const Muon::Span*>(volSpans);
    return spans;
  }
  return 0;
}


const Trk::TrackingVolume* Muon::MuonTrackingGeometryBuilder::processVolume(const Trk::Volume* vol, int etaN , int phiN, std::string volumeName) const
{
  MsgStream log(msgSvc(), name());

  const Trk::TrackingVolume* tVol = 0;

  unsigned int colorCode = 12;

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
    Trk::CylinderVolumeBounds* subBds = 0;

    double phiSect = M_PI/phiN;
    double etaSect = (cyl->halflengthZ())/etaN;

    subBds = new Trk::CylinderVolumeBounds(cyl->innerRadius(), cyl->outerRadius(), phiSect,etaSect);
    const Trk::Volume* protVol= new Trk::Volume(0, subBds);     

    // create subvolumes & BinnedArray
    std::vector<Trk::TrackingVolumeOrderPosition> subVolumes;
    std::vector<const Trk::TrackingVolume*> sVols;                // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsNeg;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsPos;             // for gluing
    for (int eta = 0; eta < etaN; eta++) {
      colorCode = 26 - colorCode;
      for (int phi = 0; phi < phiN; phi++) {
	colorCode = 26 - colorCode;
        // define subvolume
        double posZ = (vol->center())[2]+ etaSect * (2.*eta+1.-etaN) ;
        HepTransform3D  transf( HepRotateZ3D( phiSect*(2*phi+1))*HepTranslateZ3D(posZ));
        const Trk::Volume* subVol= new Trk::Volume(*protVol, transf);     
        // enclosed muon objects ?   
	std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) ; 
        std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( subVol);
        const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( *subVol,
								   m_muonMaterial,
								   m_muonMagneticField,
								   detVols,
								   volName );
        // prepare blending
        if (m_blendInertMaterial && detVols) {
          for (unsigned int id=0;id<detVols->size();id++) {
            if (!(*detVols)[id]->layerRepresentation() && (*detVols)[id]->constituents()){
	      if (!m_blendMap[(*detVols)[id]]) 
		m_blendMap[(*detVols)[id]] = new std::vector<const Trk::TrackingVolume*>;
	      m_blendMap[(*detVols)[id]]->push_back(sVol);
	    }
	  }
	}  
        //delete subVol;
        sVol->registerColorCode(colorCode); 
	// reference position 
	HepPoint3D gp(subBds->outerRadius(),0.,0.);
	subVolumes.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
                                                             new HepPoint3D(transf*gp)));
        //glue subVolumes
        sVols.push_back(sVol); 
        if (eta==0)      sVolsNeg.push_back(sVol); 
        if (eta==etaN-1) sVolsPos.push_back(sVol); 
        // in phi        
        if ( phiN>1 && phi>0) {
          m_trackingVolumeHelper->glueTrackingVolumes(*sVol, Trk::tubeSectorNegativePhi,
						      *(sVols[eta*phiN+phi-1]), Trk::tubeSectorPositivePhi);
          if ( phi==phiN-1 )  m_trackingVolumeHelper->glueTrackingVolumes(*(sVols[eta*phiN]), Trk::tubeSectorNegativePhi,
                                                                          *sVol, Trk::tubeSectorPositivePhi);
	}
        // in eta
        if ( etaN>1 && eta>0) m_trackingVolumeHelper->glueTrackingVolumes(*sVol, Trk::negativeFaceXY,
						      *(sVols[(eta-1)*phiN+phi]), Trk::positiveFaceXY);        
        //
      }
    }

    Trk::BinUtility2DPhiZ* volBinUtil=new Trk::BinUtility2DPhiZ(phiN,etaN,subBds->outerRadius(),cyl->halflengthZ(),M_PI, new HepTransform3D(vol->transform()));
    delete protVol;
    Trk::BinnedArray2D<Trk::TrackingVolume>* subVols=new Trk::BinnedArray2D<Trk::TrackingVolume>(subVolumes,volBinUtil);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    0,subVols,
				    volumeName);
    // register glue volumes
    const Trk::GlueVolumesDescriptor& volGlueVolumes = tVol->glueVolumesDescriptor();
    volGlueVolumes.registerGlueVolumes(Trk::tubeInnerCover,sVols);
    volGlueVolumes.registerGlueVolumes(Trk::tubeOuterCover,sVols);
    volGlueVolumes.registerGlueVolumes(Trk::negativeFaceXY,sVolsNeg);
    volGlueVolumes.registerGlueVolumes(Trk::positiveFaceXY,sVolsPos);

  } else {
    // enclosed muon objects ? 
    std::vector<const Trk::DetachedTrackingVolume*>* muonObjs = getDetachedObjects( vol);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    muonObjs,
				    volumeName);
    // prepare blending
    if (m_blendInertMaterial && muonObjs) {
      for (unsigned int id=0;id<muonObjs->size();id++) {
	if (!(*muonObjs)[id]->layerRepresentation() && (*muonObjs)[id]->constituents()){
	  if (!m_blendMap[(*muonObjs)[id]]) 
	    m_blendMap[(*muonObjs)[id]] = new std::vector<const Trk::TrackingVolume*>;
	  m_blendMap[(*muonObjs)[id]]->push_back(tVol);
	}
      }
    }  
  }

  return tVol;
} 

const Trk::TrackingVolume* Muon::MuonTrackingGeometryBuilder::processVolume(const Trk::Volume* vol, int mode , std::string volumeName) const
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << name() << "processing volume in mode:"<< mode << endreq;

  // mode : -1 ( adjusted z/phi partition )
  //         0 ( -"- plus barrel H binning )            
  //         0 ( -"- plus inner endcap H binning )            
  //         0 ( -"- plus outer endcap H binning )            

  const Trk::TrackingVolume* tVol = 0;

  unsigned int colorCode = 12;

  // retrieve cylinder
  const Trk::CylinderVolumeBounds* cyl=dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  if (!cyl) {
    log << MSG::ERROR << " process volume: volume cylinder boundaries not retrieved, return 0 " << endreq;
    return 0; 
  }
  // create vector of zSteps for this volume
  std::vector<double> zSteps;
  std::vector<int> zTypes;
  zSteps.clear();
  double zPos = vol->center()[2];
  double hz = cyl->halflengthZ();
  double z1 = zPos-hz;  double z2 = zPos+hz ;
  zSteps.push_back(z1);
  for (unsigned int iz=0;iz<m_zPartitions.size();iz++) {
    if ( m_zPartitions[iz]> z1 && m_zPartitions[iz] < z2 ) {
      zSteps.push_back(m_zPartitions[iz]);
      if (!zTypes.size()) { 
	if (iz==0) zTypes.push_back(0);
        else zTypes.push_back(m_zPartitionsType[iz-1]);
      }
      zTypes.push_back(m_zPartitionsType[iz]);
      z1 = m_zPartitions[iz];
    }
  }
  zSteps.push_back(z2);

  // R/H binning ?
  unsigned int etaN = zSteps.size()-1;
  unsigned int phiN = m_adjustedPhi.size();   

  if ( mode > -1 ) {
    // create z,phi bin utilities
    Trk::BinUtility1DZZ* zBinUtil = new Trk::BinUtility1DZZ(zSteps);
    Trk::BinUtility1DF* pBinUtil = new Trk::BinUtility1DF(m_adjustedPhi);
    std::vector<std::vector<Trk::BinUtility1D*> >  hBinUtil;
    for (unsigned iz=0;iz < zSteps.size()-1; iz++) {
      std::vector<Trk::BinUtility1D*> phBinUtil; 
      for (unsigned ip=0;ip < m_adjustedPhi.size(); ip++) {
        // retrieve reference phi
	double phiRef = 0.5*m_adjustedPhi[ip];
        if (ip<m_adjustedPhi.size()-1) phiRef += 0.5*m_adjustedPhi[ip+1] ;
	else phiRef += 0.5*m_adjustedPhi[0]+M_PI ;

        phBinUtil.push_back(new Trk::BinUtility1DH(phiRef,m_hPartitions[mode][zTypes[iz]][m_adjustedPhiType[ip]]));
      }
      hBinUtil.push_back(phBinUtil);
    }

    // subvolume boundaries
    Trk::BevelledCylinderVolumeBounds* subBds=0;

    // create subvolumes & BinnedArray
    std::vector<Trk::TrackingVolumeOrderPosition>  subVolumesVect;
    std::vector<std::vector<std::vector<const Trk::TrackingVolume*> > > subVolumes;
    std::vector<std::vector<Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> > > > hBins;
    std::vector<const Trk::TrackingVolume*> sVolsInn;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsOut;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsNeg;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsPos;             // for gluing
    for (unsigned int eta = 0; eta < zSteps.size()-1; eta++) {
      colorCode = 26 -colorCode;
      double posZ = 0.5*(zSteps[eta] + zSteps[eta+1]) ;
      double   hZ = 0.5*fabs(zSteps[eta+1] - zSteps[eta]) ;
      std::vector<std::vector<const Trk::TrackingVolume*> > phiSubs;
      std::vector<Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> > >  phBins;
      for (unsigned int phi = 0; phi < phiN; phi++) {
	colorCode = 26 -colorCode;
	double posPhi = 0.5*m_adjustedPhi[phi];
	double phiSect = 0.;
        if (phi<phiN-1) {
	  posPhi += 0.5*m_adjustedPhi[phi+1] ;
          phiSect = 0.5*fabs(m_adjustedPhi[phi+1]-m_adjustedPhi[phi]);
	} else {
	  posPhi += 0.5*m_adjustedPhi[0]+M_PI ;
          phiSect = 0.5*fabs(m_adjustedPhi[0]+2*M_PI-m_adjustedPhi[phi]);
	}
	std::vector<std::pair<int,double> > hSteps =  m_hPartitions[mode][zTypes[eta]][m_adjustedPhiType[phi]];
	std::vector<const Trk::TrackingVolume*> hSubs;
	std::vector<Trk::TrackingVolumeOrderPosition> hSubsTr;
        unsigned int hCode = 1; 
        for (unsigned int h = 0; h < hSteps.size()-1; h++) {
	  hCode = 1 - hCode; 
          int volType = 0;     // cylinder 
          if ( hSteps[h].first == 1 && hSteps[h+1].first == 0 ) volType = 1;  
          if ( hSteps[h].first == 0 && hSteps[h+1].first == 1 ) volType = 2;  
          if ( hSteps[h].first == 1 && hSteps[h+1].first == 1 ) volType = 3;  
	  // define subvolume
	  subBds = new Trk::BevelledCylinderVolumeBounds(hSteps[h].second,
							 hSteps[h+1].second,
							 phiSect, 
							 hZ,
							 volType);
	  HepTransform3D* transf = new HepTransform3D(HepRotateZ3D(posPhi)*HepTranslateZ3D(posZ));
	  Trk::Volume subVol(transf, subBds);
       
	  // enclosed muon objects ? also adjusts material properties in case of material blend  
	  std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) +MuonGM::buildString(h,2) ; 
	  std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( &subVol);

	  const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( subVol,
								     m_muonMaterial,
								     m_muonMagneticField,
								     detVols,
								     volName );
                                                                                                                          
	  // prepare blending
	  if (m_blendInertMaterial && detVols) {
	    for (unsigned int id=0;id<detVols->size();id++) {
	      if (!(*detVols)[id]->layerRepresentation() && (*detVols)[id]->constituents()){
		if (!m_blendMap[(*detVols)[id]]) 
		  m_blendMap[(*detVols)[id]] = new std::vector<const Trk::TrackingVolume*>;
		m_blendMap[(*detVols)[id]]->push_back(sVol);
	      }
	    }
	  }  
          //
          sVol->registerColorCode(colorCode+hCode);
	  // reference position 
	  HepPoint3D gp(subBds->mediumRadius(),0.,0.);
	  subVolumesVect.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, false),
	                                                       new HepPoint3D((*transf)*gp)));
	  hSubsTr.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
	                                                       new HepPoint3D((*transf)*gp)));
	  hSubs.push_back(sVol);

	  //glue subVolume
	  if (h==0)                sVolsInn.push_back(sVol); 
	  if (h==hSteps.size()-1)  sVolsOut.push_back(sVol); 
	  if (eta==0)      sVolsNeg.push_back(sVol); 
	  if (eta==etaN-1) sVolsPos.push_back(sVol); 
          // in R/H
          if (h>0) { // glue 'manually'
            if (volType == 1 || volType == 3 ) {  // plane surface
	      m_trackingVolumeHelper->setOutsideTrackingVolume(*sVol, Trk::tubeSectorInnerCover,hSubs[h-1]); 
	      m_trackingVolumeHelper->setOutsideTrackingVolume(*(hSubs[h-1]), Trk::tubeSectorOuterCover,sVol); 
            } else {  // cylinder surface
	      m_trackingVolumeHelper->setInsideTrackingVolume(*sVol, Trk::tubeSectorInnerCover,hSubs[h-1]); 
	      m_trackingVolumeHelper->setOutsideTrackingVolume(*(hSubs[h-1]), Trk::tubeSectorOuterCover,sVol);
	    } 
          }
	  // in phi        
	  if ( phiN>1 && phi>0) {
	    m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*sVol, Trk::tubeSectorNegativePhi,phBins[phi-1]);
	    if ( phi==phiN-1 )  m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*sVol, Trk::tubeSectorPositivePhi, phBins[0]);
	  }
	  // in eta
	  if ( etaN>1 && eta>0) m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*sVol, Trk::negativeFaceXY, hBins[eta-1][phi]);     
	  //
	}
        phiSubs.push_back(hSubs); 
	Trk::BinnedArray1D<Trk::TrackingVolume>* volBinArray = new Trk::BinnedArray1D<Trk::TrackingVolume>(hSubsTr,hBinUtil[eta][phi]->clone());
        phBins.push_back(Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> >(volBinArray));

        // finish phi gluing
        if (phiN>1 && phi>0) {
	  for (unsigned int j=0; j<phiSubs[phi-1].size(); j++) {
	    m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*(phiSubs[phi-1][j]), Trk::tubeSectorPositivePhi,phBins[phi]);
	  }
	}
        if (phi==phiN-1) {
	  for (unsigned int j=0; j<phiSubs[0].size(); j++) {
	    m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*(phiSubs[0][j]), Trk::tubeSectorNegativePhi,phBins[phi]);
	  }
	}
        // finish eta gluing
        if ( etaN>1 && eta>0) {
	  for (unsigned int j=0; j<subVolumes[eta-1][phi].size(); j++) {
	    m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*(subVolumes[eta-1][phi][j]), Trk::positiveFaceXY,phBins[phi]);
	  }
	}
      }
      subVolumes.push_back(phiSubs);
      hBins.push_back(phBins);
    }

    Trk::BinUtility3DZFH* volBinUtil=new Trk::BinUtility3DZFH(zBinUtil,pBinUtil,hBinUtil,new HepTransform3D(vol->transform()));

    Trk::BinnedArray3D<Trk::TrackingVolume>* subVols=new Trk::BinnedArray3D<Trk::TrackingVolume>(subVolumesVect,volBinUtil);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    0,subVols,
				    volumeName);
    // register glue volumes
    const Trk::GlueVolumesDescriptor& volGlueVolumes = tVol->glueVolumesDescriptor();
    volGlueVolumes.registerGlueVolumes(Trk::tubeInnerCover,sVolsInn);
    volGlueVolumes.registerGlueVolumes(Trk::tubeOuterCover,sVolsOut);
    volGlueVolumes.registerGlueVolumes(Trk::negativeFaceXY,sVolsNeg);
    volGlueVolumes.registerGlueVolumes(Trk::positiveFaceXY,sVolsPos);
  
    return tVol;
  } 

  // proceed with 2D z/phi binning
  // partitions ? include protection against wrong setup
  if (phiN < 1) {
    log << MSG::ERROR << name() << "wrong partition setup" << endreq;
    phiN = 1;
  } else {
    log << MSG::DEBUG << name() << "partition setup:(z,phi):"<<etaN<<","<<phiN << endreq;
  }

  if ( etaN * phiN > 1 ) {  // partition
    // subvolume boundaries
    Trk::CylinderVolumeBounds* subBds=0;

    // create subvolumes & BinnedArray
    std::vector<Trk::TrackingVolumeOrderPosition> subVolumes(etaN*phiN);
    std::vector<const Trk::TrackingVolume*> sVols(etaN*phiN);                // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsNeg(phiN);             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsPos(phiN);             // for gluing
    for (unsigned int eta = 0; eta < zSteps.size()-1; eta++) {
      double posZ = 0.5*(zSteps[eta] + zSteps[eta+1]) ;
      double   hZ = 0.5*fabs(zSteps[eta+1] - zSteps[eta]) ;
      colorCode = 26 -colorCode;
      for (unsigned int phi = 0; phi < phiN; phi++) {
	colorCode = 26 -colorCode;
	double posPhi = 0.5*m_adjustedPhi[phi];
	double phiSect = 0.;
        if (phi<phiN-1) {
	  posPhi += 0.5*m_adjustedPhi[phi+1] ;
          phiSect = 0.5*fabs(m_adjustedPhi[phi+1]-m_adjustedPhi[phi]);
	} else {
	  posPhi += 0.5*m_adjustedPhi[0]+M_PI ;
          phiSect = 0.5*fabs(m_adjustedPhi[0]+2*M_PI-m_adjustedPhi[phi]);
	}
	// define subvolume
	subBds = new Trk::CylinderVolumeBounds(cyl->innerRadius(),
					       cyl->outerRadius(),
					       phiSect, 
					       hZ);
        HepTransform3D* transf = new HepTransform3D(HepRotateZ3D(posPhi)*HepTranslateZ3D(posZ));
        Trk::Volume subVol(transf, subBds);     
        // enclosed muon objects ?   
	std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) ; 

	Trk::MaterialProperties mat=m_muonMaterial;
        std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( &subVol);
        const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( subVol,
								   m_muonMaterial,
								   m_muonMagneticField,
								   detVols,
								   volName );
        // prepare blending
        if (m_blendInertMaterial && detVols) {
          for (unsigned int id=0;id<detVols->size();id++) {
            if (!(*detVols)[id]->layerRepresentation() && (*detVols)[id]->constituents()){
	      if (!m_blendMap[(*detVols)[id]]) 
		m_blendMap[(*detVols)[id]] = new std::vector<const Trk::TrackingVolume*>;
	      m_blendMap[(*detVols)[id]]->push_back(sVol);
	    }
	  }
	}  
        //delete subVol;
        sVol->registerColorCode(colorCode); 
        // reference position 
	HepPoint3D gp(subBds->outerRadius(),0.,0.);
	//subVolumes.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
        //                                                     new HepPoint3D((*transf)*gp)));
	subVolumes[phi*etaN+eta] = Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, false),
                                                             new HepPoint3D((*transf)*gp));
        //glue subVolumes
        //sVols[phi*etaN+eta] = sVol; 
        sVols[phiN*eta+phi] = sVol; 
        if (eta==0)      sVolsNeg[phi]=sVol; 
        if (eta==etaN-1) sVolsPos[phi]=sVol; 
        // in phi        
        if ( phiN>1 && phi>0) {
          m_trackingVolumeHelper->glueTrackingVolumes(*sVol, Trk::tubeSectorNegativePhi,
						      *(sVols[eta*phiN+phi-1]), Trk::tubeSectorPositivePhi);
          if ( phi==phiN-1 )  m_trackingVolumeHelper->glueTrackingVolumes(*(sVols[eta*phiN]), Trk::tubeSectorNegativePhi,
                                                                          *sVol, Trk::tubeSectorPositivePhi);
	}
        // in eta
        if ( etaN>1 && eta>0) m_trackingVolumeHelper->glueTrackingVolumes(*sVol, Trk::negativeFaceXY,
						      *(sVols[(eta-1)*phiN+phi]), Trk::positiveFaceXY);        
        //
      }
    }

    Trk::BinUtility2DZF* volBinUtil=new Trk::BinUtility2DZF(zSteps,m_adjustedPhi,new HepTransform3D(vol->transform()));

    Trk::BinnedArray2D<Trk::TrackingVolume>* subVols=new Trk::BinnedArray2D<Trk::TrackingVolume>(subVolumes,volBinUtil);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    0,subVols,
				    volumeName);
    // register glue volumes
    const Trk::GlueVolumesDescriptor& volGlueVolumes = tVol->glueVolumesDescriptor();
    volGlueVolumes.registerGlueVolumes(Trk::tubeInnerCover,sVols);
    volGlueVolumes.registerGlueVolumes(Trk::tubeOuterCover,sVols);
    volGlueVolumes.registerGlueVolumes(Trk::negativeFaceXY,sVolsNeg);
    volGlueVolumes.registerGlueVolumes(Trk::positiveFaceXY,sVolsPos);

  } else {
    // enclosed muon objects ? 
    std::vector<const Trk::DetachedTrackingVolume*>* muonObjs = getDetachedObjects( vol);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    muonObjs,
				    volumeName);
    // prepare blending
    if (m_blendInertMaterial && muonObjs) {
      for (unsigned int id=0;id<muonObjs->size();id++) {
	if (!(*muonObjs)[id]->layerRepresentation() && (*muonObjs)[id]->constituents()){
	  if (!m_blendMap[(*muonObjs)[id]]) 
	    m_blendMap[(*muonObjs)[id]] = new std::vector<const Trk::TrackingVolume*>;
	  m_blendMap[(*muonObjs)[id]]->push_back(tVol);
	}
      }
    }  
  }

  return tVol;
} 

std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonTrackingGeometryBuilder::getDetachedObjects(const Trk::Volume* vol) const
{

  std::vector<const Trk::DetachedTrackingVolume*>* detTVs = 0;

  if (!m_stations && !m_inertObjs) return detTVs;
  
  // get min/max Z/Phi from volume (allways a cylinder/bevelled cylinder volume )  
  const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::BevelledCylinderVolumeBounds* bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&(vol->volumeBounds()));
   
  double rmed = 0.; double dphi = 0.; double hz = 0.; double rMin = 0.; double rMax = 0.; double rMaxc = 0.; int type = 0;
  if (cyl) {
    rmed = cyl->mediumRadius();
    dphi = cyl->halfPhiSector();
    hz = cyl->halflengthZ();
    rMin = cyl->innerRadius();
    rMax = cyl->outerRadius();
    rMaxc = rMax;
  } else if (bcyl) {
    rmed = bcyl->mediumRadius();
    dphi = bcyl->halfPhiSector();
    hz = bcyl->halflengthZ();
    rMin = bcyl->innerRadius();
    rMax = bcyl->outerRadius();    
    rMaxc = rMax;    
    type = bcyl->type();
    if (type>1) rMaxc *=1./cos(dphi);
  } else return 0;
 
  HepPoint3D center(rmed,0.,0.);
  center = vol->transform() * center;

  double zMin = center[2] - hz; 
  double zMax = center[2] + hz;
  double pMin = 0.;
  double pMax = +2*M_PI;
  bool phiLim = false;
  if (dphi < M_PI) { 
    pMin = center.phi() - dphi + M_PI; 
    pMax = center.phi() + dphi + M_PI;
    phiLim = true;
  } 
   
  std::list<const Trk::DetachedTrackingVolume*> detached;
  // active, use corrected rMax
  if (m_stationSpan) {
    for (unsigned int i=0; i<m_stationSpan->size() ; i++) {
      const Muon::Span* s = (*m_stationSpan)[i];
      bool rLimit = !m_static3d || ( (*s)[4] <= rMaxc && (*s)[5] >= rMin );
   
      if ( rLimit && (*s)[0] <= zMax && (*s)[1] >= zMin && m_stations->size()>i) {
	if (phiLim) {
	  if (pMin>=0 && pMax<=2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && (*s)[2] <= pMax && (*s)[3] >= pMin ) detached.push_back((*m_stations)[i]);
	    if ( (*s)[2]>(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin) ) detached.push_back((*m_stations)[i]);
	  } else if (pMin < 0) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin+2*M_PI) ) detached.push_back((*m_stations)[i]);
	    if ( (*s)[2]>(*s)[3]  ) detached.push_back((*m_stations)[i]);
	  } else if (pMax > 2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax-2*M_PI || (*s)[3] >= pMin) ) detached.push_back((*m_stations)[i]);
	    if ( (*s)[2]>(*s)[3]  ) detached.push_back((*m_stations)[i]);
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
      const Muon::Span* s = (*m_inertSpan)[i];
      //bool rail = ( (*m_inertObjs)[i]->name() == "Rail" ) ? true : false; 
      bool rLimit = (!m_static3d || ( (*s)[4] <= rMax && (*s)[5] >= rMin ) ); 
      if ( rLimit && (*s)[0] <= zMax && (*s)[1] >= zMin && m_inertObjs->size()>i) {
	if (phiLim) {
	  if (pMin>=0 && pMax<=2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && (*s)[2] <= pMax && (*s)[3] >= pMin )  detached.push_back((*m_inertObjs)[i]);
	    if ( (*s)[2]>(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin) ) detached.push_back((*m_inertObjs)[i]);
	  } else if (pMin < 0) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin+2*M_PI) ) detached.push_back((*m_inertObjs)[i]);
	    if ( (*s)[2]>(*s)[3]  ) detached.push_back((*m_inertObjs)[i]);
	  } else if (pMax > 2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax-2*M_PI || (*s)[3] >= pMin) ) detached.push_back((*m_inertObjs)[i]);
	    if ( (*s)[2]>(*s)[3]  ) detached.push_back((*m_inertObjs)[i]);
	  }
	} else  detached.push_back((*m_inertObjs)[i]);
      } 
    }
  }
  if (!detached.empty()) detTVs = new std::vector<const Trk::DetachedTrackingVolume*>(detached.begin(), detached.end()); 
  return detTVs;
}

bool Muon::MuonTrackingGeometryBuilder::enclosed(const Trk::Volume* vol, const Trk::Volume* cs) const
{
  MsgStream log(msgSvc(), name());

  bool encl = false;
  
  // get min/max Z/Phi from volume (allways a cylinder/bevelled cylinder volume )  
  const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::BevelledCylinderVolumeBounds* bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&(vol->volumeBounds()));
   
  double rmed = 0.; double dphi = 0.; double hz = 0.; double rMin = 0.; double rMax = 0.; double rMaxc = 0.; int type = 0;
  if (cyl) {
    rmed = cyl->mediumRadius();
    dphi = cyl->halfPhiSector();
    hz = cyl->halflengthZ();
    rMin = cyl->innerRadius();
    rMax = cyl->outerRadius();
    rMaxc = rMax;
  } else if (bcyl) {
    rmed = bcyl->mediumRadius();
    dphi = bcyl->halfPhiSector();
    hz = bcyl->halflengthZ();
    rMin = bcyl->innerRadius();
    rMax = bcyl->outerRadius();    
    rMaxc = rMax;    
    type = bcyl->type();
    if (type>1) rMaxc *=1./cos(dphi);
  } else return 0;
 
  HepPoint3D center(rmed,0.,0.);
  center = vol->transform() * center;

  double zMin = center[2] - hz; 
  double zMax = center[2] + hz;
  double pMin = 0.;
  double pMax = +2*M_PI;
  bool phiLim = false;
  if (dphi < M_PI) { 
    pMin = center.phi() - dphi + M_PI; 
    pMax = center.phi() + dphi + M_PI;
    phiLim = true;
  } 
  //
  //const Muon::Span* s = findVolumeSpan(&(cs->volumeBounds()), cs->transform(), 0.,0.) ;
  std::auto_ptr<const Muon::Span> s (findVolumeSpan(&(cs->volumeBounds()), cs->transform(), 0.,0.) );
  //log << MSG::DEBUG << "enclosing volume:z:"<< zMin<<","<<zMax<<":r:"<< rMin<<","<<rMax<<":phi:"<<pMin<<","<<pMax<< endreq;
  //log << MSG::DEBUG << "constituent:z:"<< (*s)[0]<<","<<(*s)[1]<<":r:"<< (*s)[4]<<","<<(*s)[5]<<":phi:"<<(*s)[2]<<","<<(*s)[3]<< endreq;
  //
  bool rLimit = (!m_static3d || ( (*s)[4] <= rMax && (*s)[5] >= rMin ) ); 
  if ( rLimit && (*s)[0] <= zMax && (*s)[1] >= zMin ) {
    if (phiLim) {
      if (pMin>=0 && pMax<=2*M_PI) {
	if ( (*s)[2]<=(*s)[3] && (*s)[2] <= pMax && (*s)[3] >= pMin )  return true;
	if ( (*s)[2]>(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin) ) return true;
      } else if (pMin < 0) {
	if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin+2*M_PI) ) return true;
	if ( (*s)[2]>(*s)[3]  ) return true;
      } else if (pMax > 2*M_PI) {
	if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax-2*M_PI || (*s)[3] >= pMin) ) return true;
	if ( (*s)[2]>(*s)[3]  ) return true;
      }
    } else {
      return true;
    }
  }
  return encl;
}

void Muon::MuonTrackingGeometryBuilder::checkVolume(const Trk::TrackingVolume* vol ) const
{
  MsgStream log(msgSvc(), name());

  std::cout << "MuonTrackingGeometryBuilder::checkVolume: " << vol->volumeName() << std::endl;

  const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  if (!cyl) {
    std::cout << "MuonTrackingGeometryBuilder::checkVolume:not a cylinder, return" << std::endl;
    return;
  }
  std::cout << "cylinder dimensions:innerR,outerR,halfZ,halfPhi:" << cyl->innerRadius() << "," << cyl->outerRadius()
	    << cyl->halflengthZ() << "," << cyl->halfPhiSector() << std::endl;

  const Trk::GlueVolumesDescriptor& glueDescr = vol->glueVolumesDescriptor();
  std::cout << "glue volumes:neg,pos,inn,outer" << (glueDescr.glueVolumes(Trk::negativeFaceXY)).size()<<"," 
                                                << (glueDescr.glueVolumes(Trk::positiveFaceXY)).size()<<","
                                                << (glueDescr.glueVolumes(Trk::tubeInnerCover)).size()<<","
	                                        << (glueDescr.glueVolumes(Trk::tubeOuterCover)).size()<<std::endl;
}

void Muon::MuonTrackingGeometryBuilder::getZParts() const
{
  // hardcode for the moment
  m_zPartitions.resize(40);
  m_zPartitionsType.resize(40);
  
  m_zPartitions[0]=-22100.;  m_zPartitions[1]=-21630.;  m_zPartitions[2]=-21200.;  m_zPartitions[3]=-15360.;
  m_zPartitions[4]=-15150.;  m_zPartitions[5]=-14940.;  m_zPartitions[6]=-14730.;  m_zPartitions[7]=-14520.;
  m_zPartitions[8]=-14080.;  m_zPartitions[9]=-13650.;  m_zPartitions[10]=-13440.;  m_zPartitions[11]=-13230.;
  m_zPartitions[12]=-11211.;  m_zPartitions[13]=-10479.;  m_zPartitions[14]=-8611.;  m_zPartitions[15]=-7879.;
  m_zPartitions[16]=-5503.;  m_zPartitions[17]=-4772.;  m_zPartitions[18]=-2078.;  m_zPartitions[19]=-1347.;

  for (unsigned int i = 20; i<40 ; i++) m_zPartitions[i] = - m_zPartitions[39-i];  

  for (unsigned int i = 0; i<40 ; i++) m_zPartitionsType[i] = 0;
  m_zPartitionsType[12] = 1;  m_zPartitionsType[14] = 1;   m_zPartitionsType[16] = 1; m_zPartitionsType[18] = 1;
  m_zPartitionsType[20] = 1;  m_zPartitionsType[22] = 1;   m_zPartitionsType[24] = 1; m_zPartitionsType[26] = 1;
  
  return;
}

void Muon::MuonTrackingGeometryBuilder::getPhiParts() const
{
  // hardcode for the moment
  m_adjustedPhi.resize(16);
  m_adjustedPhiType.resize(16);

  double phiSect[2];
  phiSect[0] = ( M_PI/8 - 0.105 ); 
  phiSect[1] = 0.105 ; 

  m_adjustedPhi[0]= - phiSect[0];
  m_adjustedPhiType[0]= 0;
  int ic = 0; int is = 1;

  while (ic < 15 ) {
    ic++; is = 1 - is;
    m_adjustedPhi[ic]= m_adjustedPhi[ic-1] + 2*phiSect[is] ;
    m_adjustedPhiType[ic]= 1-is;
  }

  return;
}

void Muon::MuonTrackingGeometryBuilder::getHParts() const
{
  // hardcode for the moment
  m_hPartitions.clear();              // barrel, inner endcap, outer endcap

  // barrel 2x2
  std::vector<std::pair<int,double> >  barrelZ0F0;
  barrelZ0F0.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  barrelZ0F0.push_back( std::pair<int,double>(0,4450.) );                // for DiskShieldingBackDisk
  barrelZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  barrelZ0F1;
  barrelZ0F1.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  barrelZ0F1.push_back( std::pair<int,double>(1,4700.) );
  barrelZ0F1.push_back( std::pair<int,double>(1,5900.) );
  barrelZ0F1.push_back( std::pair<int,double>(1,8900.) );
  barrelZ0F1.push_back( std::pair<int,double>(1,10100.) );
  barrelZ0F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  barrelZ1F0;
  barrelZ1F0.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  barrelZ1F0.push_back( std::pair<int,double>(1,5800.) );
  barrelZ1F0.push_back( std::pair<int,double>(1,6500.) );
  barrelZ1F0.push_back( std::pair<int,double>(1,8400.) );
  barrelZ1F0.push_back( std::pair<int,double>(1,9100.) );
  barrelZ1F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  barrelZ1F1;
  barrelZ1F1.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  barrelZ1F1.push_back( std::pair<int,double>(1,4700.) );
  barrelZ1F1.push_back( std::pair<int,double>(1,6000.) );
  barrelZ1F1.push_back( std::pair<int,double>(1,6800.) );
  barrelZ1F1.push_back( std::pair<int,double>(1,8900.) );
  barrelZ1F1.push_back( std::pair<int,double>(1,10100.) );
  barrelZ1F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  barrelZF(2);
  barrelZF[0].push_back(barrelZ0F0);
  barrelZF[0].push_back(barrelZ0F1);
  barrelZF[1].push_back(barrelZ1F0);
  barrelZF[1].push_back(barrelZ1F1);

  //inner endcap 2x2
  std::vector<std::pair<int,double> >  innerZ0F0;
  innerZ0F0.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  innerZ0F0.push_back( std::pair<int,double>(1,5300.) );
  innerZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  innerZ0F1;
  innerZ0F1.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  innerZ0F1.push_back( std::pair<int,double>(1,4700.) );
  innerZ0F1.push_back( std::pair<int,double>(1,5900.) );
  innerZ0F1.push_back( std::pair<int,double>(1,8900.) );
  innerZ0F1.push_back( std::pair<int,double>(1,10100.) );
  innerZ0F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  innerZ1F0;
  innerZ1F0.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  innerZ1F0.push_back( std::pair<int,double>(1,5300.) );
  innerZ1F0.push_back( std::pair<int,double>(1,5800.) );
  innerZ1F0.push_back( std::pair<int,double>(1,6500.) );
  innerZ1F0.push_back( std::pair<int,double>(1,8400.) );
  innerZ1F0.push_back( std::pair<int,double>(1,9100.) );
  innerZ1F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  innerZ1F1;
  innerZ1F1.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  innerZ1F1.push_back( std::pair<int,double>(1,4700.) );
  innerZ1F1.push_back( std::pair<int,double>(1,6000.) );
  innerZ1F1.push_back( std::pair<int,double>(1,6800.) );
  innerZ1F1.push_back( std::pair<int,double>(1,8900.) );
  innerZ1F1.push_back( std::pair<int,double>(1,10100.) );
  innerZ1F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  innerZF(2);
  innerZF[0].push_back(innerZ0F0);
  innerZF[0].push_back(innerZ0F1);
  innerZF[1].push_back(innerZ1F0);
  innerZF[1].push_back(innerZ1F1);

  // outer 1x1
  std::vector<std::pair<int,double> >  outerZ0F0;
  outerZ0F0.push_back( std::pair<int,double>(0,m_outerShieldRadius) );
  outerZ0F0.push_back( std::pair<int,double>(0,2275.) );
  outerZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  outerZF(2);
  outerZF[0].push_back(outerZ0F0);
  outerZF[0].push_back(outerZ0F0);
  outerZF[1].push_back(outerZ0F0);
  outerZF[1].push_back(outerZ0F0);

  // collect everything
  m_hPartitions.push_back(barrelZF);
  m_hPartitions.push_back(innerZF);
  m_hPartitions.push_back(outerZF);
   
  return;
}

double  Muon::MuonTrackingGeometryBuilder::calculateVolume( const Trk::Volume* envelope) const
{
  double envVol = 0.;
  
  if (!envelope) return 0.;
  
  const Trk::CylinderVolumeBounds*  cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::CuboidVolumeBounds*    box = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::BevelledCylinderVolumeBounds*  bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::PrismVolumeBounds* prism = dynamic_cast<const Trk::PrismVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::SimplePolygonBrepVolumeBounds* spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::CombinedVolumeBounds*  comb = dynamic_cast<const Trk::CombinedVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::SubtractedVolumeBounds*  sub = dynamic_cast<const Trk::SubtractedVolumeBounds*> (&(envelope->volumeBounds()));  

  if ( cyl ) envVol = 2*cyl->halfPhiSector()*(cyl->outerRadius()*cyl->outerRadius()-cyl->innerRadius()*cyl->innerRadius())*cyl->halflengthZ();
  if ( box ) envVol = (8*box->halflengthX()*box->halflengthY()*box->halflengthZ());
  if ( trd ) envVol = (4*(trd->minHalflengthX()+trd->maxHalflengthX())*trd->halflengthY()*trd->halflengthZ());
  if ( bcyl ) {
    int type = bcyl->type();
    if ( type<1 ) envVol = 2*bcyl->halfPhiSector()*(bcyl->outerRadius()*bcyl->outerRadius()-bcyl->innerRadius()*bcyl->innerRadius())*bcyl->halflengthZ(); 
    if ( type==1 ) envVol = 2*bcyl->halflengthZ()*( bcyl->halfPhiSector()*bcyl->outerRadius()*bcyl->outerRadius()
						    -bcyl->innerRadius()*bcyl->innerRadius()*tan(bcyl->halfPhiSector()) ); 
    if ( type==2 ) envVol = 2*bcyl->halflengthZ()*( -bcyl->halfPhiSector()*bcyl->innerRadius()*bcyl->innerRadius()
						    +bcyl->outerRadius()*bcyl->outerRadius()*tan(bcyl->halfPhiSector()) ); 
    if ( type==3 ) envVol = 2*bcyl->halflengthZ()*tan(bcyl->halfPhiSector())*( bcyl->outerRadius()*bcyl->outerRadius() 
									       -bcyl->innerRadius()*bcyl->innerRadius()); 
  }
  if ( prism ) {
    std::vector<std::pair<double,double> > v=prism->xyVertices();
    double a2 = v[1].first*v[1].first+v[1].second*v[1].second
               +v[0].first*v[0].first+v[0].second*v[0].second
	    -2*(v[0].first*v[1].first+v[0].second*v[1].second);
    double c2 = v[2].first*v[2].first+v[2].second*v[2].second
               +v[0].first*v[0].first+v[0].second*v[0].second
	    -2*(v[0].first*v[2].first+v[0].second*v[2].second);
    double ca = v[1].first*v[2].first+v[1].second*v[2].second
               +v[0].first*v[0].first+v[0].second*v[0].second
	       -v[0].first*v[1].first-v[0].second*v[1].second
	       -v[0].first*v[2].first-v[0].second*v[2].second;
    double vv = sqrt(c2-ca*ca/a2);
    envVol = vv*sqrt(a2)*prism->halflengthZ();
  }
  if ( spb ) {
    envVol = calculateVolume(spb->combinedVolume());    // exceptional use of combined volume (no intersections)
  }   
  if ( comb ) {
    envVol = calculateVolume(comb->first()) + calculateVolume(comb->second());
  }
  if ( sub ) {
    return -1;
  }

  return envVol;
}

void Muon::MuonTrackingGeometryBuilder::blendMaterial() const
{
  MsgStream log(msgSvc(), name());
  // loop over map
  std::map<const Trk::DetachedTrackingVolume*,std::vector<const Trk::TrackingVolume*>* >::iterator mIter = m_blendMap.begin();

  std::vector<std::pair<const Trk::Volume*,std::pair<double,double> > >* cs = 0;
  
  for ( ; mIter!= m_blendMap.end(); mIter++) {
    cs = (*mIter).first->constituents();     
    if (!cs) continue;
    // find material source
    const Trk::MaterialProperties* detMat = (*mIter).first->trackingVolume();
    //std::cout << "blending:"<<(*mIter).first->name()<< std::endl;
    if ( (*mIter).first->trackingVolume()->confinedDenseVolumes()) detMat = (*(*mIter).first->trackingVolume()->confinedDenseVolumes())[0];
    for (unsigned int ic=0; ic<cs->size(); ic++) {
      const Trk::Volume* nCs = new Trk::Volume(*((*cs)[ic].first),(*mIter).first->trackingVolume()->transform());
      double csVol = (*cs)[ic].second.first*calculateVolume(nCs);      
      double enVol = 0.;
      // loop over frame volumes, check if confined
      std::vector<const Trk::TrackingVolume*>::iterator fIter = (*mIter).second->begin(); 
      std::vector<bool> fEncl; 
      fEncl.clear();
      // blending factors can be saved, and not recalculated for each clone
      for ( ; fIter!=(*mIter).second->end(); fIter++) {
        fEncl.push_back(enclosed(*fIter,nCs));
        if ( fEncl.back() ) enVol += calculateVolume(*fIter);
        //if ((*cs)[ic].second.second<0. && fEncl.back() ) enVol += calculateVolume(*fIter);
        //if ((*cs)[ic].second.second<0.)  enVol += calculateVolume(*fIter);
      }
      delete nCs;
      // diluting factor
      double dil =  enVol>0. ?  csVol/enVol : 0.;
      //std::cout << "const:dil:"<< ic<<","<<dil<< std::endl;
      if (dil>0.) { 
	for ( fIter=(*mIter).second->begin(); fIter!=(*mIter).second->end(); fIter++) { 
	  if (fEncl[fIter-(*mIter).second->begin()]) (*fIter)->addMaterial(*detMat,dil); ;
	}
      } else {
	log << MSG::DEBUG << "diluting factor:"<< dil<<" for "<< (*mIter).first->name()<<","<<ic<<endreq;
      }
    }
  }
}

