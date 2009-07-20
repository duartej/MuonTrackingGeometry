//////////////////////////////////////////////////////////////////
// MuonTrackingGeometryBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonTrackingGeometryBuilder.h"
#include "MuonReadoutGeometry/GlobalUtilities.h" 
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
  m_magFieldMode(Trk::RealisticField),
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
  m_outerEndcapZ(23000.),
  m_bigWheel(15600.),
  m_outerWheel(21000.),
  m_ectZ(7880.),
  m_beamPipeRadius(70.),
  //m_innerShieldRadius(810.),
  m_innerShieldRadius(850.),
  m_outerShieldRadius(1550.),
  m_diskShieldZ(6909.),
  m_barrelEtaPartition(9),
  m_innerEndcapEtaPartition(3),
  m_outerEndcapEtaPartition(3),
  m_phiPartition(16),
  m_adjustStatic(true),
  m_static3d(true),
  m_blendInertMaterial(true),
  m_removeBlended(false),
  m_inertPerm(0),
  m_alignTolerance(0.),
  m_colorCode(0),
  m_activeAdjustLevel(2),
  m_inertAdjustLevel(1),
  m_frameNum(0),
  m_frameStat(0),
  m_entryVolume("MuonSpectrometerEntrance"),
  m_exitVolume("All::Container::CompleteDetector"),
  m_chronoStatSvc( "ChronoStatSvc", n )
{
  m_stationSpan = 0;
  m_inertSpan = 0;
  m_stations = 0;
  m_inertObjs = 0;
  m_blendMap.clear();

  declareInterface<Trk::IGeometryBuilder>(this);

  declareProperty("SimpleMuonGeometry",               m_muonSimple);  
  declareProperty("LoadMSEntry",                      m_loadMSentry);  
  declareProperty("BuildActiveMaterial",              m_muonActive);  
  declareProperty("BuildInertMaterial",               m_muonInert);
  declareProperty("InertMaterialBuilder",             m_inertBuilder);
  declareProperty("MagneticFieldMode",                m_magFieldMode);
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
  declareProperty("ActiveAdjustLevel",              m_activeAdjustLevel);
  declareProperty("InertAdjustLevel",               m_inertAdjustLevel);
  declareProperty("BlendInertMaterial",             m_blendInertMaterial);
  declareProperty("RemoveBlendedMaterialObjects",   m_removeBlended);
  declareProperty("AlignmentPositionTolerance",     m_alignTolerance);
  declareProperty("ColorCode",                      m_colorCode);
  // calo entry volume & exit volume
  declareProperty("EntryVolumeName",                   m_entryVolume);
  declareProperty("ExitVolumeName",                    m_exitVolume);  
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
    if (m_magFieldMode && m_magFieldTool.retrieve().isFailure())
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
    } else {
      m_activeAdjustLevel = 0;                // no active material to consider
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

    if ( !m_muonInert ) m_inertAdjustLevel = 0;

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
  if (m_inertObjs && m_blendInertMaterial && m_removeBlended) {
    while ( m_inertPerm<m_inertObjs->size() 
	    && (*m_inertObjs)[m_inertPerm]->name().substr((*m_inertObjs)[m_inertPerm]->name().size()-4,4)=="PERM") m_inertPerm++;
  }
  
  // find object's span with tolerance for the alignment 
  if (!m_stationSpan) m_stationSpan = findVolumesSpan(m_stations, 100.*m_alignTolerance, m_alignTolerance*deg);
  if (!m_inertSpan)   m_inertSpan = findVolumesSpan(m_inertObjs,0.,0.);
 
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
  
  // magnetic field steering
  if (!m_magFieldMode)
      m_muonMagneticField = Trk::MagneticFieldProperties();
  else      
      m_muonMagneticField = Trk::MagneticFieldProperties(&(*m_magFieldTool), (Trk::MagneticFieldMode)m_magFieldMode);
  
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
  const Trk::TrackingVolume* negativeMuonOuterWheel = 0;
  Trk::CylinderVolumeBounds* negativeOuterWheelBounds = 0;
  const Trk::TrackingVolume* negativeMuonBigWheel = 0;
  Trk::CylinderVolumeBounds* negativeBigWheelBounds = 0;
  const Trk::TrackingVolume* negativeMuonOuterBuffer = 0;
  Trk::CylinderVolumeBounds* negativeOuterBufferBounds = 0;
  const Trk::TrackingVolume* positiveMuonOuterWheel = 0;
  const Trk::TrackingVolume* negativeMuonSmallWheel = 0;
  const Trk::TrackingVolume* positiveMuonSmallWheel = 0;
  Trk::CylinderVolumeBounds* negativeSmallWheelBounds = 0;
  const Trk::TrackingVolume* negativeECT = 0;
  const Trk::TrackingVolume* positiveECT = 0;
  Trk::CylinderVolumeBounds* negativeECTBounds = 0;
  const Trk::TrackingVolume* positiveMuonBigWheel = 0;
  const Trk::TrackingVolume* positiveMuonOuterBuffer = 0;
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
    if ( tvol->volumeName() == m_entryVolume ) msEntryDefined = true; 
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
	std::string nameEncl = msEntryDefined ? "All::Gaps::Barrel" : m_entryVolume ;
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
				       m_entryVolume, bcorr);
    enclosed->registerColorCode(0); 
  }
  
  // create central volume ("enclosed" + disk shields ) - this is to allow safe gluing with 3D MS binning
  getShieldParts();
  
  negDiskShieldBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
						      0.5*(m_diskShieldZ - m_barrelZ) );
  Hep3Vector negDiskShieldPosition(0.,0.,-0.5*(m_diskShieldZ+m_barrelZ));
  Trk::Volume negDiskVol(new HepTransform3D(Trk::s_idRotation,negDiskShieldPosition),negDiskShieldBounds);
  negDiskShield = processShield(&negDiskVol,2,"Muons::Detectors::NegativeDiskShield");

  posDiskShieldBounds = new Trk::CylinderVolumeBounds(m_innerBarrelRadius,
						      0.5*(m_diskShieldZ - m_barrelZ) );
  Hep3Vector posDiskShieldPosition(0.,0., 0.5*(m_diskShieldZ+m_barrelZ));
  Trk::Volume posDiskVol(new HepTransform3D(Trk::s_idRotation,posDiskShieldPosition),posDiskShieldBounds);
  posDiskShield = processShield(&posDiskVol,2,"Muons::Detectors::PositiveDiskShield");
	   
  log << MSG::INFO << name() << "glue enclosed  + disk shields" << endreq;
  const Trk::TrackingVolume* centralP = m_trackingVolumeHelper->glueTrackingVolumeArrays(*enclosed, Trk::positiveFaceXY,
											 *posDiskShield,Trk::negativeFaceXY, 
											 "Muon::Container::CentralP");    
  const Trk::TrackingVolume* central  = m_trackingVolumeHelper->glueTrackingVolumeArrays(*centralP, Trk::negativeFaceXY,
											 *negDiskShield,Trk::positiveFaceXY, 
											 "Muon::Container::Central");    
  
// define basic volumes
  if (m_adjustStatic) { getZParts(); getHParts();}

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
// build as smallWheel+ECT
// small wheel
   double smallWheelZHalfSize = 0.5*(m_ectZ - m_diskShieldZ);
   negativeSmallWheelBounds = new Trk::CylinderVolumeBounds(m_innerShieldRadius,
							    m_outerBarrelRadius,
							    smallWheelZHalfSize);
   Hep3Vector negSmallWheelPosition(0.,0.,-m_ectZ+smallWheelZHalfSize);
   HepTransform3D* negSmallWheelTransf = new HepTransform3D(Trk::s_idRotation,negSmallWheelPosition);
   Trk::Volume negSWVol(negSmallWheelTransf,negativeSmallWheelBounds);
   if (m_adjustStatic && m_static3d) negativeMuonSmallWheel = processVolume( &negSWVol,1,"Muon::Detectors::NegativeSmallWheel" ); 
   else if (m_adjustStatic) negativeMuonSmallWheel = processVolume( &negSWVol,-1,"Muon::Detectors::NegativeSmallWheel" ); 
   else negativeMuonSmallWheel = processVolume( &negSWVol,m_innerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::NegativeSmallWheel" ); 
   //checkVolume(negativeMuonSmallWheel);
   //
   Hep3Vector posSmallWheelShift(0.,0.,2*(m_ectZ-smallWheelZHalfSize));
   Trk::Volume posSWVol(negSWVol,*(new HepTransform3D(Trk::s_idRotation,posSmallWheelShift)));
   if (m_adjustStatic && m_static3d) positiveMuonSmallWheel = processVolume( &posSWVol,1,"Muon::Detectors::PositiveSmallWheel" ); 
   else if (m_adjustStatic) positiveMuonSmallWheel = processVolume( &posSWVol,-1,"Muon::Detectors::PositiveSmallWheel" ); 
   else positiveMuonSmallWheel = processVolume( &posSWVol,m_innerEndcapEtaPartition,m_phiPartition,
						 "Muon::Detectors::PositiveSmallWheel" ); 
   //checkVolume(positiveMuonSmallWheel);
// ECT
   double ectZHalfSize = 0.5*(m_innerEndcapZ - m_ectZ);
   negativeECTBounds = new Trk::CylinderVolumeBounds(m_innerShieldRadius,
						     m_outerBarrelRadius,
						     ectZHalfSize);
   Hep3Vector negECTPosition(0.,0.,-m_ectZ-ectZHalfSize);
   HepTransform3D* negECTTransf = new HepTransform3D(Trk::s_idRotation,negECTPosition);
   Trk::Volume negECTVol(negECTTransf,negativeECTBounds);
   if (m_adjustStatic && m_static3d) negativeECT = processVolume( &negECTVol,2,"Muon::Detectors::NegativeECT" ); 
   else if (m_adjustStatic) negativeECT = processVolume( &negECTVol,-1,"Muon::Detectors::NegativeECT" ); 
   else negativeECT = processVolume( &negECTVol,m_innerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::NegativeECT" ); 
   //checkVolume(negativeECT);
   //
   Hep3Vector posECTShift(0.,0.,2*(m_ectZ+ectZHalfSize));
   Trk::Volume posECTVol(negECTVol,*(new HepTransform3D(Trk::s_idRotation,posECTShift)));
   if (m_adjustStatic && m_static3d) positiveECT = processVolume( &posECTVol,2,"Muon::Detectors::PositiveECT" ); 
   else if (m_adjustStatic) positiveECT = processVolume( &posECTVol,-1,"Muon::Detectors::PositiveECT" ); 
   else positiveECT = processVolume( &posECTVol,m_innerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::PositiveECT" ); 
   //checkVolume(positiveECT);
   // glue
   const Trk::TrackingVolume* negativeMuonInnerEndcap =
     m_trackingVolumeHelper->glueTrackingVolumeArrays(*negativeECT, Trk::positiveFaceXY,
						      *negativeMuonSmallWheel, Trk::negativeFaceXY, 
						      "Muon::Container::NegInnerEndcap");  
   const Trk::TrackingVolume* positiveMuonInnerEndcap = 
     m_trackingVolumeHelper->glueTrackingVolumeArrays(*positiveMuonSmallWheel, Trk::positiveFaceXY,
						      *positiveECT, Trk::negativeFaceXY,
						      "Muon::Container::PosInnerEndcap");  
   
// inner shields   
   double innerEndcapZHalfSize = 0.5*(m_innerEndcapZ - m_diskShieldZ);
   negInnerShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_innerShieldRadius,
							innerEndcapZHalfSize);
   Hep3Vector negInnerEndcapPosition(0.,0.,-m_diskShieldZ-innerEndcapZHalfSize);
   Trk::Volume negisVol(new HepTransform3D(Trk::s_idRotation,negInnerEndcapPosition),negInnerShieldBounds);
   negInnerShield = processShield(&negisVol,1,"Muons::Detectors::NegativeInnerShield");

   posInnerShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_innerShieldRadius,
							innerEndcapZHalfSize);
   Hep3Vector posInnerEndcapPosition(0.,0.,m_diskShieldZ+innerEndcapZHalfSize);
   Trk::Volume posisVol(new HepTransform3D(Trk::s_idRotation,posInnerEndcapPosition),posInnerShieldBounds);
   posInnerShield = processShield(&posisVol,1,"Muons::Detectors::PositiveInnerShield");

// outer Endcap
// build as bigWheel+buffer+outerWheel
// outer wheel
   double outerWheelZHalfSize = 0.5*(m_outerEndcapZ - m_outerWheel);
   negativeOuterWheelBounds = new Trk::CylinderVolumeBounds(m_outerShieldRadius,
							    m_outerBarrelRadius,
							    outerWheelZHalfSize);
   Hep3Vector negOuterWheelPosition(0.,0.,-m_outerEndcapZ+outerWheelZHalfSize);
   HepTransform3D* negOuterWheelTransf = new HepTransform3D(Trk::s_idRotation,negOuterWheelPosition);
   Trk::Volume negOWVol(negOuterWheelTransf,negativeOuterWheelBounds);
   if (m_adjustStatic && m_static3d) negativeMuonOuterWheel = processVolume( &negOWVol,3,"Muon::Detectors::NegativeOuterWheel" ); 
   else if (m_adjustStatic) negativeMuonOuterWheel = processVolume( &negOWVol,-1,"Muon::Detectors::NegativeOuterWheel" ); 
   else negativeMuonOuterWheel = processVolume( &negOWVol,m_outerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::NegativeOuterWheel" ); 
   //
   Hep3Vector posOuterWheelShift(0.,0.,2*(m_outerEndcapZ-outerWheelZHalfSize));
   Trk::Volume posOWVol(negOWVol,*(new HepTransform3D(Trk::s_idRotation,posOuterWheelShift)));
   if (m_adjustStatic && m_static3d) positiveMuonOuterWheel = processVolume( &posOWVol,3,"Muon::Detectors::PositiveOuterWheel" ); 
   else if (m_adjustStatic) positiveMuonOuterWheel = processVolume( &posOWVol,-1,"Muon::Detectors::PositiveOuterWheel" ); 
   else positiveMuonOuterWheel = processVolume( &posOWVol,m_outerEndcapEtaPartition,m_phiPartition,
						 "Muon::Detectors::PositiveOuterWheel" ); 
// outer buffer
   double outerBufferZHalfSize = 0.5*(m_outerWheel - m_bigWheel);
   negativeOuterBufferBounds = new Trk::CylinderVolumeBounds(m_outerShieldRadius,
							     m_outerBarrelRadius,
							     outerBufferZHalfSize);
   Hep3Vector negOuterBufferPosition(0.,0.,-m_bigWheel-outerBufferZHalfSize);
   HepTransform3D* negOuterBufferTransf = new HepTransform3D(Trk::s_idRotation,negOuterBufferPosition);
   Trk::Volume negBuffVol(negOuterBufferTransf,negativeOuterBufferBounds);
   if (m_adjustStatic && m_static3d) negativeMuonOuterBuffer = processVolume( &negBuffVol,3,"Muon::Detectors::NegativeOuterBuffer" ); 
   else if (m_adjustStatic) negativeMuonOuterBuffer = processVolume( &negBuffVol,-1,"Muon::Detectors::NegativeOuterBuffer" ); 
   else negativeMuonOuterBuffer = processVolume( &negBuffVol,m_outerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::NegativeOuterBuffer" ); 
   //
   Hep3Vector posOuterBufferShift(0.,0.,2*(m_bigWheel+outerBufferZHalfSize));
   Trk::Volume posBuffVol(negBuffVol,*(new HepTransform3D(Trk::s_idRotation,posOuterBufferShift)));
   if (m_adjustStatic && m_static3d) positiveMuonOuterBuffer = processVolume( &posBuffVol,3,"Muon::Detectors::PositiveOuterBuffer" ); 
   else if (m_adjustStatic) positiveMuonOuterBuffer = processVolume( &posBuffVol,-1,"Muon::Detectors::PositiveOuterBuffer" ); 
   else positiveMuonOuterBuffer = processVolume( &posBuffVol,m_outerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::PositiveOuterBuffer" ); 
// big wheel
   double bigWheelZHalfSize = 0.5*(m_bigWheel - m_innerEndcapZ);
   negativeBigWheelBounds = new Trk::CylinderVolumeBounds(m_outerShieldRadius,
							  m_outerBarrelRadius,
							  bigWheelZHalfSize);
   Hep3Vector negBigWheelPosition(0.,0.,-m_innerEndcapZ-bigWheelZHalfSize);
   HepTransform3D* negBigWheelTransf = new HepTransform3D(Trk::s_idRotation,negBigWheelPosition);
   Trk::Volume negBWVol(negBigWheelTransf,negativeBigWheelBounds);
   if (m_adjustStatic && m_static3d) negativeMuonBigWheel = processVolume( &negBWVol,3,"Muon::Detectors::NegativeBigWheel" ); 
   else if (m_adjustStatic) negativeMuonBigWheel = processVolume( &negBWVol,-1,"Muon::Detectors::NegativeBigWheel" ); 
   else negativeMuonBigWheel = processVolume( &negBWVol,m_outerEndcapEtaPartition,m_phiPartition,
						"Muon::Detectors::NegativeBigWheel" ); 
   //
   Hep3Vector posBigWheelShift(0.,0.,2*(m_innerEndcapZ+bigWheelZHalfSize));
   Trk::Volume posBWVol(negBWVol,*(new HepTransform3D(Trk::s_idRotation,posBigWheelShift)));
   if (m_adjustStatic && m_static3d) positiveMuonBigWheel = processVolume( &posBWVol,3,"Muon::Detectors::PositiveBigWheel" ); 
   else if (m_adjustStatic) positiveMuonBigWheel = processVolume( &posBWVol,-1,"Muon::Detectors::PositiveBigWheel" ); 
   else positiveMuonBigWheel = processVolume( &posBWVol,m_outerEndcapEtaPartition,m_phiPartition,
						 "Muon::Detectors::PositiveBigWheel" ); 
// glue
   const Trk::TrackingVolume* negNavOEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*negativeMuonOuterWheel, Trk::positiveFaceXY,
											       *negativeMuonOuterBuffer, Trk::negativeFaceXY, 
											       "Muon::Container::NegOEndcap");  
   const Trk::TrackingVolume* posNavOEndcap = m_trackingVolumeHelper->glueTrackingVolumeArrays(*positiveMuonOuterBuffer, Trk::positiveFaceXY,
											       *positiveMuonOuterWheel, Trk::negativeFaceXY,
											       "Muon::Container::PosOEndcap");  
   const Trk::TrackingVolume* negativeMuonOuterEndcap =
     m_trackingVolumeHelper->glueTrackingVolumeArrays(*negNavOEndcap, Trk::positiveFaceXY,
						      *negativeMuonBigWheel, Trk::negativeFaceXY, 
						      "Muon::Container::NegOuterEndcap");  
   const Trk::TrackingVolume* positiveMuonOuterEndcap =
     m_trackingVolumeHelper->glueTrackingVolumeArrays(*positiveMuonBigWheel, Trk::positiveFaceXY, 
						      *posNavOEndcap, Trk::negativeFaceXY,
						      "Muon::Container::PosOuterEndcap");  

// outer shields   
   double outerEndcapZHalfSize = 0.5*(m_outerEndcapZ-m_innerEndcapZ);
   double outerEndcapPosition  = 0.5*(m_outerEndcapZ+m_innerEndcapZ);
   Hep3Vector negOuterShieldPosition(0.,0.,-outerEndcapPosition);
   negOuterShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_outerShieldRadius,
							outerEndcapZHalfSize);
   Trk::Volume negosVol(new HepTransform3D(Trk::s_idRotation,negOuterShieldPosition),negOuterShieldBounds);
   negOuterShield = processShield(&negosVol,0,"Muons::Detectors::NegativeOuterShield");

   posOuterShieldBounds = new Trk::CylinderVolumeBounds(m_beamPipeRadius,
							m_outerShieldRadius,
							outerEndcapZHalfSize);
   Hep3Vector posOuterShieldPosition(0.,0.,outerEndcapPosition);
   Trk::Volume pososVol(new HepTransform3D(Trk::s_idRotation,posOuterShieldPosition),posOuterShieldBounds);
   posOuterShield = processShield(&pososVol,0,"Muons::Detectors::PositiveOuterShield");

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

   negBeamPipe->registerColorCode(0); 
   posBeamPipe->registerColorCode(0); 

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
											  m_exitVolume);
// blend material
   if (m_blendInertMaterial) blendMaterial();

// tracking geometry 
   Trk::TrackingGeometry* trackingGeometry = new Trk::TrackingGeometry(detector,Trk::globalSearch);

// clean-up
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
   
   for (size_t i = 0; i < m_spans.size(); i++) delete m_spans[i];
   
   for (std::map<const Trk::DetachedTrackingVolume*,std::vector<const Trk::TrackingVolume*>* >::iterator it = m_blendMap.begin();
	it != m_blendMap.end();
	++it)
     {
       delete it->second;
     }   

   m_chronoStatSvc->chronoStop("MS::build-up");

   log << MSG::INFO  << name() <<" returning tracking geometry " << endreq;    
   log << MSG::INFO  << name() <<" with "<< m_frameNum<<" subvolumes at navigation level" << endreq;    
   log << MSG::INFO  << name() <<"( mean number of enclosed detached volumes:"<< float(m_frameStat)/m_frameNum<<")" << endreq;    
   return trackingGeometry;  
}

// finalize
StatusCode Muon::MuonTrackingGeometryBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    if (m_stations) {
      for (size_t i = 0; i < m_stations->size(); i++) {
	if ((*m_stations)[i]) delete (*m_stations)[i];
        else log << MSG::DEBUG << name() << " station pointer corrupted ! " << endreq; 
      }
      delete m_stations; m_stations = 0;
    } 
    if (m_inertObjs) {
      unsigned int inLim = (m_blendInertMaterial && m_removeBlended) ? m_inertPerm : m_inertObjs->size();
      for (size_t i = 0; i < inLim; i++) {
	if ((*m_inertObjs)[i])	delete (*m_inertObjs)[i];
        else log << MSG::DEBUG << name() << " inert object pointer corrupted ! " << endreq; 
      }
      delete m_inertObjs; m_inertObjs = 0;
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
    edges.push_back( Trk::GlobalPosition(0.,0.,bcyl->halflengthZ()));
    edges.push_back( Trk::GlobalPosition(0.,0.,-bcyl->halflengthZ()));
  }
  if (cyl) {
    edges.push_back( Trk::GlobalPosition(0.,0.,cyl->halflengthZ()));
    edges.push_back( Trk::GlobalPosition(0.,0.,-cyl->halflengthZ()));
  }
  if (spb) {
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif      
    const std::vector<std::pair<double, double> > vtcs = spb->xyVertices();
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif      
    for (unsigned int i=0;i<vtcs.size();i++) {
      edges.push_back( Trk::GlobalPosition(vtcs[i].first,vtcs[i].second, spb->halflengthZ()) );
      edges.push_back( Trk::GlobalPosition(vtcs[i].first,vtcs[i].second, -spb->halflengthZ()) );
    }
  }
  if (prism) {
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif  
    const std::vector<std::pair<double, double> > vtcs = prism->xyVertices();
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif      
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
    log << MSG::DEBUG << "edges:"<< ie<<","<<gp<<","<< phi<< endreq;
    double rad = gp.perp();
    if (cyl || bcyl) {
      double radius = 0.; double hz = 0.;
      Trk::GlobalDirection dir = (transform*Trk::GlobalDirection(0.,0.,1.));
      double thAx = dir.theta();
      if (cyl) { radius = cyl->outerRadius(); hz = cyl->halflengthZ();}
      if (bcyl) { radius = bcyl->outerRadius(); hz = bcyl->halflengthZ();}
      if ( gp[2]-radius*sin(thAx) <minZ ) minZ = gp[2]-radius*sin(thAx);
      if ( gp[2]+radius*sin(thAx) >maxZ ) maxZ = gp[2]+radius*sin(thAx);
      if ( rad-radius*fabs(cos(thAx)) < minR )  minR = rad>radius ? rad-radius*fabs(cos(thAx)): 0;
      if ( rad+radius*fabs(cos(thAx)) > maxR )  maxR = rad+radius*fabs(cos(thAx));
      // distance of cylinder axis and global axis
      if (dir.perp()>0.) {
        // distance to minimal approach
        double dMA = fabs(dir[0]*gp[0]+dir[1]*gp[1])/dir.perp()/dir.perp();
        double dMD = sqrt (gp.perp()*gp.perp()-dMA*dMA);
        if (dMA<2*hz && dMD-radius < minR ) minR = fmax(0.,dMD-radius); 
      }
      double dph = rad>0.? atan(radius/rad) : M_PI;
      if ( phi-dph <M_PI && phi-dph < minP0 ) minP0 = phi-dph;
      if ( phi+dph <M_PI && phi+dph > maxP0 ) maxP0 = phi+dph;
      if ( phi-dph >M_PI && phi-dph < minP1 ) minP1 = phi-dph;
      if ( phi+dph >M_PI && phi+dph > maxP1 ) maxP1 = phi+dph;      
    } else {
      if ( gp[2]<minZ ) minZ = gp[2];
      if ( gp[2]>maxZ ) maxZ = gp[2];
      if ( phi<M_PI && phi < minP0 ) minP0 = phi; 
      if ( phi<M_PI && phi > maxP0 ) maxP0 = phi; 
      if ( phi>M_PI && phi < minP1 ) minP1 = phi; 
      if ( phi>M_PI && phi > maxP1 ) maxP1 = phi; 
      //if ( phi < minPhi ) minPhi = phi; 
      //if ( phi > maxPhi ) maxPhi = phi; 
      if ( rad < minR ) minR = rad; 
      if ( rad > maxR ) maxR = rad;
    } 
  }
  if (maxPhi<minPhi) {
    if (maxP0>=minP0 && maxP1<minP1) { minPhi = minP0; maxPhi = maxP0; }
    else if ( maxP1>=minP1 && maxP0<minP0) { minPhi = minP1; maxPhi = maxP1; }
    else if ( maxP1 - minP0 < (maxP0 - minP1+2*M_PI) ) { minPhi = minP0; maxPhi = maxP1; }
    else { minPhi = minP1 ; maxPhi = maxP0; }  
  }
  if ( box || trd || dtrd || spb ) {
    span.push_back( minZ - zTol );  
    span.push_back( maxZ + zTol );  
    span.push_back( minPhi - phiTol );  
    span.push_back( maxPhi + phiTol );  
    span.push_back( fmax(m_beamPipeRadius+0.001, minR - zTol) );  
    span.push_back( maxR + zTol );  
  } else if (bcyl || cyl ) {
    span.push_back( minZ - cylZcorr -zTol );
    span.push_back( maxZ + cylZcorr +zTol );
    span.push_back( minPhi - phiTol );
    span.push_back( maxPhi + phiTol );
    span.push_back( fmax(m_beamPipeRadius+0.001, minR - zTol) );  
    span.push_back( maxR + zTol );  
  } else {
    log << MSG::ERROR  << name() <<" volume shape not recognized: "<< endreq;
    for (int i=0; i<6; i++) span.push_back(0.);
  }
  const Muon::Span* newSpan=new Muon::Span(span);
  return newSpan;
}

const std::vector<std::vector<std::pair<const Trk::DetachedTrackingVolume*,const Muon::Span*> >* >* Muon::MuonTrackingGeometryBuilder::findVolumesSpan(const std::vector<const Trk::DetachedTrackingVolume*>*& objs, double zTol, double phiTol) const
{
  MsgStream log(msgSvc(), name());

  if (!objs || !objs->size()) return 0;
  std::vector<std::vector<std::pair<const Trk::DetachedTrackingVolume*,const Span*> >* >* spans =
   new std::vector<std::vector<std::pair<const Trk::DetachedTrackingVolume*,const Span*> >* >(9);
  // split MS into 9 blocks to speed up the build-up of geometry
  for (unsigned int i=0;i<9;i++) (*spans)[i] = new std::vector<std::pair<const Trk::DetachedTrackingVolume*,const Span*> >; 
  for (unsigned int iobj=0; iobj<objs->size(); iobj++) {
    HepTransform3D  transform = (*objs)[iobj]->trackingVolume()->transform();
    const Muon::Span* span = findVolumeSpan(&((*objs)[iobj]->trackingVolume()->volumeBounds()), transform, zTol, phiTol);
    log << MSG::DEBUG << "span:"<<(*objs)[iobj]->name()<< ","<<(*span)[0]<<","<< (*span)[1]<<","<<(*span)[2]<<","
    << (*span)[3]<<","<< (*span)[4]<<","<< (*span)[5] << endreq;  
    // negative outer wheel
    if ( (*span)[0] < -m_bigWheel ) (*spans)[0]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // negative big wheen
    if ( (*span)[0] < -m_innerEndcapZ && (*span)[1]>-m_bigWheel ) (*spans)[1]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // neg.ect
    if ( (*span)[0] < -m_ectZ && (*span)[1]>-m_innerEndcapZ )     (*spans)[2]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // neg.small whell
    if ( (*span)[0] < -m_diskShieldZ && (*span)[1]>-m_ectZ )          (*spans)[3]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // barrel  
    if ( (*span)[0] <  m_diskShieldZ && (*span)[1]> -m_diskShieldZ  &&
	 ((*span)[5]> m_innerBarrelRadius || (*span)[0]<-m_barrelZ || (*span)[1]>m_barrelZ)  )
                                                                  (*spans)[4]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // pos.small whell
    if ( (*span)[0] < m_ectZ && (*span)[1]> m_barrelZ )            (*spans)[5]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // pos.ect
    if ( (*span)[0] < m_innerEndcapZ && (*span)[1]> m_ectZ )       (*spans)[6]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // positive big wheel
    if ( (*span)[0] < m_bigWheel && (*span)[1]> m_innerEndcapZ )   (*spans)[7]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));
    // positive outer wheel
    if ( (*span)[1] >  m_bigWheel )                                (*spans)[8]->push_back(std::pair<const Trk::DetachedTrackingVolume*,const Span*>((*objs)[iobj],span));

    m_spans.push_back(span);                 // keep track of things to delete
  }

  return spans;
}


const Trk::TrackingVolume* Muon::MuonTrackingGeometryBuilder::processVolume(const Trk::Volume* vol, int etaN , int phiN, std::string volumeName) const
{
  MsgStream log(msgSvc(), name());

  const Trk::TrackingVolume* tVol = 0;

  unsigned int colorCode = m_colorCode;

  std::vector<const Trk::DetachedTrackingVolume*> blendVols; 

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
      if (colorCode>0) colorCode = 26 - colorCode;
      for (int phi = 0; phi < phiN; phi++) {
	if (colorCode>0) colorCode = 26 - colorCode;
        // define subvolume
        double posZ = (vol->center())[2]+ etaSect * (2.*eta+1.-etaN) ;
        HepTransform3D  transf( HepRotateZ3D( phiSect*(2*phi+1))*HepTranslateZ3D(posZ));
        const Trk::Volume* subVol= new Trk::Volume(*protVol, transf);     
        // enclosed muon objects ?   
	std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) ; 
        blendVols.clear();
        std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( subVol, blendVols);
        const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( *subVol,
								   m_muonMaterial,
								   m_muonMagneticField,
								   detVols,
								   volName );
        // statistics
        m_frameNum++ ; if (detVols) m_frameStat += detVols->size();  
        // prepare blending
        if (m_blendInertMaterial && blendVols.size()) {
          for (unsigned int id=0;id<blendVols.size();id++) {
	    if (!m_blendMap[blendVols[id]]) 
	      m_blendMap[blendVols[id]] = new std::vector<const Trk::TrackingVolume*>;
	    m_blendMap[blendVols[id]]->push_back(sVol);
	  }
	}  
        //
        sVol->registerColorCode(colorCode); 
	// reference position 
	HepPoint3D gp(subBds->outerRadius(),0.,0.);
	subVolumes.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
                                                             Trk::GlobalPosition(transf*gp)));
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
    blendVols.clear();
    std::vector<const Trk::DetachedTrackingVolume*>* muonObjs = getDetachedObjects( vol, blendVols);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    muonObjs,
				    volumeName);
    // statistics
    m_frameNum++ ; if (muonObjs) m_frameStat += muonObjs->size();  
    // prepare blending
    if (m_blendInertMaterial && blendVols.size()) {
      for (unsigned int id=0;id<blendVols.size();id++) {
	if (!m_blendMap[blendVols[id]]) 
	  m_blendMap[blendVols[id]] = new std::vector<const Trk::TrackingVolume*>;
	m_blendMap[blendVols[id]]->push_back(tVol);
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

  unsigned int colorCode = m_colorCode;

  std::vector<const Trk::DetachedTrackingVolume* > blendVols;

  //getPartitionFromMaterial(vol);

  // retrieve cylinder
  const Trk::CylinderVolumeBounds* cyl=dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  if (!cyl) {
    log << MSG::ERROR << " process volume: volume cylinder boundaries not retrieved, return 0 " << endreq;
    return 0; 
  }
  // create vector of zSteps for this volume
  std::vector<double> zSteps;
  std::vector<int> zTypes;
  zSteps.clear(); zTypes.clear();
  double zPos = vol->center()[2];
  double hz = cyl->halflengthZ();
  double z1 = zPos-hz;  double z2 = zPos+hz ;
  zSteps.push_back(z1);
  for (unsigned int iz=0;iz<m_zPartitions.size();iz++) {
    if ( m_zPartitions[iz]==zSteps.front()) zTypes.push_back(m_zPartitionsType[iz]);
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

  // phi binning
  if (fabs(zPos)> m_barrelZ && cyl->outerRadius()<m_outerBarrelRadius) getPhiParts(0);
  else if (fabs(zPos)<= m_ectZ) getPhiParts(2);
  else if (fabs(zPos)<= m_innerEndcapZ) getPhiParts(3);
  else if (fabs(zPos)> m_outerWheel && cyl->outerRadius()> m_outerShieldRadius ) getPhiParts(1);
  else if (fabs(zPos)> m_innerEndcapZ && fabs(zPos)<m_bigWheel && cyl->outerRadius()> m_outerShieldRadius ) getPhiParts(1);
  else getPhiParts(0);

  // R/H binning ?
  unsigned int etaN = zSteps.size()-1;
  unsigned int phiN = m_adjustedPhi.size();   
  int phiTypeMax = 0;                // count different partitions

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

        if (m_adjustedPhiType[ip]>phiTypeMax) phiTypeMax = m_adjustedPhiType[ip];

        phBinUtil.push_back(new Trk::BinUtility1DH(phiRef,m_hPartitions[mode][zTypes[iz]][m_adjustedPhiType[ip]]));
      }
      hBinUtil.push_back(phBinUtil);
    }

    // create subvolumes & BinnedArray
    std::vector<Trk::TrackingVolumeOrderPosition>  subVolumesVect;
    std::vector<std::vector<std::vector<const Trk::TrackingVolume*> > > subVolumes;
    std::vector<std::vector<Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> > > > hBins;
    std::vector<const Trk::TrackingVolume*> sVolsInn;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsOut;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsNeg;             // for gluing
    std::vector<const Trk::TrackingVolume*> sVolsPos;             // for gluing
    for (unsigned int eta = 0; eta < zSteps.size()-1; eta++) {
      if (colorCode>0) colorCode = 6 -colorCode;
      double posZ = 0.5*(zSteps[eta] + zSteps[eta+1]) ;
      double   hZ = 0.5*fabs(zSteps[eta+1] - zSteps[eta]) ;
      std::vector<std::vector<const Trk::TrackingVolume*> > phiSubs;
      std::vector<Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> > >  phBins;
      std::vector<int> phiType(phiTypeMax+1,-1);        // indication of first phi/R partition built for a given type (for cloning)
      std::vector<std::vector<Trk::Volume*> > garbVol(phiTypeMax+1);
      unsigned int pCode = 1; 
      for (unsigned int phi = 0; phi < phiN; phi++) {
	pCode = (colorCode>0) ? 3-pCode : 0;
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
        int phiP = phiType[m_adjustedPhiType[phi]]; 

        unsigned int hCode = 1; 
        for (unsigned int h = 0; h < hSteps.size()-1; h++) {
	  hCode = colorCode>0 ? 1 - hCode : 0; 
          // similar volume may exist already
	  Trk::Volume* subVol=0;
	  HepTransform3D* transf = new HepTransform3D(HepRotateZ3D(posPhi)*HepTranslateZ3D(posZ));
          //
          int volType = 0;     // cylinder 
          if ( hSteps[h].first == 1 && hSteps[h+1].first == 0 ) volType = 1;  
          if ( hSteps[h].first == 0 && hSteps[h+1].first == 1 ) volType = 2;  
          if ( hSteps[h].first == 1 && hSteps[h+1].first == 1 ) volType = 3;  
	  // define subvolume
          if (phiP>-1 ) {
            subVol = new Trk::Volume(*(phiSubs[phiP][h]),(*transf)*phiSubs[phiP][h]->transform().inverse());
	  } else if ( phiSect<0.5*M_PI) {
	    Trk::BevelledCylinderVolumeBounds* subBds = new Trk::BevelledCylinderVolumeBounds(hSteps[h].second,
							   hSteps[h+1].second,
							   phiSect, 
							   hZ,
							   volType);
	    subVol = new Trk::Volume(transf, subBds);
	  } else {
	    Trk::CylinderVolumeBounds* subBds = new Trk::CylinderVolumeBounds(hSteps[h].second,
						   hSteps[h+1].second,
						   phiSect, 
						   hZ);
	    subVol = new Trk::Volume(transf, subBds);
	  }
       
	  // enclosed muon objects ? also adjusts material properties in case of material blend  
	  std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) +MuonGM::buildString(h,2) ; 
          blendVols.clear(); 
	  std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( subVol, blendVols);

	  const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( *subVol,
								     m_muonMaterial,
								     m_muonMagneticField,
								     detVols,
								     volName );
                                                                              
	  // statistics
	  m_frameNum++ ; if (detVols) m_frameStat += detVols->size();  
	  // prepare blending
	  if (m_blendInertMaterial && blendVols.size()) {
	    for (unsigned int id=0;id<blendVols.size();id++) {
	      if (!m_blendMap[blendVols[id]]) 
		m_blendMap[blendVols[id]] = new std::vector<const Trk::TrackingVolume*>;
	      m_blendMap[blendVols[id]]->push_back(sVol);
	    }
	  }  
          //
          sVol->registerColorCode(colorCode+pCode+hCode);
	  // reference position 
	  HepPoint3D gp(0.5*(hSteps[h].second+hSteps[h+1].second),0.,0.);
	  subVolumesVect.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, false),
	                                                       Trk::GlobalPosition((*transf)*gp)));
	  hSubsTr.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
	                                                       Trk::GlobalPosition((*transf)*gp)));
	  hSubs.push_back(sVol);

          // cleanup 
          if (phiP>-1) {delete transf; delete subVol;}
          else garbVol[m_adjustedPhiType[phi]].push_back(subVol);       // don't delete before cloned

	  //glue subVolume
	  if (h==0)                sVolsInn.push_back(sVol); 
	  if (h==hSteps.size()-2)  sVolsOut.push_back(sVol); 
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
        // save link to current partition for cloning 
        if (phiP<0) phiType[m_adjustedPhiType[phi]] = phi;    

        // finish phi gluing
        if (phiN>1 && phi>0) {
	  for (unsigned int j=0; j<phiSubs[phi-1].size(); j++) {
	    m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*(phiSubs[phi-1][j]), Trk::tubeSectorPositivePhi,phBins[phi]);
	  }
	}
        if (phiN>1 && phi==phiN-1) {
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
      // get rid of the garbage
      for (unsigned int j=0;j<garbVol.size();j++)
	for (unsigned int jj=0;jj<garbVol[j].size();jj++) delete garbVol[j][jj];
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
        blendVols.clear();
        std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( &subVol, blendVols);
        const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( subVol,
								   m_muonMaterial,
								   m_muonMagneticField,
								   detVols,
								   volName );
        // statistics
        m_frameNum++ ; if (detVols) m_frameStat += detVols->size();  
        // prepare blending
        if (m_blendInertMaterial && blendVols.size()) {
          for (unsigned int id=0;id<blendVols.size();id++) {
	    if (!m_blendMap[blendVols[id]]) 
	      m_blendMap[blendVols[id]] = new std::vector<const Trk::TrackingVolume*>;
	    m_blendMap[blendVols[id]]->push_back(sVol);
	  }
	}  
        //delete subVol;
        sVol->registerColorCode(colorCode); 
        // reference position 
	HepPoint3D gp(subBds->outerRadius(),0.,0.);
	//subVolumes.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
        //                                                     Trk::GlobalPosition((*transf)*gp)));
	subVolumes[phi*etaN+eta] = Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, false),
                                                             Trk::GlobalPosition((*transf)*gp));
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
    blendVols.clear();
    std::vector<const Trk::DetachedTrackingVolume*>* muonObjs = getDetachedObjects( vol, blendVols);

    tVol = new Trk::TrackingVolume( *vol,
                                    m_muonMaterial,
				    m_muonMagneticField,
				    muonObjs,
				    volumeName);
    // statistics
    m_frameNum++ ; if (muonObjs) m_frameStat += muonObjs->size();  
    // prepare blending
    if (m_blendInertMaterial && blendVols.size()) {
      for (unsigned int id=0;id<blendVols.size();id++) {
	if (!m_blendMap[blendVols[id]]) 
	  m_blendMap[blendVols[id]] = new std::vector<const Trk::TrackingVolume*>;
	m_blendMap[blendVols[id]]->push_back(tVol);
      }
    }  
  }

  return tVol;
} 

const Trk::TrackingVolume* Muon::MuonTrackingGeometryBuilder::processShield(const Trk::Volume* vol, int type,std::string volumeName) const
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << name() << "processing shield volume in mode:"<< type << endreq;

  const Trk::TrackingVolume* tVol = 0;

  unsigned int colorCode = m_colorCode;

  std::vector<const Trk::DetachedTrackingVolume*> blendVols;

  //getPartitionFromMaterial(vol);

  // retrieve cylinder
  const Trk::CylinderVolumeBounds* cyl=dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  if (!cyl) {
    log << MSG::ERROR << " process volume: volume cylinder boundaries not retrieved, return 0 " << endreq;
    return 0; 
  }
  // create vector of zSteps for this volume
  std::vector<double> zSteps;
  zSteps.clear();
  double zPos = vol->center()[2];
  double hz = cyl->halflengthZ();
  double z1 = zPos-hz;  double z2 = zPos+hz ;
  zSteps.push_back(z1);
  for (unsigned int iz=0;iz<m_shieldZPart.size();iz++) {
    if ( m_shieldZPart[iz]> z1 && m_shieldZPart[iz] < z2 ) {
      zSteps.push_back(m_shieldZPart[iz]);
      z1 = m_shieldZPart[iz];
    }
  }
  zSteps.push_back(z2);

  // phi binning trivial
  m_adjustedPhi.clear();
  m_adjustedPhi.push_back(0.);

  unsigned int etaN = zSteps.size()-1;

  // create z,h bin utilities
  Trk::BinUtility1DZZ* zBinUtil = new Trk::BinUtility1DZZ(zSteps);
  Trk::BinUtility1DF* pBinUtil = new Trk::BinUtility1DF(m_adjustedPhi);
  std::vector<std::vector<Trk::BinUtility1D*> >  hBinUtil;
  double phiRef = 0.;
  for (unsigned iz=0;iz < zSteps.size()-1; iz++) {
    std::vector<Trk::BinUtility1D*> phBinUtil;
    phBinUtil.push_back(new Trk::BinUtility1DH(phiRef,m_shieldHPart[type]) );
    hBinUtil.push_back(phBinUtil);
  }
  // subvolume boundaries
  Trk::CylinderVolumeBounds* subBds=0;

  // create subvolumes & BinnedArray
  std::vector<Trk::TrackingVolumeOrderPosition>  subVolumesVect;
  std::vector<std::vector<std::vector<const Trk::TrackingVolume*> > > subVolumes;
  std::vector<std::vector<Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> > > > hBins;
  std::vector<const Trk::TrackingVolume*> sVolsInn;             // for gluing
  std::vector<const Trk::TrackingVolume*> sVolsOut;             // for gluing
  std::vector<const Trk::TrackingVolume*> sVolsNeg;             // for gluing
  std::vector<const Trk::TrackingVolume*> sVolsPos;             // for gluing
  for (unsigned int eta = 0; eta < zSteps.size()-1; eta++) {
    if (colorCode>0) colorCode = 26 -colorCode;
    double posZ = 0.5*(zSteps[eta] + zSteps[eta+1]) ;
    double   hZ = 0.5*fabs(zSteps[eta+1] - zSteps[eta]) ;
    std::vector<std::vector<const Trk::TrackingVolume*> > phiSubs;
    std::vector<Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> > >  phBins;
    int phi = 0;
    double posPhi = 0.;
    double phiSect = M_PI;
    std::vector<std::pair<int,double> > hSteps =  m_shieldHPart[type];
    std::vector<const Trk::TrackingVolume*> hSubs;
    std::vector<Trk::TrackingVolumeOrderPosition> hSubsTr;
    unsigned int hCode = 1; 
    for (unsigned int h = 0; h < hSteps.size()-1; h++) {
      hCode = (colorCode>0) ? 1 - hCode : 0; 
      // define subvolume
      subBds = new Trk::CylinderVolumeBounds(hSteps[h].second,
					     hSteps[h+1].second,
					     phiSect, 
					     hZ);
      HepTransform3D* transf = new HepTransform3D(HepRotateZ3D(posPhi)*HepTranslateZ3D(posZ));
      Trk::Volume subVol(transf, subBds);
      
      // enclosed muon objects ? also adjusts material properties in case of material blend  
      std::string volName = volumeName +MuonGM::buildString(eta,2) +MuonGM::buildString(phi,2) +MuonGM::buildString(h,2) ; 
      blendVols.clear();
      std::vector<const Trk::DetachedTrackingVolume*>* detVols= getDetachedObjects( &subVol, blendVols);
      
      const Trk::TrackingVolume* sVol = new Trk::TrackingVolume( subVol,
								 m_muonMaterial,
								 m_muonMagneticField,
								 detVols,
								 volName );
      
      // statistics
      m_frameNum++ ; if (detVols) m_frameStat += detVols->size();  
      // prepare blending
      if (m_blendInertMaterial && blendVols.size()) {
	for (unsigned int id=0;id<blendVols.size();id++) {
	  if (!m_blendMap[blendVols[id]]) 
	    m_blendMap[blendVols[id]] = new std::vector<const Trk::TrackingVolume*>;
	  m_blendMap[blendVols[id]]->push_back(sVol);
	}
      }  
      //
      sVol->registerColorCode(colorCode+hCode);
      // reference position 
      HepPoint3D gp(subBds->mediumRadius(),0.,0.);
      subVolumesVect.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, false),
								Trk::GlobalPosition((*transf)*gp)));
      hSubsTr.push_back(Trk::TrackingVolumeOrderPosition(Trk::SharedObject<const Trk::TrackingVolume>(sVol, true),
							 Trk::GlobalPosition((*transf)*gp)));
      hSubs.push_back(sVol);
      
      //glue subVolume
      if (h==0)                sVolsInn.push_back(sVol); 
      if (h==hSteps.size()-2)  sVolsOut.push_back(sVol); 
      if (eta==0)      sVolsNeg.push_back(sVol); 
      if (eta==etaN-1) sVolsPos.push_back(sVol); 
      // in R/H
      if (h>0) { // glue 'manually'
	m_trackingVolumeHelper->setInsideTrackingVolume(*sVol, Trk::tubeSectorInnerCover,hSubs[h-1]); 
	m_trackingVolumeHelper->setOutsideTrackingVolume(*(hSubs[h-1]), Trk::tubeSectorOuterCover,sVol);
      }
      // in eta
      if ( etaN>1 && eta>0) m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*sVol, Trk::negativeFaceXY, hBins[eta-1][phi]);     
    }	
    phiSubs.push_back(hSubs); 
    Trk::BinnedArray1D<Trk::TrackingVolume>* volBinArray = new Trk::BinnedArray1D<Trk::TrackingVolume>(hSubsTr,hBinUtil[eta][phi]->clone());
    phBins.push_back(Trk::SharedObject<Trk::BinnedArray<Trk::TrackingVolume> >(volBinArray));

    // finish eta gluing
    if ( etaN>1 && eta>0) {
      for (unsigned int j=0; j<subVolumes[eta-1][phi].size(); j++) {
	m_trackingVolumeHelper->setOutsideTrackingVolumeArray(*(subVolumes[eta-1][phi][j]), Trk::positiveFaceXY,phBins[phi]);
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

std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonTrackingGeometryBuilder::getDetachedObjects(const Trk::Volume* vol, std::vector<const Trk::DetachedTrackingVolume*>& blendVols, int mode ) const
{
  // mode : 0 all, 1 active only, 2 inert only

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

  // define detector region
  int gMode = (zMax <= -m_bigWheel) ? 0 : 1;
  if ( zMin >= m_bigWheel ) gMode = 8;
  else if ( zMin >= m_innerEndcapZ ) gMode = 7;
  else if ( zMin >= m_ectZ )         gMode = 6;
  else if ( zMin >= m_diskShieldZ )      gMode = 5;
  else if ( zMin >=-m_diskShieldZ )      gMode = 4;
  else if ( zMin >=-m_ectZ )         gMode = 3;
  else if ( zMin >=-m_innerEndcapZ ) gMode = 2;
   
  std::list<const Trk::DetachedTrackingVolume*> detached;
  // active, use corrected rMax
  if (mode <2 && m_stationSpan) {
    for (unsigned int i=0; i<(*m_stationSpan)[gMode]->size() ; i++) {
      const Muon::Span* s = (*((*m_stationSpan)[gMode]))[i].second;          // span
      const Trk::DetachedTrackingVolume* station = (*((*m_stationSpan)[gMode]))[i].first;   // station
      bool rLimit = !m_static3d || ( (*s)[4] <= rMaxc && (*s)[5] >= rMin );
   
      if ( rLimit && (*s)[0] < zMax && (*s)[1] > zMin ) {
	if (phiLim) {
	  if (pMin>=0 && pMax<=2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && (*s)[2] <= pMax && (*s)[3] >= pMin ) detached.push_back(station);
	    if ( (*s)[2]>(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin) ) detached.push_back(station);
	  } else if (pMin < 0) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin+2*M_PI) ) detached.push_back(station);
	    if ( (*s)[2]>(*s)[3]  ) detached.push_back(station);
	  } else if (pMax > 2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax-2*M_PI || (*s)[3] >= pMin) ) detached.push_back(station);
	    if ( (*s)[2]>(*s)[3]  ) detached.push_back(station);
	  }
	} else {
	  detached.push_back(station);
	}
      } 
    }
  }
  // passive 
  if (mode !=1 && m_inertSpan) {
    for (unsigned int i=0; i<(*m_inertSpan)[gMode]->size() ; i++) {
      const Muon::Span* s = (*((*m_inertSpan)[gMode]))[i].second;
      const Trk::DetachedTrackingVolume* inert = (*((*m_inertSpan)[gMode]))[i].first;
      //bool rail = ( (*m_inertObjs)[i]->name() == "Rail" ) ? true : false; 
      bool rLimit = (!m_static3d || ( (*s)[4] <= rMaxc && (*s)[5] >= rMin ) ); 
      if ( rLimit && (*s)[0] < zMax && (*s)[1] > zMin ) {
        bool accepted = false;
	if (phiLim) {
	  if (pMin>=0 && pMax<=2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && (*s)[2] <= pMax && (*s)[3] >= pMin )  accepted = true;
	    if ( (*s)[2]>(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin) ) accepted = true;
	  } else if (pMin < 0) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax || (*s)[3] >= pMin+2*M_PI) ) accepted = true;
	    if ( (*s)[2]>(*s)[3]  ) accepted = true;
	  } else if (pMax > 2*M_PI) {
	    if ( (*s)[2]<=(*s)[3] && ((*s)[2] <= pMax-2*M_PI || (*s)[3] >= pMin) ) accepted = true;
	    if ( (*s)[2]>(*s)[3]  ) accepted = true;
	  }
	} else  accepted = true;
	if (accepted) {
          bool perm = inert->name().substr(inert->name().size()-4,4)=="PERM";
          if ( !m_blendInertMaterial || !m_removeBlended || perm ) detached.push_back(inert);
          if ( m_blendInertMaterial && !perm ) blendVols.push_back(inert);
	}
      } 
    }
  }
  if (!detached.empty()) detTVs = new std::vector<const Trk::DetachedTrackingVolume*>(detached.begin(), detached.end()); 
  return detTVs;
}

bool Muon::MuonTrackingGeometryBuilder::enclosed(const Trk::Volume* vol, const Muon::Span* s) const
{
  MsgStream log(msgSvc(), name());

  bool encl = false;
  double tol = 1.;
  double ptol = 0.11;     // 0.08 for BT, 0.11 feet
  
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
  //std::auto_ptr<const Muon::Span> s (findVolumeSpan(&(cs->volumeBounds()), cs->transform(), 0.,0.) );
  log << MSG::DEBUG << "enclosing volume:z:"<< zMin<<","<<zMax<<":r:"<< rMin<<","<<rMax<<":phi:"<<pMin<<","<<pMax<< endreq;
  //
  bool rLimit = (!m_static3d || ( (*s)[4] < rMax-tol && (*s)[5] > rMin+tol ) ); 
  if ( rLimit && (*s)[0] < zMax-tol && (*s)[1] > zMin+tol ) {
    if (phiLim) {
      if (pMin>=0 && pMax<=2*M_PI) {
	if ( (*s)[2]<=(*s)[3] && (*s)[2] < pMax-ptol && (*s)[3] > pMin+ptol )  return true;
	if ( (*s)[2]>(*s)[3] && ((*s)[2] < pMax-ptol || (*s)[3] > pMin+ptol) ) return true;
      } else if (pMin < 0) {
	if ( (*s)[2]<=(*s)[3] && ((*s)[2] < pMax-ptol || (*s)[3] > pMin+ptol+2*M_PI) ) return true;
	if ( (*s)[2]>(*s)[3]  ) return true;
      } else if (pMax > 2*M_PI) {
	if ( (*s)[2]<=(*s)[3] && ((*s)[2] < pMax-ptol-2*M_PI || (*s)[3] > pMin+ptol) ) return true;
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
  // activeAdjustLevel:  1: separate MDT stations
  //                        +(inertLevel=0) barrel Z partition
  //                     2: split TGC
  //                        +(inertLevel=0) barrel R partition
  //                     3: split TGC supports
  // inertAdjustLevel:   1: BT,ECT     
 
  // hardcode for the moment
  m_zPartitions.clear();
  m_zPartitionsType.clear();
  m_zPartitions.reserve(120);
  m_zPartitionsType.reserve(120);

  // outer endcap
  m_zPartitions.push_back(-m_outerEndcapZ);  m_zPartitionsType.push_back(1);   // EO
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-21630.);  m_zPartitionsType.push_back(1); }   // EOL
  m_zPartitions.push_back(-m_outerWheel);    m_zPartitionsType.push_back(0);   // Octogon
  m_zPartitions.push_back(-17990.);          m_zPartitionsType.push_back(0);   // buffer
  m_zPartitions.push_back(-m_bigWheel);      m_zPartitionsType.push_back(1);   // TGC3 
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-15225.);  m_zPartitionsType.push_back(1);}    
  if (m_activeAdjustLevel>1) { m_zPartitions.push_back(-15172.);  m_zPartitionsType.push_back(1);}   // end supp 
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-15128.);  m_zPartitionsType.push_back(1);}   // supp 
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-15070.);  m_zPartitionsType.push_back(1);}    
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-14940.);  m_zPartitionsType.push_back(1); }  //
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-14805.);  m_zPartitionsType.push_back(1);}    
  if (m_activeAdjustLevel>1) { m_zPartitions.push_back(-14733.);  m_zPartitionsType.push_back(1);}   // end supp.
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-14708.);  m_zPartitionsType.push_back(1);}   // supp.   
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-14650.);  m_zPartitionsType.push_back(1);}   //    
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-14560.);  m_zPartitionsType.push_back(1); }   // EML 
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-14080.);  m_zPartitionsType.push_back(1); }   // EMS
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-13620.);  m_zPartitionsType.push_back(1); }   // TGC
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-13525.);  m_zPartitionsType.push_back(1); }   // TGC
  if (m_activeAdjustLevel>1) { m_zPartitions.push_back(-13448.5); m_zPartitionsType.push_back(1); }   // end supp.
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-13421.5); m_zPartitionsType.push_back(1); }   // supp.
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-13346);   m_zPartitionsType.push_back(1); }   // TGC

  // inner endcap
  m_zPartitions.push_back(-m_innerEndcapZ);     m_zPartitionsType.push_back(0);   // 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-12790);  m_zPartitionsType.push_back(0); }   // ECT 
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-12100.); m_zPartitionsType.push_back(0); }   // 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-12000.); m_zPartitionsType.push_back(0); }   // 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-11210.); m_zPartitionsType.push_back(1); }   // BT 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-10480.); m_zPartitionsType.push_back(0); }   // 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-9700.);  m_zPartitionsType.push_back(0); }   // 
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-9300.);  m_zPartitionsType.push_back(0); }   // rib
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-8800.);  m_zPartitionsType.push_back(0); }   // ect
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-8610.);  m_zPartitionsType.push_back(1); }   // BT
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-8000.);  m_zPartitionsType.push_back(1); }   // BT
  m_zPartitions.push_back(-m_ectZ);  m_zPartitionsType.push_back(0);              // ECT/small whell
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-7450.); m_zPartitionsType.push_back(0); }   // EIS 
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-7364.); m_zPartitionsType.push_back(0); }   // EIS 
  if (m_activeAdjustLevel>0 || m_inertAdjustLevel>0) {m_zPartitions.push_back(-7170.); m_zPartitionsType.push_back(0);}   // cone assembly,TGC 
  if (m_activeAdjustLevel>0) { m_zPartitions.push_back(-7030.); m_zPartitionsType.push_back(0); }   // TGC
  if (m_activeAdjustLevel>2) { m_zPartitions.push_back(-6978.); m_zPartitionsType.push_back(0); }   // TGC

  // barrel
  m_zPartitions.push_back(-m_diskShieldZ);   m_zPartitionsType.push_back(0);   // disk 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-6829.);  m_zPartitionsType.push_back(0); }   // back disk 
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-6600.);  m_zPartitionsType.push_back(0); }   // 
  if (m_activeAdjustLevel>0){ m_zPartitions.push_back(-6100.);  m_zPartitionsType.push_back(0); }  
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-5503.);  m_zPartitionsType.push_back(1); }   // BT 
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-4772.);  m_zPartitionsType.push_back(0); }   //  
  if (m_activeAdjustLevel>0){ m_zPartitions.push_back(-4300.);  m_zPartitionsType.push_back(0); }  //  
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-3700.);  m_zPartitionsType.push_back(0); }   //  
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-3300.);  m_zPartitionsType.push_back(0); }   //  
  if (m_activeAdjustLevel>0){ m_zPartitions.push_back(-2600.);  m_zPartitionsType.push_back(0); }   //  
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-2078.);  m_zPartitionsType.push_back(1); }   // BT  
  if (m_inertAdjustLevel>0) { m_zPartitions.push_back(-1347.);  m_zPartitionsType.push_back(1); }   //  cryoring 
  if (m_activeAdjustLevel>0){ m_zPartitions.push_back(-800.);   m_zPartitionsType.push_back(1); }  //  cryoring 
  if (m_inertAdjustLevel>1) { m_zPartitions.push_back(-300.);   m_zPartitionsType.push_back(0); }   //   
  if (m_inertAdjustLevel+m_activeAdjustLevel<1) { m_zPartitions.push_back(-0.7*m_diskShieldZ);   m_zPartitionsType.push_back(0); }   //

  unsigned int zSiz = m_zPartitions.size();
  for (unsigned int i = 0; i<zSiz ; i++) {
    m_zPartitions.push_back(- m_zPartitions[zSiz-1-i]);  
    if (i<zSiz-1) m_zPartitionsType.push_back(m_zPartitionsType[zSiz-2-i]);  
  } 
  
  return;
}

void Muon::MuonTrackingGeometryBuilder::getPhiParts(int mode) const
{
  if (mode==0) {             // trivial
    m_adjustedPhi.clear();
    m_adjustedPhiType.clear();
    m_adjustedPhi.push_back(0.);
    m_adjustedPhiType.push_back(0);

  } else if (mode==1) {  
    int phiNum = 1;
    if ( m_activeAdjustLevel>0 ) phiNum = m_phiPartition;
    m_adjustedPhi.resize(phiNum);
    m_adjustedPhiType.resize(phiNum);
    m_adjustedPhi[0] = 0.;
    m_adjustedPhiType[0] = 0;
    int ic = 0;
    while (ic < phiNum-1 ) {
      ic++; 
      m_adjustedPhi[ic]= m_adjustedPhi[ic-1] + 2.*M_PI/phiNum ;
      m_adjustedPhiType[ic]= 0;
    }

  } else if (mode==2) {        // barrel(BT)
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

  } else if (mode==3) {      // ECT(+BT)
    // hardcode for the moment
    m_adjustedPhi.resize(32);
    m_adjustedPhiType.resize(32);
    
    double phiSect[3];
    phiSect[0] = 0.126;
    phiSect[2] = 0.105;
    phiSect[1] = 0.5*( M_PI/8. - phiSect[0]-phiSect[2] ); 
    
    m_adjustedPhi[0]= - phiSect[0];
    m_adjustedPhiType[0]= 0;
    m_adjustedPhi[1]= m_adjustedPhi[0]+2*phiSect[0];
    m_adjustedPhiType[1]= 1;
    int ic = 1; int is = 0;

    while (ic < 31 ) {
      ic++; is = 2 - is;
      m_adjustedPhi[ic]= m_adjustedPhi[ic-1] + 2*phiSect[1] ;
      m_adjustedPhiType[ic]= is;
      ic++;
      m_adjustedPhi[ic]= m_adjustedPhi[ic-1] + 2*phiSect[is] ;
      m_adjustedPhiType[ic]= 1;
    }
    //for (unsigned int ic=0;ic<m_adjustedPhi.size();ic++) std::cout<<"adjustedPhi:"<<ic<<","<<m_adjustedPhi[ic]<<","<<m_adjustedPhiType[ic]<< std::endl;
  }

  return;
}

void Muon::MuonTrackingGeometryBuilder::getHParts() const
{
  // hardcode for the moment
  m_hPartitions.clear();              // barrel, inner endcap, outer endcap

  // 0: barrel 2x2
  // non BT sector
  std::vector<std::pair<int,double> >  barrelZ0F0;
  barrelZ0F0.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  if (m_activeAdjustLevel>0) {
    barrelZ0F0.push_back( std::pair<int,double>(0,4450.) );                // for DiskShieldingBackDisk
    barrelZ0F0.push_back( std::pair<int,double>(0,6500.) );                // BI/BM
    barrelZ0F0.push_back( std::pair<int,double>(0,8900.) );                // BM/BO
  }
  barrelZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  barrelZ0F1;
  barrelZ0F1.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  if (m_inertAdjustLevel>0) {
    barrelZ0F1.push_back( std::pair<int,double>(1,4500.) );
    barrelZ0F1.push_back( std::pair<int,double>(1,5900.) );
  } else if (m_activeAdjustLevel>0) barrelZ0F1.push_back( std::pair<int,double>(0,4450.) );                
  if (m_activeAdjustLevel>0) barrelZ0F1.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0)  barrelZ0F1.push_back( std::pair<int,double>(1,8900.) );
  else if (m_activeAdjustLevel>0)  barrelZ0F1.push_back( std::pair<int,double>(0,8900.) );
  if (m_inertAdjustLevel>0)  barrelZ0F1.push_back( std::pair<int,double>(1,10100.) );
  barrelZ0F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // BT sector
  std::vector<std::pair<int,double> >  barrelZ1F0;
  barrelZ1F0.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  if (m_activeAdjustLevel+m_inertAdjustLevel>0) barrelZ1F0.push_back( std::pair<int,double>(0,4450.) );                
  if (m_inertAdjustLevel>0) {
    barrelZ1F0.push_back( std::pair<int,double>(1,5800.) );
    barrelZ1F0.push_back( std::pair<int,double>(1,6500.) );
  } else if (m_activeAdjustLevel>0) barrelZ1F0.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0) {
    barrelZ1F0.push_back( std::pair<int,double>(1,6750.) );
    barrelZ1F0.push_back( std::pair<int,double>(1,8400.) );
  }
  if (m_activeAdjustLevel>0) barrelZ1F0.push_back( std::pair<int,double>(0,8750.) );  // adapted for cryoring (from 8900)
  if (m_inertAdjustLevel>0) barrelZ1F0.push_back( std::pair<int,double>(1,9850.) );   // adapted for cryoring (from 9600)
  barrelZ1F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  barrelZ1F1;
  barrelZ1F1.push_back( std::pair<int,double>(0,m_innerBarrelRadius) );
  if (m_inertAdjustLevel>0) {
    barrelZ1F1.push_back( std::pair<int,double>(1,4500.) );
    barrelZ1F1.push_back( std::pair<int,double>(1,6000.) );
  } else if (m_activeAdjustLevel>0) barrelZ1F1.push_back( std::pair<int,double>(0,4450.) );
  if (m_activeAdjustLevel>0) barrelZ1F1.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0) barrelZ1F1.push_back( std::pair<int,double>(1,6800.) );
  if (m_inertAdjustLevel>0) {
    barrelZ1F1.push_back( std::pair<int,double>(1,8900.) );
    barrelZ1F1.push_back( std::pair<int,double>(1,10100.) );
  } else if (m_activeAdjustLevel>0) barrelZ1F1.push_back( std::pair<int,double>(0,8900.) );
  barrelZ1F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  barrelZF(2);
  barrelZF[0].push_back(barrelZ0F0);
  barrelZF[0].push_back(barrelZ0F1);
  barrelZF[1].push_back(barrelZ1F0);
  barrelZF[1].push_back(barrelZ1F1);

  // small wheel 1x2 ( no z BT sector) 
  // non BT sector
  std::vector<std::pair<int,double> >  swZ0F0;
  swZ0F0.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_activeAdjustLevel>1) {
    swZ0F0.push_back( std::pair<int,double>(0,2700.) );                
  }
  if (m_activeAdjustLevel+m_inertAdjustLevel>0) swZ0F0.push_back( std::pair<int,double>(0,4450.) );                
  if (m_activeAdjustLevel>0) {
    swZ0F0.push_back( std::pair<int,double>(0,6560.) );                // BI/BM
    swZ0F0.push_back( std::pair<int,double>(0,8900.) );                // BM/BO
  }
  swZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // phi BT sector
  std::vector<std::pair<int,double> >  swZ0F1;
  swZ0F1.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_activeAdjustLevel>1) swZ0F1.push_back( std::pair<int,double>(0,2700.) );               
  if (m_inertAdjustLevel+m_activeAdjustLevel>0) swZ0F1.push_back( std::pair<int,double>(0,4450.) );
  if (m_inertAdjustLevel>0) swZ0F1.push_back( std::pair<int,double>(1,5900.) );
  if (m_activeAdjustLevel>0) swZ0F1.push_back( std::pair<int,double>(0,6560.) );
  if (m_inertAdjustLevel>0) {
    swZ0F1.push_back( std::pair<int,double>(1,8900.) );
    swZ0F1.push_back( std::pair<int,double>(1,10100.) );
  } else if (m_activeAdjustLevel>0) swZ0F1.push_back( std::pair<int,double>(0,8900.) );
  swZ0F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  swZF(1);
  swZF[0].push_back(swZ0F0);
  swZF[0].push_back(swZ0F1);

  // inner endcap/ECT 2x3
  // ect coil, non-BT z
  std::vector<std::pair<int,double> >  innerZ0F0;
  innerZ0F0.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_inertAdjustLevel>0) innerZ0F0.push_back( std::pair<int,double>(0,1100.) );
  if (m_inertAdjustLevel>1) innerZ0F0.push_back( std::pair<int,double>(1,5150.) );
  if (m_inertAdjustLevel>0) innerZ0F0.push_back( std::pair<int,double>(1,5300.) );
  if (m_activeAdjustLevel>0) {
    innerZ0F0.push_back( std::pair<int,double>(0,6500.) );
    innerZ0F0.push_back( std::pair<int,double>(0,8900.) );
  }
  innerZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // coil gap, non-BT z
  std::vector<std::pair<int,double> >  innerZ0F1;
  innerZ0F1.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_inertAdjustLevel>0) innerZ0F1.push_back( std::pair<int,double>(0,1100.) );
  if (m_inertAdjustLevel>1) {
    innerZ0F1.push_back( std::pair<int,double>(1,1400.) );
    innerZ0F1.push_back( std::pair<int,double>(1,1685.) );
  }
  if (m_inertAdjustLevel>0) {
    innerZ0F1.push_back( std::pair<int,double>(1,4700.) );
    innerZ0F1.push_back( std::pair<int,double>(1,5900.) );
  }
  if (m_activeAdjustLevel>0) {
    innerZ0F1.push_back( std::pair<int,double>(0,6500.) );
    innerZ0F1.push_back( std::pair<int,double>(0,8900.) );
  }
  innerZ0F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // BT coil, no-BT z 
  std::vector<std::pair<int,double> >  innerZ0F2;
  innerZ0F2.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_inertAdjustLevel>0) innerZ0F2.push_back( std::pair<int,double>(0,1100.) );
  if (m_inertAdjustLevel>1) {
    innerZ0F2.push_back( std::pair<int,double>(1,1400.) );
    innerZ0F2.push_back( std::pair<int,double>(1,1685.) );
  }
  if (m_inertAdjustLevel>0) {
    innerZ0F2.push_back( std::pair<int,double>(1,4450.) );
    innerZ0F2.push_back( std::pair<int,double>(1,5900.) );
  }
  if (m_activeAdjustLevel>0) innerZ0F2.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0) {
    innerZ0F2.push_back( std::pair<int,double>(1,8900.) );
    innerZ0F2.push_back( std::pair<int,double>(1,10100.) );
  } else if (m_activeAdjustLevel>0) innerZ0F2.push_back( std::pair<int,double>(0,8900.) );
  innerZ0F2.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // ect coil, z BT sector
  std::vector<std::pair<int,double> >  innerZ1F0;
  innerZ1F0.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_inertAdjustLevel>0) innerZ1F0.push_back( std::pair<int,double>(0,1100.) );
  if (m_inertAdjustLevel>1) innerZ1F0.push_back( std::pair<int,double>(1,5150.) );
  if (m_inertAdjustLevel>0) innerZ1F0.push_back( std::pair<int,double>(1,5300.) );
  if (m_inertAdjustLevel>0) innerZ1F0.push_back( std::pair<int,double>(1,5800.) );
  if (m_inertAdjustLevel>0) innerZ1F0.push_back( std::pair<int,double>(1,6750.) );
  else if (m_activeAdjustLevel>0) innerZ1F0.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0) {
    innerZ1F0.push_back( std::pair<int,double>(1,8400.) );
    innerZ1F0.push_back( std::pair<int,double>(1,9600.) );
  } else if (m_activeAdjustLevel>0) innerZ1F0.push_back( std::pair<int,double>(0,8900.) );
  innerZ1F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // coil gap, BT z sector
  std::vector<std::pair<int,double> >  innerZ1F1;
  innerZ1F1.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_inertAdjustLevel>0)  innerZ1F1.push_back( std::pair<int,double>(0,1100.) );
  if (m_inertAdjustLevel>1) {
    innerZ1F1.push_back( std::pair<int,double>(1,1400.) );
    innerZ1F1.push_back( std::pair<int,double>(1,1685.) );
  }
  if (m_inertAdjustLevel>0) {
    innerZ1F1.push_back( std::pair<int,double>(1,4700.) );
    innerZ1F1.push_back( std::pair<int,double>(1,5800.) );
    innerZ1F1.push_back( std::pair<int,double>(1,6750.) );
  } else if (m_activeAdjustLevel>0) innerZ1F1.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0) {
    innerZ1F1.push_back( std::pair<int,double>(1,8400.) );
    innerZ1F1.push_back( std::pair<int,double>(1,9600.) );
  } else if (m_activeAdjustLevel>0) innerZ1F1.push_back( std::pair<int,double>(0,8900.) );
  innerZ1F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  // BT coil, BT z sector
  std::vector<std::pair<int,double> >  innerZ1F2;
  innerZ1F2.push_back( std::pair<int,double>(0,m_innerShieldRadius) );
  if (m_inertAdjustLevel>0) innerZ1F2.push_back( std::pair<int,double>(0,1100.) );
  if (m_inertAdjustLevel>1) {
    innerZ1F2.push_back( std::pair<int,double>(1,1400.) );
    innerZ1F2.push_back( std::pair<int,double>(1,1685.) );
  }
  innerZ1F2.push_back( std::pair<int,double>(0,4150.) );
  if (m_inertAdjustLevel>0) {
    innerZ1F2.push_back( std::pair<int,double>(1,4700.) );
    innerZ1F2.push_back( std::pair<int,double>(1,5900.) );
    innerZ1F2.push_back( std::pair<int,double>(1,6800.) );
  } else if (m_activeAdjustLevel>0) innerZ1F2.push_back( std::pair<int,double>(0,6500.) );
  if (m_inertAdjustLevel>0) {
    innerZ1F2.push_back( std::pair<int,double>(1,8900.) );
    innerZ1F2.push_back( std::pair<int,double>(1,10100.) );
  } else if (m_activeAdjustLevel>0) innerZ1F2.push_back( std::pair<int,double>(0,8900.) );
  innerZ1F2.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  innerZF(2);
  innerZF[0].push_back(innerZ0F0);
  innerZF[0].push_back(innerZ0F1);
  innerZF[0].push_back(innerZ0F2);
  innerZF[1].push_back(innerZ1F0);
  innerZF[1].push_back(innerZ1F1);
  innerZF[1].push_back(innerZ1F2);

  // outer 1x1
  std::vector<std::pair<int,double> >  outerZ0F0;
  outerZ0F0.push_back( std::pair<int,double>(0,m_outerShieldRadius) );
  outerZ0F0.push_back( std::pair<int,double>(0,2275.) );
  outerZ0F0.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::pair<int,double> >  outerZ0F1;
  outerZ0F1.push_back( std::pair<int,double>(0,m_outerShieldRadius) );
  if ( m_activeAdjustLevel>0 ) {
    outerZ0F1.push_back( std::pair<int,double>(0,3600.) );
    outerZ0F1.push_back( std::pair<int,double>(0,5300.) );
    outerZ0F1.push_back( std::pair<int,double>(0,7000.) );
    outerZ0F1.push_back( std::pair<int,double>(0,8500.) );
    outerZ0F1.push_back( std::pair<int,double>(0,10000.) );
    outerZ0F1.push_back( std::pair<int,double>(0,12000.) );
  }
  outerZ0F1.push_back( std::pair<int,double>(0,m_outerBarrelRadius) );

  std::vector<std::vector<std::vector<std::pair<int,double> > > >  outerZF(2);
  outerZF[0].push_back(outerZ0F0);
  outerZF[0].push_back(outerZ0F0);
  outerZF[1].push_back(outerZ0F1);
  outerZF[1].push_back(outerZ0F1);

  // collect everything
  m_hPartitions.push_back(barrelZF);
  m_hPartitions.push_back(swZF);
  m_hPartitions.push_back(innerZF);
  m_hPartitions.push_back(outerZF);
   
  return;
}

void Muon::MuonTrackingGeometryBuilder::getShieldParts() const
{
  m_shieldZPart.resize(18);

  m_shieldZPart[0]=-21900.;      // elm2
  m_shieldZPart[1]=-21500.;      // elm1
  m_shieldZPart[2]=-21000.;      // octogon
  m_shieldZPart[3]=-18000.;      // tube
  m_shieldZPart[4]=-12790.;      // ect
  m_shieldZPart[5]= -8000.;      // ect
  m_shieldZPart[6]= -7169.;      // cone
  m_shieldZPart[7]= -6829.;      // disk
  m_shieldZPart[8]= -6779.;       //

  for (unsigned int i = 9; i<18 ; i++) m_shieldZPart[i] = - m_shieldZPart[17-i];  

  m_shieldHPart.clear();

  std::vector<std::pair<int,double> >  outerShield;
  outerShield.push_back(std::pair<int,double>(0,m_beamPipeRadius));
  outerShield.push_back(std::pair<int,double>(0,600.));
  outerShield.push_back(std::pair<int,double>(0,1150.));
  outerShield.push_back(std::pair<int,double>(0,m_outerShieldRadius));
  m_shieldHPart.push_back(outerShield);

  std::vector<std::pair<int,double> >  innerShield;
  innerShield.push_back(std::pair<int,double>(0,m_beamPipeRadius));
  innerShield.push_back(std::pair<int,double>(0,530.));
  innerShield.push_back(std::pair<int,double>(0,m_innerShieldRadius));
  m_shieldHPart.push_back(innerShield);

  std::vector<std::pair<int,double> >  diskShield;
  diskShield.push_back(std::pair<int,double>(0,0.));
  diskShield.push_back(std::pair<int,double>(0,540.));
  diskShield.push_back(std::pair<int,double>(0,750.));
  diskShield.push_back(std::pair<int,double>(0,2700.));
  diskShield.push_back(std::pair<int,double>(0,m_innerBarrelRadius));
  m_shieldHPart.push_back(diskShield);

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
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif    
    std::vector<std::pair<double,double> > v=prism->xyVertices();
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif      
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
    if ( (*mIter).first->trackingVolume()->confinedDenseVolumes()) detMat = (*(*mIter).first->trackingVolume()->confinedDenseVolumes())[0];
    for (unsigned int ic=0; ic<cs->size(); ic++) {
      const Trk::Volume* nCs = new Trk::Volume(*((*cs)[ic].first),(*mIter).first->trackingVolume()->transform());
      double fraction = (*cs)[ic].second.second>0 ? 1. : (*cs)[ic].second.first;
      double csVol = fraction*calculateVolume(nCs);      
      const Muon::Span* s = findVolumeSpan(&(nCs->volumeBounds()), nCs->transform(), 0.,0.) ;
      if (s) log << MSG::DEBUG << "constituent:"<<ic<<":z:"<< (*s)[0]<<","<<(*s)[1]<<":r:"<< (*s)[4]<<","<<(*s)[5]
	    <<":phi:"<<(*s)[2]<<","<<(*s)[3]<< endreq;      
      double enVol = 0.;
      // loop over frame volumes, check if confined
      std::vector<const Trk::TrackingVolume*>::iterator fIter = (*mIter).second->begin(); 
      std::vector<bool> fEncl; 
      fEncl.clear();
      // blending factors can be saved, and not recalculated for each clone
      for ( ; fIter!=(*mIter).second->end(); fIter++) {
        fEncl.push_back(enclosed(*fIter,s));
        if ( fEncl.back() ) enVol += calculateVolume(*fIter);
      }
      delete nCs; delete s;
      // diluting factor
      double dil =  enVol>0. ?  csVol/enVol : 0.;
      //std::cout << "const:dil:"<< ic<<","<<dil<< std::endl;
      if (dil>0.) { 
	for ( fIter=(*mIter).second->begin(); fIter!=(*mIter).second->end(); fIter++) { 
	  if (fEncl[fIter-(*mIter).second->begin()]) { (*fIter)->addMaterial(*detMat,dil); if (m_colorCode==0) (*fIter)->registerColorCode(12) ; 
	  log << MSG::DEBUG << (*fIter)->volumeName()<<" acquires material from "<<  (*mIter).first->name()<< endreq;  }
	}
	log << MSG::DEBUG << "diluting factor:"<< dil<<" for "<< (*mIter).first->name()<<","<<ic<<endreq;
      } else {
	log << MSG::DEBUG << "diluting factor:"<< dil<<" for "<< (*mIter).first->name()<<","<<ic<<endreq;
      }
    }
    if ( m_removeBlended ) {  log << MSG::DEBUG << "deleting "<< (*mIter).first->name()<< endreq; delete (*mIter).first; }
  }
}
