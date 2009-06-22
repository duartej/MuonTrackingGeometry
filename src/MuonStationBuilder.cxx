///////////////////////////////////////////////////////////////////
// MuonStationBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonStationBuilder.h"
#include "MuonTrackingGeometry/MuonStationTypeBuilder.h"
//MuonSpectrometer include
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/MuonStation.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "MuonReadoutGeometry/CscReadoutElement.h"
#include "MuonReadoutGeometry/TgcReadoutElement.h"
#include "MuonIdHelpers/MdtIdHelper.h"
#include "MuonIdHelpers/RpcIdHelper.h"
// Trk
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeHelper.h"
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkDetDescrUtils/BinUtility1DX.h"
#include "TrkDetDescrUtils/BinUtility1DY.h"
#include "TrkDetDescrUtils/BinnedArray.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/DoubleTrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkVolumes/BoundarySurfaceFace.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkMagFieldInterfaces/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/TrackingGeometry.h"
#include "TrkGeometry/Layer.h"
#include<fstream>
#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoShapeShift.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoTubs.h"
#include "GeoModelKernel/GeoCons.h"
//mw
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoShapeUnion.h"
#include "GeoModelKernel/GeoShapeIntersection.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoTrap.h"
#include "GeoModelKernel/GeoPgon.h"
#include "GeoModelKernel/GeoPara.h"
#include "GeoModelKernel/GeoVolumeCursor.h"
// StoreGate
#include "StoreGate/StoreGateSvc.h"

// BField
#include "BFieldAth/MagFieldAthena.h"

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

// STD
#include <map>

// Gaudi
#include "GaudiKernel/MsgStream.h"

#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoShapeShift.h"
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoVolumeCursor.h"

// constructor
Muon::MuonStationBuilder::MuonStationBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgTool(t,n,p),
  m_muonMgrLocation("MuonMgr"),
  m_magFieldMode(Trk::RealisticField),
  m_magFieldTool("Trk::MagneticFieldTool/AtlasMagneticFieldTool"),
  m_muonStationTypeBuilder("Muon::MuonStationTypeBuilder/MuonStationTypeBuilder"),
  m_trackingVolumeHelper("Trk::TrackingVolumeHelper/TrackingVolumeHelper"),
  m_buildBarrel(true),
  m_buildEndcap(true),
  m_buildCsc(true),
  m_buildTgc(true),
  m_resolveActiveLayers(true)
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("StationTypeBuilder",               m_muonStationTypeBuilder);
  declareProperty("MagneticFieldMode",                m_magFieldMode);  
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
  declareProperty("BuildBarrelStations",              m_buildBarrel);
  declareProperty("BuildEndcapStations",              m_buildEndcap);
  declareProperty("BuildCSCStations",                 m_buildCsc);
  declareProperty("BuildTGCStations",                 m_buildTgc);
  declareProperty("ResolveActiveLayers",              m_resolveActiveLayers);
}

// destructor
Muon::MuonStationBuilder::~MuonStationBuilder()
{}

// Athena standard methods
// initialize
StatusCode Muon::MuonStationBuilder::initialize()
{
    
    MsgStream log(msgSvc(), name());

    StatusCode s = AlgTool::initialize();
    if (s.isFailure()) log << MSG::INFO << "failing to initialize?" << endreq; 

    // Get DetectorStore service
    //
    StoreGateSvc* m_detStore=0;
    StatusCode ds = service("DetectorStore",m_detStore);
    if (ds.isFailure()) {
        log << MSG::FATAL << "DetectorStore service not found !" << endreq;
    }
    // get Muon Spectrometer Description Manager
    ds = m_detStore->retrieve(m_muonMgr);
    if (ds.isFailure()) {
        log << MSG::ERROR << "Could not get MuonDetectorManager, no layers for muons will be built. " << endreq;
    }
    
    m_mdtIdHelper = m_muonMgr-> mdtIdHelper(); 
    m_rpcIdHelper = m_muonMgr-> rpcIdHelper();
    m_cscIdHelper = m_muonMgr-> cscIdHelper();
    m_tgcIdHelper = m_muonMgr-> tgcIdHelper();

    log << MSG::INFO << m_muonMgr->geometryVersion() << endreq; 
    
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

    // Retrieve muon station builder tool   -------------------------------------------------    
   if (m_muonStationTypeBuilder.retrieve().isFailure())
    {
      log << MSG::FATAL << "Failed to retrieve tool " << m_muonStationTypeBuilder << endreq;
      return StatusCode::FAILURE;
    } else
      log << MSG::INFO << "Retrieved tool " << m_muonStationTypeBuilder << endreq;

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
    // set the magnetic field 
    if (!m_magFieldMode)
        m_muonMagneticField = Trk::MagneticFieldProperties();
    else                                            
        m_muonMagneticField = Trk::MagneticFieldProperties(&(*m_magFieldTool), Trk::RealisticField);    

    log << MSG::INFO  << name() <<" initialize() successful" << endreq;    
    
  return StatusCode::SUCCESS;
}

const std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonStationBuilder::buildDetachedTrackingVolumes()
 const
{
  MsgStream log(msgSvc(), name());

  std::vector<const Trk::DetachedTrackingVolume*> mStations;

  if (m_muonMgr) { 
    // retrieve muon station prototypes from GeoModel
    const std::vector<const Trk::DetachedTrackingVolume*>* msTypes = buildDetachedTrackingVolumeTypes();
    std::vector<const Trk::DetachedTrackingVolume*>::const_iterator msTypeIter = msTypes->begin();
    
    // position MDT chambers by repeating loop over muon tree
    // link to top tree
    const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
    GeoVolumeCursor vol(top);
    while (!vol.atEnd())
    {
      const GeoVPhysVol* cv = &(*(vol.getVolume()));
      const GeoLogVol* clv = cv->getLogVol();
      std::string vname = clv->getName();
      if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && 
            ( (m_buildBarrel && vname.substr(0,1) =="B")
            ||(m_buildEndcap && vname.substr(0,1) =="E")
	    ||(m_buildCsc    && vname.substr(0,1) =="C")
	    ||(m_buildTgc    && vname.substr(0,1) =="T") ) ) {

        int etaphi = vol.getId();        // retrive eta/phi indexes
        int sign =( etaphi < 0 ) ? -1 : 1 ;
        etaphi = sign*etaphi;
        int is_mirr = etaphi/1000;
        etaphi = etaphi - is_mirr*1000;
        int eta = etaphi/100;
        int phi = etaphi - eta*100;
        eta = eta*sign;
	MuonGM::MuonStation* gmStation = m_muonMgr->getMuonStation(vname.substr(0,3),eta,phi);
        // try to retrieve 
	if ( !gmStation) {
          gmStation = m_muonMgr->getMuonStation(vname.substr(0,4),eta,phi);
        }
        // assembly ?
	if ( !gmStation) {
	  int etaphi = vol.getId();        // retrieve eta/phi indexes
          int a_etaphi = static_cast<int> (etaphi/100000);
          int sideC = static_cast<int> (a_etaphi/10000);
          a_etaphi -= sideC*10000; 
          is_mirr = static_cast<int> (a_etaphi/1000);
          a_etaphi -= is_mirr*1000; 
          eta = static_cast<int> (a_etaphi/100);
          phi = a_etaphi - eta*100;
          if (sideC) eta *=-1;
          gmStation = m_muonMgr->getMuonStation(vname.substr(0,3),eta,phi);
        }
        //
        if (!gmStation) log << MSG::WARNING << "Muon station not found! "<<vname<<","<<eta<<","<<phi  <<std::endl; 
        std::string stName = (clv->getName()).substr(0,vname.size()-8);
        if (stName.substr(0,1)=="B" && eta < 0 ) {
          stName = (clv->getName()).substr(0,vname.size()-8) + "-";
        }
        if (stName.substr(0,1)=="T" || stName.substr(0,1)=="C") {
          //std::string tgc_name = cv->getChildVol(0)->getLogVol()->getName();
          stName = vname.substr(0,4);
        }
        // loop over prototypes
        const Trk::DetachedTrackingVolume* msTV = 0;
        for (msTypeIter = msTypes->begin(); msTypeIter != msTypes->end(); ++msTypeIter) { 
          std::string msTypeName = (*msTypeIter)->name();
          if (  (stName.substr(0,1)=="T" && stName == msTypeName.substr(0,stName.size()) )
              ||(stName.substr(0,1)!="T" && stName == msTypeName) ) {
            msTV = *msTypeIter;
	    if (msTV && gmStation) {
	      HepTransform3D transf = gmStation->getTransform(); 
              Identifier stId(0);
              if (stName.substr(0,1)=="C") {
		stId = m_cscIdHelper->elementID(vname.substr(0,3),eta,phi);
              }
              // adjust eta,phi
              if (msTypeName.substr(0,1)=="C") {
                eta = 1;
		if (transf.getTranslation().z() < 0 ) eta = 0;
		double phic = transf.getTranslation().phi() + 0.1 ;  
                phi = static_cast<int> (phic<0 ? 4*phic/M_PI+8 : 4*phic/M_PI);
              } 
	      if (msTypeName.substr(0,1)=="T") {
		bool az = true;
		std::string msName = msTV->trackingVolume()->volumeName();
		if (transf.getTranslation().z() < 0 ) az = false;
		if (msName.substr(7,2)=="01") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="02") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="03") eta = az ? 6 : 3;
		if (msName.substr(7,2)=="04") eta = az ? 7 : 2;
		if (msName.substr(7,2)=="05") eta = az ? 8 : 1;
		if (msName.substr(7,2)=="06") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="07") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="08") eta = az ? 6 : 3;
		if (msName.substr(7,2)=="09") eta = az ? 7 : 2;
		if (msName.substr(7,2)=="10") eta = az ? 8 : 1;
		if (msName.substr(7,2)=="11") eta = az ? 9 : 0;
		if (msName.substr(7,2)=="12") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="13") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="14") eta = az ? 6 : 3;
		if (msName.substr(7,2)=="15") eta = az ? 7 : 2;
		if (msName.substr(7,2)=="16") eta = az ? 8 : 1;
		if (msName.substr(7,2)=="17") eta = az ? 9 : 0;
		if (msName.substr(7,2)=="18") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="19") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="20") eta = az ? 5 : 4;
		if (msName.substr(7,2)=="21") eta = az ? 5 : 4;
	      }     
              if (stName.substr(0,1)=="T") {
		int etaSt = eta - 4;
		if (eta < 5) etaSt = eta - 5; 
		double phic = transf.getTranslation().phi(); 
		if (msTypeName.substr(2,1)=="E" && msTypeName.substr(0,3)!="T4E")
		  phi = static_cast<int> (phic<0 ? 24*phic/M_PI+48 : 24*phic/M_PI);
		else
                  phi = static_cast<int> (phic<0 ? 12*phic/M_PI+24 : 12*phic/M_PI);
                phi++;
                stId = m_tgcIdHelper->elementID(vname.substr(0,3),etaSt,phi);
              } else if (stName.substr(0,3)=="BML") {
		stId = m_rpcIdHelper->elementID(vname.substr(0,3),eta,phi,1);
              } else if (stName.substr(0,1)!="C" ) {
		stId = m_mdtIdHelper->elementID(vname.substr(0,3),eta,phi);
              }
              if (!stId) log << MSG::WARNING << "identifier of the station not found:"<<vname <<","<<eta<<","<<phi<<endreq;
              unsigned int iD = stId;
	      // clone station from prototype
	      const Trk::DetachedTrackingVolume* newStat = msTV->clone(vname,transf);
              // identify layer representation
              newStat->layerRepresentation()->setLayerType(iD);
              // resolved stations only:
              if (m_resolveActiveLayers) {
		// glue components
		glueComponents(newStat);
		// identify layers
		identifyLayers(newStat,eta,phi); 
	      } 
	      mStations.push_back(newStat);
            }
	  }
	}  
	if (!msTV)  log << MSG::DEBUG  << name() <<" this station has no prototype: " << vname << endreq;    
      }
      vol.next();      
    }
    // clean up prototypes
    for (unsigned int it = 0; it < msTypes->size(); it++) delete (*msTypes)[it];
    delete msTypes;   
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonStations=new std::vector<const Trk::DetachedTrackingVolume*>(mStations);

  //
  log << MSG::INFO  << name() << "returns " << (*muonStations).size() << " stations" << endreq;	 
  return muonStations; 
}

const std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonStationBuilder::buildDetachedTrackingVolumeTypes() const 
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building station types" << endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////
      std::vector<const Trk::DetachedTrackingVolume*> stations;

   if (m_muonMgr){

      // link to top tree
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      GeoVolumeCursor vol (top);
      while (!vol.atEnd())
      {
        const GeoVPhysVol* cv = &(*(vol.getVolume()));
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
	if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && 
            ( (m_buildBarrel && vname.substr(0,1) =="B")
            ||(m_buildEndcap && vname.substr(0,1) =="E")
            ||(m_buildCsc && vname.substr(0,1) =="C")
	    ||(m_buildTgc && vname.substr(0,1) =="T") ) )
	{
          int etaphi = vol.getId();        // retrieve eta/phi indexes
          int sign  = ( etaphi < 0 ) ? -1 : 1;
          etaphi = sign*etaphi;
          int is_mirr = etaphi/1000;
          etaphi = etaphi - is_mirr*1000;
          int eta = etaphi/100;
          int phi = etaphi - eta*100;
          eta = eta*sign;
	  // std::cout << vname.substr(0,3) <<","<<is_mirr<<","<<eta<<","<<phi<<std::endl;
	  MuonGM::MuonStation* gmStation = m_muonMgr->getMuonStation(vname.substr(0,3),eta,phi);
	  if ( !gmStation) {
            gmStation = m_muonMgr->getMuonStation(vname.substr(0,4),eta,phi);
          }
	  // assembly ?
	  if ( !gmStation) {
	    int etaphi = vol.getId();        // retrieve eta/phi indexes
	    int a_etaphi = static_cast<int> (etaphi/100000);
	    int sideC = static_cast<int> (a_etaphi/10000);
	    a_etaphi -= sideC*10000; 
	    is_mirr = static_cast<int> (a_etaphi/1000);
	    a_etaphi -= is_mirr*1000; 
	    eta = static_cast<int> (a_etaphi/100);
	    phi = a_etaphi - eta*100;
	    if (sideC) eta *=-1;
	    gmStation = m_muonMgr->getMuonStation(vname.substr(0,3),eta,phi);
	  }
          //
          std::string name = (clv->getName()).substr(0,vname.size()-8);
          // is this station known ?        
          // if TGC station, look for 1 component instead
          if (name.substr(0,1)=="T") {
            std::string tgc_name = cv->getChildVol(0)->getLogVol()->getName();
            name = tgc_name;
          }
          if (name.substr(0,1)=="B" && eta < 0 ) {
            name = (clv->getName()).substr(0,vname.size()-8) + "-";
          }
          unsigned is=0; 
          for (unsigned in=0; in< stations.size(); in++) 
          {
            if (stations[in]!=0 && name == stations[in]->name()) is++;
          }
          // if (is!=0 ) std::cout << "prototype exists" << std::endl;
          // if (is==0 ) std::cout << "station shape type:"<< name<< ","<<clv->getShape()->type()<<std::endl; 
          if (is==0 )          
          {
            log << MSG::DEBUG <<" new station type " << name << "," << clv->getShape()->type() << endreq;    
            log << MSG::DEBUG <<" prototype built from eta, phi:" << eta << "," << phi << endreq;    
             
            if (name.substr(0,2)=="CS" || name.substr(0,1)=="T") {
              if (m_muonStationTypeBuilder) {
                if (name.substr(0,2)=="CS") { 
                  const Trk::TrackingVolume* csc_station = m_muonStationTypeBuilder->processCscStation(cv, name);   
		  // create layer representation
		  std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>*> layerRepr =
		    m_muonStationTypeBuilder->createLayerRepresentation(csc_station);
		  // create prototype as detached tracking volume
		  const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(name,csc_station,layerRepr.first,layerRepr.second);
                  if (!m_resolveActiveLayers) typeStat->trackingVolume()->clear();
		  stations.push_back(typeStat); 
                } else {
                  std::vector<const Trk::TrackingVolume*> tgc_stations = m_muonStationTypeBuilder->processTgcStation(cv);   
                  for (unsigned int i=0;i<tgc_stations.size();i++) {
		    // create layer representation
		    std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>*> layerRepr =
		      m_muonStationTypeBuilder->createLayerRepresentation(tgc_stations[i]);
		    // create prototype as detached tracking volume
		    const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(name,tgc_stations[i],layerRepr.first,layerRepr.second);
		    if (!m_resolveActiveLayers) typeStat->trackingVolume()->clear();
                    stations.push_back(typeStat); 
                  }
                }
              }
            } else {    
              
              const GeoShape* shapeS = clv->getShape();
              while ( shapeS->type() != "Trd" ){
		if (shapeS->type()=="Shift") {
		  const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (shapeS); 
                  shapeS = shift->getOp();
		} else if (shapeS->type()=="Subtraction") {
		  const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (shapeS); 
                  shapeS = sub->getOpA();
		} else if (shapeS->type()=="Union") {
		  const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion*> (shapeS); 
                  shapeS = uni->getOpA();
		} else {
		  log << MSG::WARNING << "unexpected station shape ? "<< shapeS->type() << ", station not built" << std::endl;
		  break; 
		}
	      }
	      const GeoTrd* trd=dynamic_cast<const GeoTrd*> ( shapeS );

	      double halfX1=0.;
	      double halfX2=0.;
	      double halfY1=0.;
	      double halfY2=0.;
	      double halfZ=0.;   
	      if (trd) {
		//
		halfX1 = trd->getXHalfLength1();
		halfX2 = trd->getXHalfLength2();
		halfY1 = trd->getYHalfLength1();
		halfY2 = trd->getYHalfLength2();
		halfZ  = trd->getZHalfLength();              
	      
		// define enveloping volume
		const Trk::TrackingVolumeArray* confinedVolumes = 0; 
		Trk::Volume* envelope = 0;
		std::string shape = "Trd";
		if (halfX1==halfX2 && halfY1==halfY2) shape = "Box";
		if (shape=="Box") { 
		  Trk::CuboidVolumeBounds* envBounds = new Trk::CuboidVolumeBounds(halfX1,halfY1,halfZ);
		  // station components
		  if (m_muonStationTypeBuilder) confinedVolumes = 
						  m_muonStationTypeBuilder->processBoxStationComponents(cv,envBounds); 
		  // enveloping volume
		  envelope= new Trk::Volume(new HepTransform3D(),envBounds);
		}
		if (shape=="Trd") {
		  Trk::TrapezoidVolumeBounds* envBounds = 0;
		  HepTransform3D* transf =new HepTransform3D(); 
		  if (halfY1==halfY2) {
		    envBounds = new Trk::TrapezoidVolumeBounds(halfX1,halfX2,halfY1,halfZ);
		    std::cout << "CAUTION!!!: this trapezoid volume does not require XY -> YZ switch" << std::endl;
		  }
		  if (halfY1!=halfY2 && halfX1 == halfX2 ) {
		    delete transf;
		    transf = new HepTransform3D( HepRotateY3D(90*deg)* HepRotateZ3D(90*deg) );
		    envBounds = new Trk::TrapezoidVolumeBounds(halfY1,halfY2,halfZ,halfX1); 
		  }
		  if (halfX1!=halfX2 && halfY1!=halfY2 ) std::cout << "station envelope arbitrary trapezoid?" << std::endl;
		  // station components
		  if (m_muonStationTypeBuilder) confinedVolumes = 
						  m_muonStationTypeBuilder->processTrdStationComponents(cv,envBounds); 
		  // enveloping volume
		  envelope= new Trk::Volume(transf,envBounds);
		}
		// hack to verify BI/BM stations
		//  if (name.substr(0,2)=="BI") name = gmStation->getKey();
		//  if (name.substr(0,2)=="BM") name = gmStation->getKey();
		
		// ready to build the station prototype
		const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *envelope,
									     m_muonMaterial,
									     m_muonMagneticField,
									     0,confinedVolumes,
									     name);         
		delete envelope; 
		
		// identify prototype
		if (m_resolveActiveLayers && ( name.substr(0,1)=="B" || name.substr(0,1)=="E") ) identifyPrototype(newType,eta,phi,gmStation->getTransform());
		
		// create layer representation
		std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>*> layerRepr 
		  = m_muonStationTypeBuilder->createLayerRepresentation(newType);
		
		// create prototype as detached tracking volume
		const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(name,newType,layerRepr.first,layerRepr.second);
		
		if (!m_resolveActiveLayers) typeStat->trackingVolume()->clear();
		stations.push_back(typeStat); 
	      }
	    }
	  } // end new station type 
	} // end if "Shift" (station)
	vol.next();      
      }      
      log << MSG::INFO  << name() << stations.size() <<" station prototypes built " << endreq;    
   }
   
///////////////////////////////////////////////////////////////////////////////////////
   const std::vector<const Trk::DetachedTrackingVolume*>* mStations = new std::vector<const Trk::DetachedTrackingVolume*>(stations); 
   return mStations;  
}

// finalize
StatusCode Muon::MuonStationBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    return StatusCode::SUCCESS;
}
//
void Muon::MuonStationBuilder::glueComponents(const Trk::DetachedTrackingVolume* stat) const
{
   const Trk::TrackingVolumeArray* volArray = stat->trackingVolume()->confinedVolumes();
   if (volArray) {
     if (volArray->arrayObjectsNumber() > 1) {
       const std::vector< const Trk::TrackingVolume* > components = volArray->arrayObjects();
       const Trk::BinUtility1DX* binUtilityX = dynamic_cast<const Trk::BinUtility1DX*>(volArray->binUtility()); 
       const Trk::CuboidVolumeBounds* cubVolBounds = 
	 dynamic_cast<const Trk::CuboidVolumeBounds* >( &(components[0]->volumeBounds())); 
       /*
       const Trk::TrapezoidVolumeBounds* trdVolBounds = 
	 dynamic_cast<const Trk::TrapezoidVolumeBounds* >( &(components[0]->volumeBounds())); 
       const Trk::DoubleTrapezoidVolumeBounds* dtrdVolBounds = 
	 dynamic_cast<const Trk::DoubleTrapezoidVolumeBounds* >( &(components[0]->volumeBounds())); 
       */
       // identify 'lower' and 'upper' boundary surface
       Trk::BoundarySurfaceFace low = Trk::negativeFaceXY;
       Trk::BoundarySurfaceFace up = Trk::positiveFaceXY;

       // rectangular station in x ordering (MDT barrel)
       if ( cubVolBounds && binUtilityX ) {
	 low = Trk::negativeFaceYZ;
	 up  = Trk::positiveFaceYZ;
       }
       
       /*
       // rotated trapezoid station in x ordering (MDT endcap)
       if ( trdVolBounds && binUtilityX ) {
	 low = Trk::negativeFaceXY;
	 up  = Trk::positiveFaceXY;
       }
       // rotated diamond station in x ordering (MDT endcap)
       if ( dtrdVolBounds && binUtilityX ) {
	 low = Trk::negativeFaceXY;
	 up  = Trk::positiveFaceXY;
       }
       */

       if (low>=0 && up >=0) {
         // glue volumes 
         for (unsigned int i = 0; i< components.size()-1; i++) {
           m_trackingVolumeHelper -> glueTrackingVolumes(*(components[i]),up,
                                                         *(components[i+1]),low );
         }
       }
     }
   }
}


void Muon::MuonStationBuilder::identifyLayers(const Trk::DetachedTrackingVolume* station, int eta, int phi ) const
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG  << name() <<" identifying layers " << endreq;    

  std::string stationName = station->trackingVolume()->volumeName();
  log << MSG::DEBUG  << " in station " << station->name() << endreq;    

  
  if (stationName.substr(0,1)=="C") { 
    int st = stationName.substr(0,3)=="CSS" ? 0 : 1;
    const MuonGM::CscReadoutElement* cscRE = m_muonMgr->getCscReadoutElement(st,eta,phi,0);    
    int cLay = cscRE ? 0 : 1;
    if (!cscRE) cscRE = m_muonMgr->getCscReadoutElement(st,eta,phi,cLay);
    if (cscRE) {
      for (int gasgap = 0; gasgap < cscRE->Ngasgaps(); gasgap++) {
	Identifier idi = m_cscIdHelper->channelID(cscRE->identify(),cscRE->ChamberLayer(),gasgap+1,0,1);   
        const Trk::PlaneSurface* stripSurf = dynamic_cast<const Trk::PlaneSurface*> (&(cscRE->surface(idi)));
        const HepPoint3D gpi = stripSurf->center();
        const Trk::TrackingVolume* assocVol = station->trackingVolume()->associatedSubVolume(gpi);
        const Trk::Layer* assocLay = 0;
        if (assocVol) assocLay = assocVol->associatedLayer(gpi);
        unsigned int iD = idi;
        if (assocVol && assocLay) assocLay->setLayerType(iD);        
	if (assocLay) assocLay->setRef((assocLay->surfaceRepresentation().transform().inverse()*gpi)[1]);
      }
    } else {
      log<< MSG::DEBUG << "cscRE not found:"<<st<<","<<eta<<","<<phi<<std::endl;
    }
  }
  

  if (stationName.substr(0,1)=="T") {
    int st = 7;
    if (stationName.substr(0,3)=="T1F") {
      st = 0;
    } else if (stationName.substr(0,3)=="T1E") {
      st = 1;
    } else if (stationName.substr(0,3)=="T2F") {
      st = 2;
    } else if (stationName.substr(0,3)=="T2E") {
      st = 3;
    } else if (stationName.substr(0,3)=="T3F") {
      st = 4;
    } else if (stationName.substr(0,3)=="T3E") {
      st = 5;
    } else if (stationName.substr(0,3)=="T4F") {
      st = 6;
    }
  
    const MuonGM::TgcReadoutElement* tgc = m_muonMgr->getTgcReadoutElement(st,eta,phi-1);
    
    if (!tgc || !(station->trackingVolume()->inside(tgc->center(),0.)) || 
        (station->trackingVolume()->center()-tgc->center()).mag()>0.1 ) {
      unsigned int phit=0;
      while ( phit<48 ) {
	const MuonGM::TgcReadoutElement* tgct = m_muonMgr->getTgcReadoutElement(st,eta,phit);
        if (tgct && station->trackingVolume()->inside(tgct->center(),0.)) {
          tgc = tgct;
          phi = phit;
          // update station identity
          Identifier oldId(station->layerRepresentation()->layerType());
          int stationName = m_tgcIdHelper->stationName(oldId);
          int stationEta  = m_tgcIdHelper->stationEta(oldId);
          Identifier stId = m_tgcIdHelper->elementID(stationName,stationEta,phi);
          station->layerRepresentation()->setLayerType(stId);
          break;
        }
	phit++;  
      }
    }
    
    if (tgc) {
      int etaSt = eta - 4;
      if (eta < 5) etaSt = eta - 5; 
      Identifier wireId  = m_tgcIdHelper->channelID(stationName.substr(0,3),etaSt,phi,1,0,1);
      const HepPoint3D gp = tgc->channelPos(wireId);
      const Trk::TrackingVolume* assocVol = station->trackingVolume()->associatedSubVolume(gp);
      if (!assocVol) log << MSG::DEBUG << "wrong tgcROE?" << stationName <<"," << eta <<"," << phi << endreq;
      if (assocVol && assocVol->confinedLayers()) {
	const std::vector<const Trk::Layer*> layers = assocVol->confinedLayers()->arrayObjects();           
	for (unsigned int il=0;il<layers.size();il++) {
	  Identifier wireId  = m_tgcIdHelper->channelID(stationName.substr(0,3),etaSt,phi,il+1,0,1);
	  unsigned int id = wireId;
          layers[il]->setLayerType(id); 
	  // strip plane surface
	  const Trk::PlaneSurface* stripSurf = dynamic_cast<const Trk::PlaneSurface*> (&(tgc->surface(wireId)));
	  if ( (layers[il]->surfaceRepresentation().transform().inverse()*stripSurf->center()).mag()>0.001)   
	    log<< MSG::INFO << "TGC strip plane shifted:"<<st<<","<<eta<<","<<phi<<":" <<
	      layers[il]->surfaceRepresentation().transform().inverse()*stripSurf->center()<< std::endl;
          /*
	  Identifier s1 = m_tgcIdHelper->channelID(stationName.substr(0,3),etaSt,phi,il+1,1,1); 
	  HepPoint3D lw1 = (layers[il]->surfaceRepresentation().transform().inverse()) * (tgc->channelPos(wireId));
	  HepPoint3D ls1 = (layers[il]->surfaceRepresentation().transform().inverse()) * (tgc->channelPos(s1));
          int wireRe = 1000*lw1[Trk::locY];                       
	  layers[il]->setRef( 10e4+ls1[Trk::locX] + 10e5*(wireRe+10e6) );                      
          */
        }
      }
    } else {
      log << MSG::WARNING << name() << "tgcROE not found for :" << stationName <<","<<eta<<","<<phi<<endreq;         
      /*
      for (unsigned int e=0; e < 10 ; e++) {
	for (unsigned int p=0; p < 48; p++) {
          tgc = m_muonMgr->getTgcReadoutElement(st,e,p);
          if (tgc) std::cout << "tgc exists:"<<e<<"," << p <<"," << station->trackingVolume()->inside(tgc->center(),0.) << std::endl;
        }
      }
      */
    }
  }

  if (stationName.substr(0,1)=="B" || stationName.substr(0,1)=="E" ) { 
    // recalculate id
    Identifier stId(station->layerRepresentation()->layerType());
    int nameIndex = m_mdtIdHelper->stationNameIndex( stationName.substr(0,3) ); 
    if (station->trackingVolume()->confinedVolumes()) {
      const std::vector<const Trk::TrackingVolume*> cVols = station->trackingVolume()->confinedVolumes()->arrayObjects();
      for (unsigned int i=0; i<cVols.size() ; i++) {
	if (cVols[i]->confinedLayers()) {
	  const std::vector<const Trk::Layer*> cLays = cVols[i]->confinedLayers()->arrayObjects();
	  for (unsigned int il=0; il<cLays.size() ; il++) {
	    Identifier id(cLays[il]->layerType());
	    if (id>0 && m_mdtIdHelper->is_mdt(id)) {
              Identifier newId = m_mdtIdHelper->channelID(nameIndex,eta,phi,
							  m_mdtIdHelper->multilayer(id),m_mdtIdHelper->tubeLayer(id),m_mdtIdHelper->tube(id));
              int newid = newId;
              cLays[il]->setLayerType(newid); 
            }
	  }
        }
	if (cVols[i]->confinedArbitraryLayers()) {
	  const std::vector<const Trk::Layer*>* cLays = cVols[i]->confinedArbitraryLayers();
	  for (unsigned int il=0; il<cLays->size() ; il++) {
	    Identifier id((*cLays)[il]->layerType());
	    if (id>0 && m_rpcIdHelper->is_rpc(id)) {
              Identifier newId = m_rpcIdHelper->channelID(nameIndex,eta,phi,m_rpcIdHelper->doubletR(id),
							  m_rpcIdHelper->doubletZ(id),m_rpcIdHelper->doubletPhi(id),m_rpcIdHelper->gasGap(id),
							  m_rpcIdHelper->measuresPhi(id),m_rpcIdHelper->strip(id));
              int newid = newId;
              (*cLays)[il]->setLayerType(newid);  
            }
	  }
	}
      }
    }
  }

  /*
  if (stationName.substr(0,1)=="B" || stationName.substr(0,1)=="E" ) { 

    // MDT
    int nameIndex = m_mdtIdHelper->stationNameIndex( stationName.substr(0,3) ); 
    int nameIndexC = nameIndex;
    if (stationName.substr(0,3)=="EIS") nameIndexC = 22; 
    if (stationName.substr(0,3)=="BIM") nameIndexC = 23; 
    for (int multi = 0; multi < 2; multi++ ) {
      const MuonGM::MdtReadoutElement* multilayer = m_muonMgr->getMdtReadoutElement(nameIndexC,eta+8,phi-1,multi);
      if (multilayer) {
        const Trk::TrackingVolume* assocVol = station->trackingVolume()->associatedSubVolume(multilayer->center());
        if (!assocVol ) log << MSG::WARNING << "valid multilayer outside station:" << stationName <<endreq;
        if (assocVol) {
	  int nLayers = multilayer->getNLayers();
	  for (int layer =1; layer <= nLayers ; layer++) {
	    Identifier id = m_mdtIdHelper->channelID(nameIndex,eta,phi,multi+1,layer,1);           
	    if (id>0) {
	      // retrieve associated layer
	      HepPoint3D gp = multilayer->tubePos(id);
	      const Trk::Layer* assocLay = assocVol->associatedLayer(gp);
	      unsigned int iD = id;
	      if (assocLay && assocLay->layerType()!= iD) {
		std::cout <<"ERROR IDENTIFICATION:"<<stationName<<","<<assocLay->layerType()<<","<<id<<std::endl;
                Identifier newId(assocLay->layerType());  
		std::cout << m_mdtIdHelper->elementID(newId) <<","<< m_mdtIdHelper->elementID(id)<<std::endl;
              }
	      if (assocLay) assocLay->setLayerType(iD); 
	      if (assocLay) {
		//const Trk::LocalPosition* locPos = (assocLay->surfaceRepresentation()).globalToLocal(gp,0.001);
                Trk::GlobalPosition locPos = assocLay->surfaceRepresentation().transform().inverse()*gp;
		if (fabs(assocLay->getRef()-locPos[1])>0.01) 
		  std::cout << "ERROR REFERENCE:"<< stationName<<","<<assocLay->getRef()<<","<<locPos[1]<<std::endl;
		if (fabs(locPos[2])>0.01) 
		  std::cout << "ERROR Alignment:"<< stationName<<","<<assocLay->getRef()<<","<<locPos[2]<<std::endl;
		assocLay->setRef(locPos[1]);
	      }
	    } 
          }
        }
      }     
    }
    
    // RPC ?
    const Trk::BinnedArray< Trk::TrackingVolume >* confinedVolumes = station->trackingVolume()->confinedVolumes();
    if (confinedVolumes){
      const std::vector<const Trk::TrackingVolume*>& vols = confinedVolumes->arrayObjects();
      for (unsigned int iv=0;iv<vols.size();iv++) if (vols[iv]->volumeName() == "RPC") {
        // for active layers do a search of associated ROE
        const std::vector<const Trk::Layer*>* layers = vols[iv]->confinedArbitraryLayers();
        int nameIndex = m_rpcIdHelper->stationNameIndex( stationName.substr(0,3) ); 
        // loop over doubletR, doubletZ 
	for (int doubletR = 0; doubletR < 2; doubletR++ ) {
	for (int doubletZ = 0; doubletZ < 3; doubletZ++ ) {
	for (int doubletPhi = 0; doubletPhi < 1; doubletPhi++ ) {
	  const MuonGM::RpcReadoutElement* rpc = m_muonMgr->getRpcReadoutElement(nameIndex-2,eta+8,phi-1,doubletR,doubletZ);
	  if (rpc) {
            if (doubletZ < rpc->getDoubletZ() ) {  
	      for (int gasGap=0; gasGap<2; gasGap++) {
		Identifier etaId = m_rpcIdHelper->channelID(nameIndex,eta,phi,
							    doubletR+1,doubletZ+1,doubletPhi+1,gasGap+1,0,1); 
		Identifier phiId = m_rpcIdHelper->channelID(nameIndex,eta,phi,
							    doubletR+1,doubletZ+1,doubletPhi+1,gasGap+1,1,1); 
		if (m_rpcIdHelper->valid(etaId)){
		  for (unsigned int il=0;il<layers->size();il++) {
		    if ((*layers)[il]->layerType() != 0 && (*layers)[il]->isOnLayer(rpc->stripPos(etaId)) ) {
		      unsigned int id = etaId;
		      if ((*layers)[il]->layerType()!= id) {
                        std::cout <<"ERROR IDENTIFICATION:RPC:"<<stationName<<","<<(*layers)[il]->layerType()<<","<<id<<std::endl;
                        std::cout <<"ERROR IDENTIFICATION:RPC:"<<il<<","<<(*layers)[il]->surfaceRepresentation().center() <<std::endl;
			Identifier newId((*layers)[il]->layerType());  
			std::cout << m_rpcIdHelper->elementID(newId) <<","<< m_rpcIdHelper->elementID(etaId)<<std::endl;
			std::cout << m_rpcIdHelper->doubletR(newId) <<","<< m_rpcIdHelper->doubletR(etaId)<<std::endl;
                        std::cout << m_rpcIdHelper->doubletZ(newId)<<","<<m_rpcIdHelper->doubletPhi(newId)<<","<<m_rpcIdHelper->gasGap(newId)
				  <<","<<m_rpcIdHelper->measuresPhi(newId)<<","<<m_rpcIdHelper->strip(newId)<<std::endl;
                        std::cout << m_rpcIdHelper->doubletZ(etaId)<<","<<m_rpcIdHelper->doubletPhi(etaId)<<","<<m_rpcIdHelper->gasGap(etaId)
				  <<","<<m_rpcIdHelper->measuresPhi(etaId)<<","<<m_rpcIdHelper->strip(etaId)<<std::endl;
                      } else {
                        std::cout <<"RPC IDENTIFICATION OK:"<<stationName<<","<<(*layers)[il]->layerType()<<","<<id<<std::endl;
                        std::cout <<"RPC IDENTIFICATION OK:"<<il<<","<<(*layers)[il]->surfaceRepresentation().center() <<std::endl;
                      }
		      (*layers)[il]->setLayerType(id);
                      //turn eta position into integer to truncate
		      HepPoint3D locPosEta = ((*layers)[il]->surfaceRepresentation().transform().inverse()) * (rpc->stripPos(etaId));
		      HepPoint3D locPosPhi = ((*layers)[il]->surfaceRepresentation().transform().inverse()) * (rpc->stripPos(phiId));
                      int etaRefi = int(1000*locPosEta[Trk::locY]);                       
		      (*layers)[il]->setRef( 10e4+locPosPhi[Trk::locX] + 10e5*(etaRefi+10e6));
		    } 
		  }
		}
	 }}}}}}                  
      }
    }
  } 
  */ 

  // by now, all the layers should be identified - verify
  if (station->trackingVolume()->confinedVolumes()) {
    const std::vector<const Trk::TrackingVolume*> cVols = station->trackingVolume()->confinedVolumes()->arrayObjects();
    for (unsigned int i=0; i<cVols.size() ; i++) {
      if (cVols[i]->confinedLayers()) {
        const std::vector<const Trk::Layer*> cLays = cVols[i]->confinedLayers()->arrayObjects();
        for (unsigned int il=0; il<cLays.size() ; il++) {
          Identifier id(cLays[il]->layerType());
          if (id==1) log << MSG::DEBUG << station->name()<<","<< cVols[i]->volumeName()<<", unidentified active layer:"<<il<<endreq;
        }
      }
      if (cVols[i]->confinedArbitraryLayers()) {
        const std::vector<const Trk::Layer*>* cLays = cVols[i]->confinedArbitraryLayers();
        for (unsigned int il=0; il<cLays->size() ; il++) {
          Identifier id((*cLays)[il]->layerType());
          if (id==1) log << MSG::DEBUG << station->name()<<","<< cVols[i]->volumeName()<<", unidentified active layer:"<<il<<endreq;
        }
      }
    }
  } 
  if (station->trackingVolume()->confinedLayers()) {
    const std::vector<const Trk::Layer*> cLays = station->trackingVolume()->confinedLayers()->arrayObjects();
    for (unsigned int il=0; il<cLays.size() ; il++) {
      Identifier id(cLays[il]->layerType());
      if (id==1) log << MSG::DEBUG << station->name()<<","<< station->name()<<", unidentified active layer:"<<il<<endreq;
    }
  }
  // end identification check

}

void Muon::MuonStationBuilder::identifyPrototype(const Trk::TrackingVolume* station, int eta, int phi, HepTransform3D transf ) const
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG  << name() <<" identifying prototype " << endreq;    

  std::string stationName = station->volumeName();
  log << MSG::DEBUG  << " for station " << stationName << endreq;    

  if (stationName.substr(0,1)=="B" || stationName.substr(0,1)=="E" ) { 
    // MDT
    int nameIndex = m_mdtIdHelper->stationNameIndex( stationName.substr(0,3) ); 
    int nameIndexC = nameIndex;
    if (stationName.substr(0,3)=="EIS") nameIndexC = 22; 
    if (stationName.substr(0,3)=="BIM") nameIndexC = 23; 
    for (int multi = 0; multi < 2; multi++ ) {
      const MuonGM::MdtReadoutElement* multilayer = m_muonMgr->getMdtReadoutElement(nameIndexC,eta+8,phi-1,multi);
      if (multilayer) {
        const Trk::TrackingVolume* assocVol = station->associatedSubVolume(transf.inverse()*multilayer->center());
        if (!assocVol ) log << MSG::WARNING << "valid multilayer outside station:" << stationName <<endreq;
        if (assocVol) {
	  int nLayers = multilayer->getNLayers();
	  for (int layer =1; layer <= nLayers ; layer++) {
	    Identifier id = m_mdtIdHelper->channelID(nameIndex,eta,phi,multi+1,layer,1);           
	    if (id>0) {
	      // retrieve associated layer
	      HepPoint3D gp = multilayer->tubePos(id);
	      const Trk::Layer* assocLay = assocVol->associatedLayer(transf.inverse()*gp);
	      unsigned int iD = id;
	      if (assocLay) assocLay->setLayerType(iD); 
	      //if (assocLay) {
	      //	const Trk::LocalPosition* locPos = (assocLay->surfaceRepresentation()).globalToLocal(gp,0.001);
	      //	if (fabs(assocLay->getRef()-(*locPos)[Trk::locY])>0.01) 
	      //  std::cout << "ERROR REFERENCE:"<< stationName<<","<<assocLay->getRef()<<","<<(*locPos)[Trk::locY]<<std::endl;
	      //	if (locPos) assocLay->setRef((*locPos)[Trk::locY]);
	      //}
	    } 
          }
        }
      }     
    }
    
    // RPC ?
    const Trk::BinnedArray< Trk::TrackingVolume >* confinedVolumes = station->confinedVolumes();
    if (confinedVolumes){
      const std::vector<const Trk::TrackingVolume*>& vols = confinedVolumes->arrayObjects();
      for (unsigned int iv=0;iv<vols.size();iv++) if (vols[iv]->volumeName() == "RPC") {
        // for active layers do a search of associated ROE
        const std::vector<const Trk::Layer*>* layers = vols[iv]->confinedArbitraryLayers();
        int nameIndex = m_rpcIdHelper->stationNameIndex( stationName.substr(0,3) ); 
        // loop over doubletR, doubletZ 
	for (int doubletR = 0; doubletR < 2; doubletR++ ) {
	for (int doubletZ = 0; doubletZ < 3; doubletZ++ ) {
	for (int doubletPhi = 0; doubletPhi < 1; doubletPhi++ ) {
	  const MuonGM::RpcReadoutElement* rpc = m_muonMgr->getRpcReadoutElement(nameIndex-2,eta+8,phi-1,doubletR,doubletZ);
	  if (rpc) {
            if (doubletZ < rpc->getDoubletZ() ) {  
	      for (int gasGap=0; gasGap<2; gasGap++) {
		Identifier etaId = m_rpcIdHelper->channelID(nameIndex,eta,phi,
							    doubletR+1,doubletZ+1,doubletPhi+1,gasGap+1,0,1); 
		//Identifier phiId = m_rpcIdHelper->channelID(nameIndex,eta,phi,
		//					    doubletR+1,doubletZ+1,doubletPhi+1,gasGap+1,1,1); 
		if (1/*m_rpcIdHelper->valid(etaId)*/){
		  for (unsigned int il=0;il<layers->size();il++) {
		    if ((*layers)[il]->layerType() != 0 && (*layers)[il]->isOnLayer(transf.inverse()*rpc->stripPos(etaId)) ) {
                      const Trk::GlobalPosition locPos1 = (*layers)[il]->surfaceRepresentation().transform().inverse()*transf.inverse()*rpc->stripPos(etaId);
                      const Trk::GlobalPosition locPos2 = rpc->surface(etaId).transform().inverse()*rpc->stripPos(etaId);
                      double swap = ( fabs( locPos1[1] - locPos2[0] ) > 0.001 ) ? 20000. : 0. ;
		      unsigned int id = etaId;
		      (*layers)[il]->setLayerType(id);
                      const Trk::GlobalPosition locPos = (*layers)[il]->surfaceRepresentation().transform().inverse()
			*transf.inverse()*rpc->surface(etaId).center(); 
		      (*layers)[il]->setRef(swap + locPos[0]);
		      //std::cout <<"identifying RPC:"<<stationName<<","<<iv<<","<<il<<":"<<id <<std::endl;
                      //turn eta position into integer to truncate
		      //HepPoint3D locPosEta = ((*layers)[il]->surfaceRepresentation().transform().inverse()) * (rpc->stripPos(etaId));
		      //HepPoint3D locPosPhi = ((*layers)[il]->surfaceRepresentation().transform().inverse()) * (rpc->stripPos(phiId));
                      //int etaRefi = int(1000*locPosEta[Trk::locY]);                       
		      //(*layers)[il]->setRef( 10e4+locPosPhi[Trk::locX] + 10e5*(etaRefi+10e6));
		    } 
		  }
		}
	 }}}}}}                  
      }
    }
  } 

  // by now, all the layers should be identified - verify
  if (station->confinedVolumes()) {
    const std::vector<const Trk::TrackingVolume*> cVols = station->confinedVolumes()->arrayObjects();
    for (unsigned int i=0; i<cVols.size() ; i++) {
      if (cVols[i]->confinedLayers()) {
        const std::vector<const Trk::Layer*> cLays = cVols[i]->confinedLayers()->arrayObjects();
        for (unsigned int il=0; il<cLays.size() ; il++) {
          Identifier id(cLays[il]->layerType());
          if (id==1) log << MSG::DEBUG << station->volumeName()<<","<< cVols[i]->volumeName()<<", unidentified active layer:"<<il<<endreq;
        }
      }
      if (cVols[i]->confinedArbitraryLayers()) {
        const std::vector<const Trk::Layer*>* cLays = cVols[i]->confinedArbitraryLayers();
        for (unsigned int il=0; il<cLays->size() ; il++) {
          Identifier id((*cLays)[il]->layerType());
          if (id==1) log << MSG::DEBUG << station->volumeName()<<","<< cVols[i]->volumeName()<<", unidentified active layer:"<<il<<endreq;
        }
      }
    }
  } 
  if (station->confinedLayers()) {
    const std::vector<const Trk::Layer*> cLays = station->confinedLayers()->arrayObjects();
    for (unsigned int il=0; il<cLays.size() ; il++) {
      Identifier id(cLays[il]->layerType());
      if (id==1) log << MSG::DEBUG << station->volumeName()<<","<< station->volumeName()<<", unidentified active layer:"<<il<<endreq;
    }
  }
  // end identification check
}

