///////////////////////////////////////////////////////////////////
// MuonStationBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonStationBuilder.h"
#include "MuonTrackingGeometry/MuonStationTypeBuilder.h"
//MuonSpectrometer include
#include "MuonGeoModel/MuonDetectorManager.h"
#include "MuonGeoModel/MuonStation.h"
#include "MuonGeoModel/MdtReadoutElement.h"
#include "MuonGeoModel/RpcReadoutElement.h"
#include "MuonGeoModel/CscReadoutElement.h"
#include "MuonGeoModel/TgcReadoutElement.h"
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

// constructor
Muon::MuonStationBuilder::MuonStationBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgTool(t,n,p),
  m_muonMgrLocation("MuonMgr"),
  m_magFieldTool(0),
  m_magFieldToolName("Trk::MagneticFieldTool"),
  m_magFieldToolInstanceName("ATLAS_TrackingMagFieldTool"),
  m_muonStationTypeBuilder(0),
  m_muonStationTypeBuilderName("Muon::MuonStationTypeBuilder"),
  m_muonStationTypeBuilderInstanceName("MuonStationTypeBuilder"),
  m_trackingVolumeHelper(0),
  m_trackingVolumeHelperName("Trk::TrackingVolumeHelper"),
  m_trackingVolumeHelperInstanceName("TrackingVolumeHelper"),
  m_buildBarrel(true),
  m_buildEndcap(true),
  m_buildCsc(true),
  m_buildTgc(true),
  m_identifyActive(false)
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
  declareProperty("BuildBarrelStations",              m_buildBarrel);
  declareProperty("BuildEndcapStations",              m_buildEndcap);
  declareProperty("BuildCSCStations",                 m_buildCsc);
  declareProperty("BuildTGCStations",                 m_buildTgc);
  declareProperty("IdentifyActiveLayers",             m_identifyActive);
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
    
    s = toolSvc()->retrieveTool(m_magFieldToolName, m_magFieldToolInstanceName, m_magFieldTool);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_magFieldToolName << " from ToolSvc. MagneticField will be 0. " << endreq;
    }

    s = toolSvc()->retrieveTool(m_muonStationTypeBuilderName, m_muonStationTypeBuilderInstanceName, m_muonStationTypeBuilder);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_muonStationTypeBuilderName << " from ToolSvc. ";
      log <<" Creation of stations might fail." << endreq;
    }
    // Retrieve the tracking volume helper tool
    s = toolSvc()->retrieveTool(m_trackingVolumeHelperName, m_trackingVolumeHelperInstanceName, m_trackingVolumeHelper);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_trackingVolumeHelperName << " from ToolSvc. ";
      log <<" Creation of Gap Volumes will fail." << endreq;
    }
 
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
    m_muonMagneticField = muonMagneticFieldProperties;

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
    const std::vector<const Trk::TrackingVolume*>* msTypes = buildDetachedTrackingVolumeTypes();
    std::vector<const Trk::TrackingVolume*>::const_iterator msTypeIter = msTypes->begin();

    // position MDT chambers by repeating loop over muon tree
    // link to top tree
    const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
    for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++) 
    {
      const GeoVPhysVol* cv = &(*(top->getChildVol(ichild))); 
      const GeoLogVol* clv = cv->getLogVol();
      std::string vname = clv->getName();
      if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" ) {

        int etaphi = top->getIdOfChildVol(ichild);        // retrive eta/phi indexes
        int sign =( etaphi < 0 ) ? -1 : 1 ;
        etaphi = sign*etaphi;
        int is_mirr = etaphi/1000;
        etaphi = etaphi - is_mirr*1000;
        int eta = etaphi/100;
        int phi = etaphi - eta*100;
        eta = eta*sign;
	MuonGM::MuonStation* gmStation = m_muonMgr->getMuonStation(vname.substr(0,3),eta,phi);
	if ( !gmStation) {
          gmStation = m_muonMgr->getMuonStation(vname.substr(0,4),eta,phi);
        }
        if (!gmStation) log << MSG::ERROR << "muon station not found! "<<vname<<","<<eta<<","<<phi  <<std::endl; 
        std::string stName = (clv->getName()).substr(0,vname.size()-8);
        if (stName.substr(0,1)=="B" && eta < 0 ) {
          stName = (clv->getName()).substr(0,vname.size()-8) + "-";
        }
        if (stName.substr(0,1)=="T" || stName.substr(0,1)=="C") {
          //std::string tgc_name = cv->getChildVol(0)->getLogVol()->getName();
          stName = vname.substr(0,4);
        }
        // loop over prototypes
        const Trk::TrackingVolume* msTV = 0;
        for (msTypeIter = msTypes->begin(); msTypeIter != msTypes->end(); ++msTypeIter) { 
          std::string msTypeName = (*msTypeIter)->volumeName();
          if ( stName == msTypeName.substr(0,stName.size()) ) {
            msTV = *msTypeIter;
	    if (msTV && gmStation) {
	      const Trk::Layer* layerRepresentation = m_muonStationTypeBuilder -> createLayerRepresentation(msTV);
	      const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(stName,msTV,layerRepresentation);
	      HepTransform3D transf = gmStation->getTransform(); 
	      //const Trk::DetachedTrackingVolume* newStat = typeStat->clone(gmStation->getKey(),transf);
	      const Trk::DetachedTrackingVolume* newStat = typeStat->clone(vname,transf);
	      // glue components
	      glueComponents(newStat);
              // eta,phi
              if (msTypeName.substr(0,1)=="C") {
                eta = 1;
		if (transf.getTranslation().z() < 0 ) eta = 0;
		double phic = transf.getTranslation().phi();  
		phi = phic<0 ? int(4*phic/M_PI)+8 : int(4*phic/M_PI);
              } 
	      if (msTypeName.substr(0,1)=="T") {
		bool az = true;
		if (transf.getTranslation().z() < 0 ) az = false;
		double phic = transf.getTranslation().phi();
		if (msTypeName.substr(2,1)=="E" && msTypeName.substr(0,3)!="T4E")
		  phi = phic<0 ? 24*phic/M_PI+48 : 24*phic/M_PI;
		else
		  phi = phic<0 ? 12*phic/M_PI+24 : 12*phic/M_PI;
		if (msTypeName.substr(7,2)=="01") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="02") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="03") eta = az ? 6 : 3;
		if (msTypeName.substr(7,2)=="04") eta = az ? 7 : 2;
		if (msTypeName.substr(7,2)=="05") eta = az ? 8 : 1;
		if (msTypeName.substr(7,2)=="06") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="07") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="08") eta = az ? 6 : 3;
		if (msTypeName.substr(7,2)=="09") eta = az ? 7 : 2;
		if (msTypeName.substr(7,2)=="10") eta = az ? 8 : 1;
		if (msTypeName.substr(7,2)=="11") eta = az ? 9 : 0;
		if (msTypeName.substr(7,2)=="12") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="13") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="14") eta = az ? 6 : 3;
		if (msTypeName.substr(7,2)=="15") eta = az ? 7 : 2;
		if (msTypeName.substr(7,2)=="16") eta = az ? 8 : 1;
		if (msTypeName.substr(7,2)=="17") eta = az ? 9 : 0;
		if (msTypeName.substr(7,2)=="18") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="19") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="20") eta = az ? 5 : 4;
		if (msTypeName.substr(7,2)=="21") eta = az ? 5 : 4;
	      }     
	      // identify layers
	      if (m_identifyActive) identifyLayers(newStat,eta,phi);  
	      mStations.push_back(newStat);
            }
	  }
	}  
	if (!msTV)  log << MSG::INFO  << name() <<" this station has no prototype: " << vname << endreq;    
      }
    }

    /*          
    // position CSC and TGC via loop over prototypes; check !    
    for (msTypeIter = msTypes->begin(); msTypeIter != msTypes->end(); ++msTypeIter) { 
      std::string msTypeName = (*msTypeIter)->volumeName();
      if ( msTypeName.substr(0,1)=="C" ||  msTypeName.substr(0,1)=="T" ) {
        const Trk::TrackingVolume* msTV = *msTypeIter;
        const Trk::Layer* layerRepresentation = m_muonStationTypeBuilder -> createLayerRepresentation(msTV);
        const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(msTypeName,msTV,layerRepresentation);
        // loop over stations to position 
        for (unsigned is=0; is < m_muonMgr->nMuonStation() ;is++) {
          std::string nameStat = ((m_muonMgr->getStation(is)).second)->getName();
          if ( nameStat.substr(0,4) == msTypeName.substr(0,4) ) {
            HepTransform3D transf = ((m_muonMgr->getStation(is)).first)->getTransform();
	    // define station name to indicate eta/phi
            int eta = 1;
            int phi = 0;
            if (msTypeName.substr(0,1)=="C") {
	      if (transf.getTranslation().z() < 0 ) eta = 0;
	      double phic = transf.getTranslation().phi();  
	      phi = phic<0 ? int(4*phic/M_PI)+8 : int(4*phic/M_PI);
            } else {
              bool az = true;
	      if (transf.getTranslation().z() < 0 ) az = false;
              double phic = transf.getTranslation().phi();
              if (msTypeName.substr(2,1)=="E" && msTypeName.substr(0,3)!="T4E")
                phi = phic<0 ? 24*phic/M_PI+48 : 24*phic/M_PI;
              else
                phi = phic<0 ? 12*phic/M_PI+24 : 12*phic/M_PI;
              if (msTypeName.substr(7,2)=="01") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="02") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="03") eta = az ? 6 : 3;
              if (msTypeName.substr(7,2)=="04") eta = az ? 7 : 2;
              if (msTypeName.substr(7,2)=="05") eta = az ? 8 : 1;
              if (msTypeName.substr(7,2)=="06") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="07") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="08") eta = az ? 6 : 3;
              if (msTypeName.substr(7,2)=="09") eta = az ? 7 : 2;
              if (msTypeName.substr(7,2)=="10") eta = az ? 8 : 1;
              if (msTypeName.substr(7,2)=="11") eta = az ? 9 : 0;
              if (msTypeName.substr(7,2)=="12") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="13") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="14") eta = az ? 6 : 3;
              if (msTypeName.substr(7,2)=="15") eta = az ? 7 : 2;
              if (msTypeName.substr(7,2)=="16") eta = az ? 8 : 1;
              if (msTypeName.substr(7,2)=="17") eta = az ? 9 : 0;
              if (msTypeName.substr(7,2)=="18") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="19") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="20") eta = az ? 5 : 4;
              if (msTypeName.substr(7,2)=="21") eta = az ? 5 : 4;
            }  
            const Trk::DetachedTrackingVolume* newStat = typeStat->clone(nameStat,transf);
            // glue components
            glueComponents(newStat);
            // identify layers
            if (m_identifyActive) identifyLayers(newStat,eta,phi);  
            mStations.push_back(newStat);
        }}
      }
    }
    // end CSC/TGC
    */
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonStations=new std::vector<const Trk::DetachedTrackingVolume*>(mStations);

  log << MSG::INFO  << name() << "returns " << (*muonStations).size() << " stations" << endreq;	 
  return muonStations; 
}

const std::vector<const Trk::TrackingVolume*>* Muon::MuonStationBuilder::buildDetachedTrackingVolumeTypes() const 
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building station types" << endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////
      std::vector<const Trk::TrackingVolume*> stations;

   if (m_muonMgr){

      // link to top tree
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++) 
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild))); 
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
	if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && 
            ( (m_buildBarrel && vname.substr(0,1) =="B")
            ||(m_buildEndcap && vname.substr(0,1) =="E")
            ||(m_buildCsc && vname.substr(0,1) =="C")
	    ||(m_buildTgc && vname.substr(0,1) =="T") ) )
	{
          int etaphi = top->getIdOfChildVol(ichild);        // retrive eta/phi indexes
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
            if (stations[in]!=0 && name == stations[in]->volumeName()) is++;
          }
          // if (is!=0 ) std::cout << "prototype exists" << std::endl;
          // if (is==0 ) std::cout << "station shape type:"<< name<< ","<<clv->getShape()->type()<<std::endl; 
          if (is==0 )          
          {
            log << MSG::INFO <<" new station type " << name << "," << clv->getShape()->type() << endreq;    
            log << MSG::INFO <<" prototype built from eta, phi:" << eta << "," << phi << endreq;    
             
            if (name.substr(0,2)=="CS" || name.substr(0,1)=="T") {
              if (m_muonStationTypeBuilder) {
                if (name.substr(0,2)=="CS") { 
                  const Trk::TrackingVolume* csc_station = m_muonStationTypeBuilder->processCscStation(cv, name);   
                  stations.push_back(csc_station); 
                } else {
                  std::vector<const Trk::TrackingVolume*> tgc_stations = m_muonStationTypeBuilder->processTgcStation(cv);   
                  for (unsigned int i=0;i<tgc_stations.size();i++) stations.push_back(tgc_stations[i]); 
                }
              }
            } else {    
	      const GeoTrd* trd=dynamic_cast<const GeoTrd*> (clv->getShape());
	      if (clv->getShape()->type()=="Shift") {
		const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (clv->getShape());
		// std::cout << shift->getOp()->type()<<std::endl; 
		trd = dynamic_cast<const GeoTrd*> (shift->getOp());
	      } 
	      if (clv->getShape()->type()=="Subtraction") {
	       const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (clv->getShape());
               while (  sub->getOpA()->type() =="Subtraction" ) {
                  sub = dynamic_cast<const GeoShapeSubtraction*> (sub->getOpA()); 
               }       
	       trd = dynamic_cast<const GeoTrd*> (sub->getOpA());
	      } 
	      double halfX1=0.;
	      double halfX2=0.;
	      double halfY1=0.;
	      double halfY2=0.;
	      double halfZ=0.;   
	      if (trd) {
		//
              /*
	      std::cout << "dimensions:"<< trd->getXHalfLength1() <<","
                                        << trd->getXHalfLength2() <<","  
                                        << trd->getYHalfLength1() <<","  
                                        << trd->getYHalfLength2() <<","  
			                << trd->getZHalfLength() <<std::endl; 
	      */
	      //
		halfX1 = trd->getXHalfLength1();
		halfX2 = trd->getXHalfLength2();
		halfY1 = trd->getYHalfLength1();
		halfY2 = trd->getYHalfLength2();
		halfZ  = trd->getZHalfLength();              
	      } else {
		std::cout << "station other than Trd shape ?" << std::endl;
	      }
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
		//std::cout << "Trd station dimensions:"<<halfX1<<","<<halfX2<<","<<halfY1<<","<<halfY2<<","<<halfZ<<std::endl;
		Trk::TrapezoidVolumeBounds* envBounds = 0;
		HepTransform3D* transf =new HepTransform3D(); 
		if (halfY1==halfY2) {
		  envBounds = new Trk::TrapezoidVolumeBounds(halfX1,halfX2,halfY1,halfZ);
		  std::cout << "CAUTION!!!: this trapezoid volume does not require XY -> YZ switch" << std::endl;
		}
		if (halfY1!=halfY2 && halfX1 == halfX2 ) {
		  //transf = new HepTransform3D( HepRotateY3D(-90*deg)*HepRotateZ3D(90*deg) );
		  transf = new HepTransform3D( HepRotateY3D(90*deg)* HepRotateZ3D(90*deg) );
		  //envBounds = new Trk::TrapezoidVolumeBounds(halfY1,halfY2,halfX1,halfZ); 
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
	      
	      stations.push_back(newType); 
	    }
	  } // end new station type 
	} // end if "Shift" (station)
      }      
      log << MSG::INFO  << name() << stations.size() <<" station prototypes built " << endreq;    
   }
   
///////////////////////////////////////////////////////////////////////////////////////
   const std::vector<const Trk::TrackingVolume*>* mStations = new std::vector<const Trk::TrackingVolume*>(stations); 
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
  log << MSG::INFO  << name() <<" identifying layers " << endreq;    

  std::string stationName = station->trackingVolume()->volumeName();
  log << MSG::INFO  << " in station " << station->name() << endreq;    

  /*
  if (stationName.substr(0,1)=="C") { 
    int st = stationName.substr(0,3)=="CSS" ? 0 : 1;
    const MuonGM::CscReadoutElement* cscRE = m_muonMgr->getCscReadoutElement(st,eta,phi,0);    
    if (cscRE) {
      for (int gasgap = 0; gasgap < cscRE->Ngasgaps(); gasgap++) {
        int etaId = eta;
        if (etaId < 1 ) etaId = -1;
	Identifier idi = m_cscIdHelper->channelID(st-50,etaId,phi+1,1,gasgap+1,0,cscRE->NetaStrips(gasgap));          
        const HepPoint3D gpi = cscRE->stripPos(idi);
        //const HepPoint3D gp = cscRE->stripPos(eta,0,gasgap+1,0,cscRE->NetaStrips(gasgap));
        const Trk::TrackingVolume* assocVol = station->trackingVolume()->associatedSubVolume(gpi);
        const Trk::Layer* assocLay = 0;
        if (assocVol) assocLay = assocVol->associatedLayer(gpi);
        unsigned int iD = idi;
        if (assocVol && assocLay) assocLay->setLayerType(iD); 
      }
    } 
  }
  */

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
  
    const MuonGM::TgcReadoutElement* tgc = m_muonMgr->getTgcReadoutElement(st,eta,phi);
    if (!tgc) {
      unsigned int phit=0;
      while ( phit<48 ) {
	const MuonGM::TgcReadoutElement* tgct = m_muonMgr->getTgcReadoutElement(st,eta,phit);
        if (tgct && station->trackingVolume()->inside(tgct->center(),0.)) {
          tgc = tgct;
          phi = phit;
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
      if (!assocVol) std::cout << "wrong tgcROE?" << stationName <<"," << eta <<"," << phi << std::endl;
      if (assocVol && assocVol->confinedLayers()) {
	const std::vector<const Trk::Layer*> layers = assocVol->confinedLayers()->arrayObjects();           
	for (unsigned int il=0;il<layers.size();il++) {
	  Identifier wireId  = m_tgcIdHelper->channelID(stationName.substr(0,3),etaSt,phi,il+1,0,1);
	  unsigned int id = wireId;
          layers[il]->setLayerType(id); 
	  // turn wire position into integer to truncate
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
      log << MSG::ERROR << name() << "tgcROE not found for :" << stationName <<","<<eta<<","<<phi<<endreq;         
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
	      if (assocLay) assocLay->setLayerType(iD); 
	      if (assocLay) {
		const Trk::LocalPosition* locPos = (assocLay->surfaceRepresentation()).globalToLocal(gp,0.001);
		if (locPos) assocLay->setRef((*locPos)[Trk::locY]);
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

  /*
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
  */
}
