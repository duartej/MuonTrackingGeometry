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
#include "TrkMagFieldTools/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/TrackingGeometry.h"
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
  m_buildTgc(true)
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
  declareProperty("BuildBarrelStations",              m_buildBarrel);
  declareProperty("BuildEndcapStations",              m_buildEndcap);
  declareProperty("BuildCSCStations",                 m_buildCsc);
  declareProperty("BuildTGCStations",                 m_buildTgc);
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
      if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && 
           ( (m_buildBarrel && vname.substr(0,1) =="B")
           ||(m_buildEndcap && vname.substr(0,1) =="E") ) ) {

        int etaphi = top->getIdOfChildVol(ichild);        // retrive eta/phi indexes
        int sign =( etaphi < 0 ) ? -1 : 1 ;
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
        std::string stName = (clv->getName()).substr(0,4);
        if (stName.substr(0,1)=="T") {
          std::string tgc_name = cv->getChildVol(0)->getLogVol()->getName();
          stName = tgc_name;
        }
        // loop over prototypes
        const Trk::TrackingVolume* msTV = 0;
        for (msTypeIter = msTypes->begin(); msTypeIter != msTypes->end(); ++msTypeIter) { 
          std::string msTypeName = (*msTypeIter)->volumeName();
          if ( stName == msTypeName ) msTV = *msTypeIter;
          /*
          if ( stName.substr(0,2)!="BM" && stName.substr(0,2)!="BI" && stName == msTypeName ) msTV = *msTypeIter;
          if ( stName.substr(0,2)=="BI" && gmStation->getKey()==msTypeName ) msTV = *msTypeIter;
          if ( stName.substr(0,2)=="BM" && gmStation->getKey()==msTypeName ) msTV = *msTypeIter;
          */
        }
        if (msTV && gmStation) {
          const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(stName,msTV);
          HepTransform3D transf = gmStation->getTransform(); 
          const Trk::DetachedTrackingVolume* newStat = typeStat->clone(gmStation->getKey(),transf);
          // glue components
          glueComponents(newStat);
          // identify layers
          identifyLayers(newStat,eta,phi);
	  // std::cout << "new station built at:" << transf.getTranslation() << std::endl;  
          mStations.push_back(newStat);  
        } else {
          log << MSG::INFO  << name() <<" this station has no prototype: " << vname << endreq;    
        }
      }
    }
          
    // position CSC and TGC via loop over prototypes; check !    
    for (msTypeIter = msTypes->begin(); msTypeIter != msTypes->end(); ++msTypeIter) { 
      std::string msTypeName = (*msTypeIter)->volumeName();
      if ( msTypeName.substr(0,1)=="C" ||  msTypeName.substr(0,1)=="T" ) {
        const Trk::TrackingVolume* msTV = *msTypeIter;
        const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(msTypeName,msTV);
        // loop over stations to position 
        for (unsigned is=0; is < m_muonMgr->nMuonStation() ;is++) {
          std::string nameStat = ((m_muonMgr->getStation(is)).second)->getName();
          if ( nameStat.substr(0,4) == msTypeName.substr(0,4) ) {
            HepTransform3D transf = ((m_muonMgr->getStation(is)).first)->getTransform();
	    // define station name to indicate eta/phi
	    //std::cout << "current name of the station:" << nameStat <<"," << msTypeName << std::endl;  
            const Trk::DetachedTrackingVolume* newStat = typeStat->clone(nameStat,transf);
            // glue components
            glueComponents(newStat);
	    // std::cout << "new station built at:" << transf.getTranslation() << std::endl;  
            mStations.push_back(newStat);
        }}
      }
    }
    // end CSC/TGC

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
          // if (!gmStation) std::cout <<" gmStation not found ?" << vname << std::endl; 

          std::string name = (clv->getName()).substr(0,4);
	  // std::cout << clv->getName() << std::endl;
          // is this station known ?        
          // if TGC station, look for 1 component instead
          if (name.substr(0,1)=="T") {
            std::string tgc_name = cv->getChildVol(0)->getLogVol()->getName();
            name = tgc_name;
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
//
void Muon::MuonStationBuilder::identifyLayers(const Trk::DetachedTrackingVolume* station, int eta, int phi ) const
{
  MsgStream log(msgSvc(), name());
  //log << MSG::INFO  << name() <<" identifying layers " << endreq;    

  std::string stationName = station->trackingVolume()->volumeName();
  //log << MSG::INFO  << " in station " << station->name() << endreq;    

  if (stationName.substr(0,1)=="B" || stationName.substr(0,1)=="E" ) { 
    // MDT
    for (int multi = 0; multi < 2; multi++ ) {
      int nameIndex = m_mdtIdHelper->stationNameIndex( stationName.substr(0,3) ); 
      int nameIndexC = nameIndex;
      if (stationName.substr(0,3)=="EIS") nameIndexC = 22; 
      if (stationName.substr(0,3)=="BIM") nameIndexC = 23; 
      const MuonGM::MdtReadoutElement* multilayer = m_muonMgr->getMdtReadoutElement(nameIndexC,eta+8,phi-1,multi);
      if (multilayer) {
	// std::cout << "multilayer found " << multi << "," << multilayer->center() << std::endl;
	if ( !station->trackingVolume()->inside(multilayer->center()) ) std::cout << "multilayer "<< multi << 
                       "outside current volume" << std::endl;
        int nLayers = multilayer->getNLayers();
        int nTubes  = multilayer->getNtubesperlayer();
        for (int layer =1; layer <= nLayers ; layer++) {
          Identifier id = m_mdtIdHelper->channelID(nameIndex,eta,phi,multi+1,layer,nTubes);           
          if (id==0) {
	    //std::cout << "identifier retrieval failed for :"<<stationName<<","<<eta<<","<<phi
	    //		     <<","<<multi<<","<<layer<<std::endl;
          } else {
            // retrieve associated layer
	    Trk::GlobalPosition gp = multilayer->tubePos(id);
            if (!station->trackingVolume()->inside(gp))  std::cout << "associated layer outside the volume: " << layer <<"," << gp << ", identifier of last tube: " <<id <<std::endl;
            const Trk::TrackingVolume* assocVol = station->trackingVolume()->associatedSubVolume(gp);
            const Trk::Layer* assocLay = 0;
            if (assocVol) assocLay = assocVol->associatedLayer(gp);
            unsigned int iD = id;
            if (assocVol && assocLay) assocLay->setLayerType(iD); 
          }    
        }
      } else {
        // std::cout << "multilayer:" << multi << " not found in chamber " << stationName << std::endl;
      }     
    }
    // RPC ?
    bool hasRpc = false;
    const Trk::BinnedArray< Trk::TrackingVolume >* confinedVolumes = station->trackingVolume()->confinedVolumes();
    if (confinedVolumes){
      const std::vector<const Trk::TrackingVolume*>& vols = confinedVolumes->arrayObjects();
      for (unsigned int i=0;i<vols.size();i++) if (vols[i]->volumeName() == "RPC") hasRpc=true;
    }
    if ( hasRpc ) {
      int nameIndex = m_rpcIdHelper->stationNameIndex( stationName.substr(0,3) ); 
      for (int doubletR = 0; doubletR < 2; doubletR++ ) {
      for (int doubletZ = 0; doubletZ < 4; doubletZ++ ) {
        //Identifier iRpcMod = m_rpcIdHelper->elementID(nameIndex,eta+8,phi-1,doubletR);
        //std::cout << "RPC identifier? " << iRpcMod << std::endl;
        //if (iRpcMod==0) {
        //	 std::cout << "RPC identifier retrieval failed for :"<<stationName<<","<<eta<<","<<phi
        //		   <<","<<doubletR << std::endl;
        //} else {
        const MuonGM::RpcReadoutElement* rpc = m_muonMgr->getRpcReadoutElement(nameIndex-2,eta+8,phi-1,doubletR,doubletZ);
        if (rpc) {
	  //std::cout << stationName <<","<<doubletR <<"," << doubletZ << " RPC is inside the current volume? " << 
          //               station->trackingVolume()->inside(rpc->center()) << std::endl;
          if (rpc && station->trackingVolume()->inside(rpc->center())) {
              // retrieve associated layer
	      Trk::GlobalPosition gp = rpc->center();
              const Trk::TrackingVolume* assocVol = station->trackingVolume()->associatedSubVolume(gp);
              const Trk::Layer* assocLay = 0;
              if (assocVol) assocLay = assocVol->associatedLayer(gp);
              Identifier idRpc = rpc->identify();
              unsigned int iD = idRpc;
              if (assocVol && assocLay) assocLay->setLayerType(iD);
          }
        }     
      }}
    }
    //
  } 

  // by now, all the layers should be identified - verify
  if (station->trackingVolume()->confinedVolumes()) {
    const std::vector<const Trk::TrackingVolume*> cVols = station->trackingVolume()->confinedVolumes()->arrayObjects();
    for (unsigned int i=0; i<cVols.size() ; i++) {
      if (cVols[i]->confinedLayers()) {
        const std::vector<const Trk::Layer*> cLays = cVols[i]->confinedLayers()->arrayObjects();
        for (unsigned int il=0; il<cLays.size() ; il++) {
          Identifier id(cLays[il]->layerType());
          if (id==1) std::cout << station->name() << "," << cVols[i]->volumeName() << ", unidentified active layer:" << il << std::endl;
        }
      }
      if (cVols[i]->confinedArbitraryLayers()) {
        const std::vector<const Trk::Layer*>* cLays = cVols[i]->confinedArbitraryLayers();
        for (unsigned int il=0; il<cLays->size() ; il++) {
          Identifier id((*cLays)[il]->layerType());
          if (id==1) std::cout << station->name() << "," << cVols[i]->volumeName() << ", unidentified active layer:" << il << std::endl;
        }
      }
    }
  } 
  // end identification check
 
}