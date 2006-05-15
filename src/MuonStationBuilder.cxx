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
#include "MuonIdHelpers/MuonIdHelper.h"
// Trk
//#include "TrkDetDescrInterfaces/ITrackingVolumeBuilder.h"
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkDetDescrUtils/BinUtility1DR.h"
#include "TrkDetDescrUtils/BinUtility1DZ.h"
#include "TrkDetDescrUtils/BinnedArray.h"
//#include "TrkDetDescrUtils/BinningType.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
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
// StoreGate
#include "StoreGate/StoreGateSvc.h"

// BField
#include "BFieldAth/MagFieldAthena.h"
// temporary
#include "TrkParameters/Perigee.h"

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

// STD
#include <map>

// Gaudi
#include "GaudiKernel/MsgStream.h"

#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoShapeShift.h"
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
  m_muonStationTypeBuilderInstanceName("MuonStationTypeBuilder")
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
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
    // s = m_detStore->retrieve(m_muonMgr, m_muonMgrLocation);

    ds = m_detStore->retrieve(m_muonMgr);

    if (ds.isFailure()) {
        log << MSG::ERROR << "Could not get MuonDetectorManager, no layers for muons will be built. " << endreq;
    }
    
    const MdtIdHelper* mdtHelp = m_muonMgr-> mdtIdHelper(); 
    const RpcIdHelper* rpcHelp = m_muonMgr-> rpcIdHelper();

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

    /*
    s = toolSvc()->retrieveTool(m_trackingVolumeArrayCreatorName, m_trackingVolumeArrayCreatorInstanceName, m_trackingVolumeArrayCreator);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_trackingVolumeArrayCreatorName << " from ToolSvc. ";
      log <<" Creation of LayerArrays might fail." << endreq;
    }   
    */     
 
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
  std::vector<const Trk::DetachedTrackingVolume*> mStations;

  if (m_muonMgr) { 
    // retrieve muon station prototypes from GeoModel
    const std::vector<const Trk::TrackingVolume*>* msTypes = buildDetachedTrackingVolumeTypes();
    std::vector<const Trk::TrackingVolume*>::const_iterator msTypeIter = msTypes->begin();
    for (; msTypeIter != msTypes->end(); ++msTypeIter) { 
      std::string msTypeName = (*msTypeIter)->volumeName();
      const Trk::TrackingVolume* msTV = *msTypeIter;
      const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(msTypeName,msTV);
      // loop over list of muon stations to position them 
      for (unsigned is=0; is < m_muonMgr->nMuonStation() ;is++) {
        std::string nameStat = ((m_muonMgr->getStation(is)).second)->getName();
        if ( nameStat.substr(0,4) == msTypeName ) {
          HepTransform3D transf = ((m_muonMgr->getStation(is)).first)->getTransform();
	  std::cout << "new station to be built at:" << transf.getTranslation() << std::endl;  
          const Trk::DetachedTrackingVolume* newStat = typeStat->clone(nameStat,transf);
	  std::cout << "new station built at:" << transf.getTranslation() << std::endl;  
          mStations.push_back(newStat);
      }}
    } // msType   
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonStations=new std::vector<const Trk::DetachedTrackingVolume*>(mStations);
  std::cout << "muonStationBuilder returns:" << (*muonStations).size() << std::endl;
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
      std::cout << "number of top tree:" << m_muonMgr->getNumTreeTops() << std::endl;
      GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++) 
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild))); 
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
        // if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" &&  vname.substr(0,1) =="B" )     // muon stations
        if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && vname.substr(0,1) !="C" && vname.substr(0,1) !="T" )     // muon stations
	// if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && vname.substr(0,1) =="C")
	//	 if (vname.size()>7 && vname.substr(vname.size()-7,7) =="Station" && vname.substr(0,1) !="T")
	{
          int etaphi = top->getIdOfChildVol(ichild);        // retrive eta/phi indexes
          int sign  = etaphi/fabs(etaphi);
          etaphi = sign*etaphi;
          int is_mirr = etaphi/1000;
          etaphi = etaphi - is_mirr*1000;
          int eta = etaphi/100;
          int phi = etaphi - eta*100;
          eta = eta*sign;
	  std::cout << vname.substr(0,3) <<","<<is_mirr<<","<<eta<<","<<phi<<std::endl;
	  MuonGM::MuonStation* gmStation = m_muonMgr->getMuonStation(vname.substr(0,3),eta,phi);
	  if (gmStation) {
            std::cout << " GM muon station found:" << gmStation->getKey() << std::endl;            
          } else {
            gmStation = m_muonMgr->getMuonStation(vname.substr(0,4),eta,phi);
            if (gmStation) std::cout << " GM muon station found:" << gmStation->getKey() << std::endl;                      } 
          if (!gmStation) std::cout <<" gmStation not found ?" << std::endl; 

          std::string name = (clv->getName()).substr(0,4);
	  std::cout << clv->getName() << std::endl;
          // is this station known ?        
          unsigned is=0; 
          for (unsigned in=0; in< stations.size(); in++) 
          {
            if (name == stations[in]->volumeName()) is++;
          }
          if (is!=0 ) std::cout << "prototype exists" << std::endl;
          if (is==0 ) std::cout << "station shape type:"<< name<< ","<<clv->getShape()->type()<<std::endl; 
          if (is==0 )          
          {
            std::cout <<"new station type:" << name <<","<<clv->getShape()->type()<<std::endl;
            if (name.substr(0,2)=="CS") {
              if (m_muonStationTypeBuilder) {
                const Trk::TrackingVolume* csc_station = m_muonStationTypeBuilder->processCscStation(cv, name); 
                stations.push_back(csc_station); 
              }
            } else {    
            const GeoTrd* trd=dynamic_cast<const GeoTrd*> (clv->getShape());
            if (clv->getShape()->type()=="Shift") {
	       const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (clv->getShape());
               std::cout << shift->getOp()->type()<<std::endl; 
	       trd = dynamic_cast<const GeoTrd*> (shift->getOp());
	    } 
            double halfX1;
            double halfX2;
            double halfY1;
            double halfY2;
            double halfZ;   
            if (trd ) 
            {
	      std::cout << "dimensions:"<< trd->getXHalfLength1() <<","
                                        << trd->getXHalfLength2() <<","  
                                        << trd->getYHalfLength1() <<","  
                                        << trd->getYHalfLength2() <<","  
			                << trd->getZHalfLength() <<std::endl; 
             
              halfX1 = trd->getXHalfLength1();
              halfX2 = trd->getXHalfLength2();
              halfY1 = trd->getYHalfLength1();
              halfY2 = trd->getYHalfLength2();
              halfZ  = trd->getZHalfLength();              
            } else {
	      std::cout << "station other than Trd shape ?"<< std::endl;
            }
            // define enveloping volume
            Trk::TrackingVolumeArray* confinedVolumes = 0; 
	    Trk::Volume* envelope;
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
	      std::cout << "Trd station dimensions:"<<halfX1<<","<<halfX2<<","<<halfY1<<","<<halfY2<<","<<halfZ<<std::endl;
	      Trk::TrapezoidVolumeBounds* envBounds;
              HepTransform3D* transf =new HepTransform3D(); 
              if (halfY1==halfY2) {
                envBounds = new Trk::TrapezoidVolumeBounds(halfX1,halfX2,halfY1,halfZ);
		std::cout << "CAUTION!!!: this trapezoid volume does not require XY -> YZ switch" << std::endl;
              }
              if (halfY1!=halfY2 && halfX1 == halfX2 ) {
	        transf = new HepTransform3D( HepRotateZ3D(90*deg) );
	        envBounds = new Trk::TrapezoidVolumeBounds(halfY1,halfY2,halfX1,halfZ); 
              }
              if (halfX1!=halfX2 && halfY1!=halfY2 ) std::cout << "station envelope arbitrary trapezoid?" << std::endl;
              // station components
              if (m_muonStationTypeBuilder) confinedVolumes = 
	          m_muonStationTypeBuilder->processTrdStationComponents(cv,envBounds); 
              // enveloping volume
	      envelope= new Trk::Volume(transf,envBounds);
            }
	    unsigned int ngc = cv->getNChildVols();
	    std::vector<Trk::TrackingVolume*> components;
            for (unsigned ic =0; ic < ngc; ic++){
              const GeoVPhysVol* gcv = &(*(cv->getChildVol(ic))); 
              const GeoLogVol* gclv = gcv->getLogVol();
	      std::string gcname=gclv->getName();
              HepTransform3D transf = cv->getXToChildVol(ic);
            } // station components

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
      std::cout << stations.size() << "station types read" << std::endl;
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

