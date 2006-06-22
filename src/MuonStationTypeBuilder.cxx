///////////////////////////////////////////////////////////////////
// MuonStationTypeBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonStationTypeBuilder.h"
//MuonSpectrometer include
#include "MuonGeoModel/MuonDetectorManager.h"
//#include "MuonGeoModel/StationSelector.h"
//#include "MuonGeoModel/Station.h"
#include "MuonGeoModel/MuonStation.h"
#include "MuonGeoModel/MdtReadoutElement.h"
#include "MuonIdHelpers/MuonIdHelper.h"
// Trk
//#include "TrkDetDescrInterfaces/ITrackingVolumeBuilder.h"
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkDetDescrUtils/BinUtility1DX.h"
#include "TrkDetDescrUtils/BinUtility1DY.h"
#include "TrkDetDescrUtils/BinnedArray.h"
#include "TrkDetDescrUtils/NavBinnedArray1D.h"
//#include "TrkDetDescrUtils/BinningType.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/DoubleTrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkSurfaces/TrapezoidBounds.h"
#include "TrkSurfaces/DiamondBounds.h"
#include "TrkGeometry/CylinderLayer.h"
#include "TrkGeometry/DiscLayer.h"
#include "TrkGeometry/PlaneLayer.h"
#include "TrkGeometry/MaterialProperties.h"
#include "TrkGeometry/LayerMaterialProperties.h"
#include "TrkGeometry/HomogenousLayerMaterial.h"
#include "TrkGeometry/OverlapDescriptor.h"
//#include "TrkGeometry/SimplifiedMaterialProperties.h"
#include "TrkMagFieldTools/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/DetachedTrackingVolume.h"
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
#include "CLHEP/Geometry/Transform3D.h"

// STD
#include <map>

// Gaudi
#include "GaudiKernel/MsgStream.h"

#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoShapeShift.h"
#include "GeoModelKernel/GeoShapeUnion.h"
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoTube.h"

static const InterfaceID IID_IMuonStationTypeBuilder("MuonStationTypeBuilder", 1, 0);

const InterfaceID& Muon::MuonStationTypeBuilder::interfaceID()
{
  return IID_IMuonStationTypeBuilder;
}


// constructor
Muon::MuonStationTypeBuilder::MuonStationTypeBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgTool(t,n,p),
  m_muonMgrLocation("MuonMgr"),
  m_trackingVolumeArrayCreator(0),
  m_trackingVolumeArrayCreatorName("Trk::TrackingVolumeArrayCreator"),
  m_trackingVolumeArrayCreatorInstanceName("TrackingVolumeArrayCreator"),
  m_magFieldTool(0),
  m_magFieldToolName("Trk::MagneticFieldTool"),
  m_magFieldToolInstanceName("ATLAS_TrackingMagFieldTool"),
  m_mdtTubeMat(0),
  m_mdtFoamMat(0),
  m_rpc46(0),
  m_rpcDed50(0),
  m_matCSC01(0),
  m_matCSCspacer1(0),
  m_matCSC02(0),
  m_matCSCspacer2(0),
  m_matTGC01(0),
  m_matTGC06(0)
{
  declareInterface<Muon::MuonStationTypeBuilder>(this);

  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
}

// destructor
Muon::MuonStationTypeBuilder::~MuonStationTypeBuilder()
{}

// Athena standard methods
// initialize
StatusCode Muon::MuonStationTypeBuilder::initialize()
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

   m_materialConverter= new Trk::GeoMaterialConverter();
     
    log << MSG::INFO  << name() <<" initialize() successful" << endreq;    
    
  return StatusCode::SUCCESS;
}


const Trk::TrackingVolumeArray* Muon::MuonStationTypeBuilder::processBoxStationComponents(const GeoVPhysVol* mv, Trk::CuboidVolumeBounds* envelope) const
{

    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" processing station components for " <<mv->getLogVol()->getName()<< endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////
   
   double tolerance = 0.0001;   

   // loop over children volumes; ( make sure they do not exceed enveloping volume boundaries ?)
   // split into connected subvolumes ( assume ordering along X unless otherwise )
      std::vector<Trk::Volume*> compVol;
      std::vector<std::string> compName;
      std::vector<const GeoVPhysVol*> compGeo;
      std::vector<HepTransform3D*> compTransf;
      double halfX;
      double halfY;
      double halfZ;   
      double halfX1;
      double halfX2;
      double halfY1;
      double halfY2;
      for (unsigned int ich =0; ich< mv->getNChildVols(); ++ich) 
      {
	// std::cout << "next component:"<< ich << std::endl;
        const GeoVPhysVol* cv = &(*(mv->getChildVol(ich))); 
        const GeoLogVol* clv = cv->getLogVol();
        HepTransform3D transf = mv->getXToChildVol(ich);        
      // TEMPORARY CORRECTION 
        if ( (mv->getLogVol()->getName()).substr(0,3)=="BMF" && (clv->getName()).substr(0,2)=="LB" ) {
            std::cout << "TEMPORARY MANUAL CORRECTION OF BMF SPACER LONG BEAM POSITION" << std::endl;
            transf = transf * HepTranslate3D(-37.5,0.,0.);
        } 
        // std::cout << "component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<"," <<transf.getTranslation()<<std::endl;
        // retrieve volumes for components
	Trk::VolumeBounds* volBounds; 
	Trk::Volume* vol; 
        if ( clv->getShape()->type()=="Trd")
	{
          const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
          halfX1 = trd->getXHalfLength1();
          halfX2 = trd->getXHalfLength2();
          halfY1 = trd->getYHalfLength1();
          halfY2 = trd->getYHalfLength2();
          halfZ  = trd->getZHalfLength();
          volBounds = new Trk::CuboidVolumeBounds(fmax(halfX1,halfX2),fmax(halfY1,halfY2),halfZ);
        } 
        if ( clv->getShape()->type()=="Box")
        {
	  const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
	  halfX1 = box->getXHalfLength();
          halfX2 = halfX1;
          halfY1 = box->getYHalfLength();  
          halfY2 = halfY1;
	  halfZ  = box->getZHalfLength(); 
          volBounds = new Trk::CuboidVolumeBounds(halfX1,halfY1,halfZ);
        }
	// std::cout << "dimensions:"<<halfX1<<","<<halfX2<<","<<halfY1<<","<<halfY2<<","<<halfZ<<std::endl;
        if ( clv->getShape()->type()!="Trd" && clv->getShape()->type()!="Box" ) {
 	  std::cout << "WARNING:component shape not Box nor Trapezoid, determining the x size from subcomponents" << std::endl; 
          double xSize = get_x_size(cv);
          volBounds = new Trk::CuboidVolumeBounds(xSize,envelope->halflengthY(),envelope->halflengthZ());
        }
	vol = new Trk::Volume(new HepTransform3D(transf),volBounds);
	//std::cout <<"volume center:"<< vol->center() << ","<< ich << std::endl;
	std::string cname = clv->getName();
        if (cname.substr(0,4) == (mv->getLogVol()->getName()).substr(0,4)) cname = cname.substr(4,cname.size()-4);  
        // order in X
        if (compVol.size()==0 || vol->center()[0]>=compVol.back()->center()[0]){
          compVol.push_back(vol);
          compName.push_back(cname);
          compGeo.push_back(cv);
          compTransf.push_back(new HepTransform3D(transf));
        } else {
	  std::vector<Trk::Volume*>::iterator volIter0=compVol.begin();
	  std::vector<Trk::Volume*>::iterator volIter=compVol.begin();
	  std::vector<std::string>::iterator  nameIter=compName.begin();
	  std::vector<const GeoVPhysVol*>::iterator  geoIter=compGeo.begin();
	  std::vector<HepTransform3D*>::iterator  transfIter=compTransf.begin();
          while ( vol->center()[0]>= (*volIter)->center()[0]) {volIter++;nameIter++;geoIter++;transfIter++;}
          compVol.insert(volIter,vol);
          compName.insert(nameIter,cname);
          compGeo.insert(geoIter,cv);
          compTransf.insert(transfIter,new HepTransform3D(transf));
        } 
      } // loop over components
      /* 
      // check components ordering
      for (unsigned i=0; i<compVol.size();i++){
         std::cout << compVol[i]->center()[0]<<" "<<compName[i]<<std::endl; 
      } 
      */
      // define enveloping volumes for each "technology"
      std::vector<const Trk::TrackingVolume*> trkVols;
      double envX = envelope->halflengthX();
      double envY = envelope->halflengthY();
      double envZ = envelope->halflengthZ();
      double currX = -envX;
      double maxX = envX;
      bool openSpacer = false;
      bool openRpc = false;
      std::vector<const GeoVPhysVol*> geoSpacer;
      std::vector<const GeoVPhysVol*> geoRpc;
      std::vector<HepTransform3D*> transfSpacer;
      std::vector<HepTransform3D*> transfRpc;
      double spacerlowXsize=0; 
      double spaceruppXsize=0; 
      double rpclowXsize=0; 
      double rpcuppXsize=0; 
      std::vector<double>* volSteps = new std::vector<double>;
      for (unsigned i=0; i<compVol.size();i++){
        bool comp_processed = false;
        const Trk::CuboidVolumeBounds* compBounds = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(compVol[i]->volumeBounds()));
        double lowX = compVol[i]->center()[0]-compBounds->halflengthX();
        double uppX = compVol[i]->center()[0]+compBounds->halflengthX();
	if ( lowX < currX ) std::cout<<"Warning: we have a clash between components here!"<< std::endl;
	if ( uppX > maxX ) std::cout<<"Warning: we have a clash between component and envelope!"<< std::endl;
        // close Rpc if no futher components
        if (openRpc  && compName[i].substr(0,3) != "RPC" && compName[i].substr(0,3) != "Ded"){
          // low edge of current volume
          double Xcurr = compVol[i]->center()[0]-compBounds->halflengthX();
          if (Xcurr>= currX+rpclowXsize+rpcuppXsize) {
            Trk::CuboidVolumeBounds* rpcBounds = new Trk::CuboidVolumeBounds(0.5*(Xcurr-currX),envY,envZ);
	    Trk::Volume* rpcVol =new Trk::Volume(new HepTranslate3D(currX+rpcBounds->halflengthX(),0.,0.),rpcBounds);
  	    //std::cout << "new Rpc volume:position:" << (rpcVol->transform()).getTranslation() << std::endl; 
	    //std::cout << "new Rpc volume:dimensions" <<0.5*(Xcurr-currX)<<","<<envY<<","<<envZ << std::endl; 
            const Trk::TrackingVolume* rpcTrkVol = processRpc(rpcVol ,geoRpc,transfRpc);
            trkVols.push_back(rpcTrkVol);  
            volSteps->push_back(Xcurr-currX); 
            currX = Xcurr;
            openRpc = false;
          } else {
	    std::cout << "clash in Rpc definition!" << std::endl; 
          } 
        }   
        // close spacer if no further components
        if (openSpacer &&  compName[i].substr(0,1) != "C" && compName[i].substr(0,2) != "LB"){
          // low edge of current volume
          double Xcurr = compVol[i]->center()[0]-compBounds->halflengthX();
          if (Xcurr-currX-(spacerlowXsize+spaceruppXsize)>= -tolerance ) {
            Trk::CuboidVolumeBounds* spacerBounds = new Trk::CuboidVolumeBounds(0.5*(Xcurr-currX),envY,envZ);
	    Trk::Volume* spacerVol =new Trk::Volume(new HepTranslate3D(currX+spacerBounds->halflengthX(),0.,0.),spacerBounds);
  	    //std::cout << "new Spacer volume:position:" << (spacerVol->transform()).getTranslation() << std::endl; 
	    //std::cout << "new Spacer volume:dimensions:" <<0.5*(Xcurr-currX)<<","<<envY<<","<<envZ << std::endl; 
            const Trk::TrackingVolume* spacerTrkVol = processSpacer(spacerVol ,geoSpacer, transfSpacer);
            trkVols.push_back(spacerTrkVol);  
            volSteps->push_back(Xcurr-currX); 
            currX = Xcurr;
            openSpacer = false;
          } else {
	    std::cout << "currX,Xcurr,lowX,uppX" << currX <<"," << Xcurr<<","<<spacerlowXsize <<"," << spaceruppXsize<<std::endl;
	    std::cout << Xcurr-currX << "," << spacerlowXsize+spaceruppXsize << std::endl;  
	    std::cout << "clash in spacer definition!" << std::endl; 
          } 
        }   
        if (compName[i].substr(0,3) == "RPC" || compName[i].substr(0,3) == "Ded" ) {
          if (!openRpc) {
            openRpc = true;
            geoRpc.clear();
	    geoRpc.push_back(compGeo[i]);
            transfRpc.clear();
	    transfRpc.push_back(compTransf[i]);
            //establish temporary volume size
            rpclowXsize = compVol[i]->center()[0]-currX;
            rpcuppXsize = compBounds->halflengthX();
            // check clash at low edge
            if (rpclowXsize < compBounds->halflengthX()) std::cout << "WARNING at rpc low edge - not enough space" << std::endl;
          } else {
            geoRpc.push_back(compGeo[i]);
	    transfRpc.push_back(compTransf[i]);
            // check temporary volume size
            if ( compVol[i]->center()[0]-currX < compBounds->halflengthX()) std::cout << "WARNING at rpc low edge - not enough space" << std::endl;
            if ( compVol[i]->center()[0]+compBounds->halflengthX() > currX + rpclowXsize+rpcuppXsize) 
	      rpcuppXsize += (compVol[i]->center()[0]+compBounds->halflengthX())-( currX + rpclowXsize+rpcuppXsize);
          }
          comp_processed = true;
        }
        if (compName[i].substr(0,1) == "C" || compName[i].substr(0,2) == "LB" ) {
          if (!openSpacer) {
            openSpacer = true;
            geoSpacer.clear();
	    geoSpacer.push_back(compGeo[i]);
            transfSpacer.clear();
	    transfSpacer.push_back(compTransf[i]);
            //establish temporary volume size
            spacerlowXsize = compVol[i]->center()[0]-currX;
            spaceruppXsize = compBounds->halflengthX();
            // check clash at low edge
            if (spacerlowXsize < compBounds->halflengthX()) std::cout << "WARNING at spacer low edge - not enough space" << std::endl;
          } else {
            geoSpacer.push_back(compGeo[i]);
	    transfSpacer.push_back(compTransf[i]);
            // check temporary volume size
            if ( compVol[i]->center()[0]-currX < compBounds->halflengthX()) std::cout << "WARNING at spacer low edge - not enough space" << std::endl;
            if ( compVol[i]->center()[0]+compBounds->halflengthX() > currX + spacerlowXsize+spaceruppXsize) 
	      spaceruppXsize += (compVol[i]->center()[0]+compBounds->halflengthX())-( currX + spacerlowXsize+spaceruppXsize);
          }
          comp_processed = true;
        }
        if (compName[i].substr(0,3) == "MDT") {
          Trk::Volume* mdtVol;
	  Trk::CuboidVolumeBounds* mdtBounds; 
          if (lowX == currX) {
	    mdtBounds = new Trk::CuboidVolumeBounds(compBounds->halflengthX(),envY,envZ);
            mdtVol = new Trk::Volume(new HepTransform3D(compVol[i]->transform()),mdtBounds);
	  } else {
	    std::cout << "lowX,currX:" << lowX <<","<<currX << std::endl;
	    std::cout << "increase the basic Mdt volume" << std::endl;
          }
	  //std::cout << "new Mdt volume:position:" << (mdtVol->transform()).getTranslation() << std::endl; 
	  //std::cout << "new Mdt volume:dimensions:" <<compBounds->halflengthX()<<","<<envY<<","<<envZ << std::endl; 
          const Trk::TrackingVolume* mdtTrkVol = processMdtBox(mdtVol,compGeo[i],compTransf[i]);
          trkVols.push_back(mdtTrkVol); 
          volSteps->push_back(2.*mdtBounds->halflengthX()); 
          currX += 2.*mdtBounds->halflengthX();
          comp_processed = true;
        }
        if ( !comp_processed ) std::cout << "unknown technology:" <<compName[i]<<std::endl;
      } // end loop over station children

      // there may be a spacer still open
      if (openSpacer) {
        if (maxX >= currX+spacerlowXsize+spaceruppXsize) {
          Trk::CuboidVolumeBounds* spacerBounds = new Trk::CuboidVolumeBounds(0.5*(maxX-currX),envY,envZ);
	  Trk::Volume* spacerVol =new Trk::Volume(new HepTranslate3D(currX+spacerBounds->halflengthX(),0.,0.),spacerBounds);
	  //std::cout << "new Spacer volume:position:" << (spacerVol->transform()).getTranslation() << std::endl; 
	  //std::cout << "new Spacer volume:dimensions:" <<0.5*(maxX-currX)<<","<<envY<<","<<envZ << std::endl; 
          const Trk::TrackingVolume* spacerTrkVol = processSpacer(spacerVol ,geoSpacer, transfSpacer);
          trkVols.push_back(spacerTrkVol);  
          volSteps->push_back(maxX-currX); 
          currX = maxX;
          openSpacer = false;
        } else {
	  //std::cout <<"currX,maxX,lowX,uppX:"<< currX<<"," << maxX <<","<<spacerlowXsize<<"," <<spaceruppXsize<< std::endl; 
	  //std::cout << "clash in spacer definition!(last volume)" << std::endl; 
        }          
      }
      // there may be an Rpc still open
      if (openRpc) {
	//std::cout << "maxX, other:"<< maxX <<"," << currX+rpclowXsize+rpcuppXsize << std::endl;
        if (maxX >= currX+rpclowXsize+rpcuppXsize) {
          Trk::CuboidVolumeBounds* rpcBounds = new Trk::CuboidVolumeBounds(0.5*(maxX-currX),envY,envZ);
	  Trk::Volume* rpcVol =new Trk::Volume(new HepTranslate3D(currX+rpcBounds->halflengthX(),0.,0.),rpcBounds);
	  //std::cout << "new Rpc volume:position:" << (rpcVol->transform()).getTranslation() << std::cout; 
	  //std::cout << "new Rpc volume:dimensions:" <<0.5*(maxX-currX)<<","<<envY<<","<<envZ << std::cout; 
          const Trk::TrackingVolume* rpcTrkVol = processRpc(rpcVol ,geoRpc,transfRpc);
          trkVols.push_back(rpcTrkVol);  
          volSteps->push_back(maxX-currX); 
          currX = maxX;
          openRpc = false;
        } else {
	  std::cout << "clash in Rpc definition!(last volume)" << std::endl; 
        }          
      }
      // create VolumeArray (1DX) 
      const Trk::TrackingVolumeArray* components = 0;
      const std::vector<const Trk::TrackingVolume*>* compVols = new std::vector<const Trk::TrackingVolume*>( trkVols );
      Trk::BinUtility* binUtility = new Trk::BinUtility1DX( -( envelope->halflengthX() ), volSteps);
      if (m_trackingVolumeArrayCreator)  components = m_trackingVolumeArrayCreator->cuboidVolumesArrayNav( *compVols, binUtility, false);
      // std::cout << "tracking volume array created" << std::endl;
     
   return components;  
}

const Trk::TrackingVolumeArray* Muon::MuonStationTypeBuilder::processTrdStationComponents(const GeoVPhysVol* mv, Trk::TrapezoidVolumeBounds* envelope ) const
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" processing station components" << endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////

    double tolerance = 0.0001;   
   
   // loop over children volumes; ( make sure they do not exceed enveloping volume boundaries ?)
   // split into connected subvolumes ( assume ordering along X unless otherwise )
      std::vector<Trk::Volume*> compVol;
      std::vector<std::string> compName;
      std::vector<const GeoVPhysVol*> compGeo;
      std::vector<HepTransform3D*> compTransf;
      double halfX;
      double halfY;
      double halfZ;   
      double halfX1;
      double halfX2;
      double halfY1;
      double halfY2;
      for (unsigned int ich =0; ich< mv->getNChildVols(); ++ich) 
      {
	// std::cout << "next component:"<< ich << std::endl;
        const GeoVPhysVol* cv = &(*(mv->getChildVol(ich))); 
        const GeoLogVol* clv = cv->getLogVol();
        HepTransform3D transf = mv->getXToChildVol(ich);        
        // std::cout << "component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","  <<transf.getTranslation()<<std::endl;
        // retrieve volumes for components
	Trk::VolumeBounds* volBounds; 
	Trk::Volume* vol; 
        if ( clv->getShape()->type()=="Trd")
	{
          const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
          halfX1 = trd->getXHalfLength1();
          halfX2 = trd->getXHalfLength2();
          halfY1 = trd->getYHalfLength1();
          halfY2 = trd->getYHalfLength2();
          halfZ  = trd->getZHalfLength();
          if (halfX1==halfX2 && halfY1==halfY2) volBounds = new Trk::CuboidVolumeBounds(fmax(halfX1,halfX2),fmax(halfY1,halfY2),halfZ);
          if (halfX1==halfX2 && halfY1!=halfY2 ) {
             transf = transf*HepRotateY3D(90*deg)*HepRotateZ3D(90*deg);
             volBounds = new Trk::TrapezoidVolumeBounds(halfY1,halfY2,halfZ,halfX1);
          }
          if (halfX1!=halfX2 && halfY1==halfY2 ) {
             volBounds = new Trk::TrapezoidVolumeBounds(halfX1,halfX2,halfY1,halfZ);
          }
          if (!volBounds) std::cout << "volume shape for component not recognized" << std::endl;
        } 
        if ( clv->getShape()->type()=="Box")
        {
	  const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
	  halfX1 = box->getXHalfLength();
          halfX2 = halfX1;
          halfY1 = box->getYHalfLength();  
          halfY2 = halfY1;
	  halfZ  = box->getZHalfLength(); 
          volBounds = new Trk::CuboidVolumeBounds(halfX1,halfY1,halfZ);
        }
	// std::cout << "dimensions:"<<halfX1<<","<<halfX2<<","<<halfY1<<","<<halfY2<<","<<halfZ<<std::endl;
        if ( clv->getShape()->type()!="Trd" && clv->getShape()->type()!="Box" ) {
 	  std::cout << "WARNING:component shape not Box nor Trapezoid, determining the x size from subcomponents" << std::endl;                  double xSize = get_x_size(cv);
          // printChildren(cv);
          transf = transf*HepRotateY3D(90*deg)*HepRotateZ3D(90*deg);
          volBounds = new Trk::TrapezoidVolumeBounds(envelope->minHalflengthX(),envelope->maxHalflengthX(),envelope->halflengthY(),xSize);
        }
	vol = new Trk::Volume(new HepTransform3D(transf),volBounds);
	// std::cout <<"volume center:"<< vol->center() << ","<< ich << std::endl;
	std::string cname = clv->getName();
        if (cname.substr(0,4) == (mv->getLogVol()->getName()).substr(0,4)) cname = cname.substr(4,cname.size()-4);  
        // order in X
        if (compVol.size()==0 || vol->center()[0]>=compVol.back()->center()[0]){
          compVol.push_back(vol);
          compName.push_back(cname);
          compGeo.push_back(cv);
          compTransf.push_back(new HepTransform3D(transf));
        } else {
	  std::vector<Trk::Volume*>::iterator volIter0=compVol.begin();
	  std::vector<Trk::Volume*>::iterator volIter=compVol.begin();
	  std::vector<std::string>::iterator  nameIter=compName.begin();
	  std::vector<const GeoVPhysVol*>::iterator  geoIter=compGeo.begin();
	  std::vector<HepTransform3D*>::iterator  transfIter=compTransf.begin();
          while ( vol->center()[0]>= (*volIter)->center()[0]) {volIter++;nameIter++;geoIter++;transfIter++;}
          compVol.insert(volIter,vol);
          compName.insert(nameIter,cname);
          compGeo.insert(geoIter,cv);
          compTransf.insert(transfIter,new HepTransform3D(transf));
        } 
      } // loop over components
      /*
      // check components ordering
      for (unsigned i=0; i<compVol.size();i++){
	std::cout << compVol[i]->center()[0]<<" "<<compName[i]<<std::endl; 
      } 
      */
      // define enveloping volumes for each "technology"
      std::vector<const Trk::TrackingVolume*> trkVols;
      double envX1 = envelope->minHalflengthX();
      double envX2 = envelope->maxHalflengthX();
      double envY = envelope->halflengthY();
      double envZ = envelope->halflengthZ();
      //
      double currX = -envZ;
      double maxX = envZ;
      //
      bool openSpacer = false;
      bool openRpc = false;
      std::vector<const GeoVPhysVol*> geoSpacer;
      std::vector<const GeoVPhysVol*> geoRpc;
      std::vector<HepTransform3D*> transfSpacer;
      std::vector<HepTransform3D*> transfRpc;
      double spacerlowXsize=0; 
      double spaceruppXsize=0; 
      double rpclowXsize=0; 
      double rpcuppXsize=0; 
      double Xcurr=0;
      double lowX;
      double uppX;
      std::vector<double>* volSteps = new std::vector<double>;
      for (unsigned i=0; i<compVol.size();i++){
        bool comp_processed = false;
        const Trk::CuboidVolumeBounds* compCubBounds = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(compVol[i]->volumeBounds()));
        const Trk::TrapezoidVolumeBounds* compTrdBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(compVol[i]->volumeBounds()));
        if (compCubBounds) {
          lowX = compVol[i]->center()[0]-compCubBounds->halflengthX();
          uppX = compVol[i]->center()[0]+compCubBounds->halflengthX();
	  if ( lowX < currX ) std::cout<<"Warning: we have a clash between components here!"<< std::endl;
	  if ( uppX > maxX ) std::cout<<"Warning: we have a clash between component and envelope!"<< std::endl;
          // low edge of current volume
          Xcurr = compVol[i]->center()[0]-compCubBounds->halflengthX();
        } 
        if (compTrdBounds) {
          lowX = compVol[i]->center()[0]-compTrdBounds->halflengthZ();
          uppX = compVol[i]->center()[0]+compTrdBounds->halflengthZ();
	  if ( lowX < currX ) std::cout<<"Warning: we have a clash between components here!"<< std::endl;
	  if ( uppX > maxX ) std::cout<<"Warning: we have a clash between component and envelope!"<< std::endl;
          // low edge of current volume
          Xcurr = compVol[i]->center()[0]-compTrdBounds->halflengthZ();
        }
        if (!compCubBounds && !compTrdBounds ) {
           std::cout << "unknown volume shape" << std::endl; 
           return 0;
        }
        // close Rpc if no futher components
        if (openRpc  && compName[i].substr(0,3) != "RPC" && compName[i].substr(0,3) != "Ded"){
	  std::cout << " RPC components for endcaps not coded " << std::endl;
          if (Xcurr>= currX+rpclowXsize+rpcuppXsize) {
            Trk::TrapezoidVolumeBounds* rpcBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,0.5*(Xcurr-currX));
	    Trk::Volume* rpcVol =new Trk::Volume(new HepTranslate3D(0.,currX+rpcBounds->halflengthZ(),0.),rpcBounds);
  	    //std::cout << "new Rpc volume:position:" << (rpcVol->transform()).getTranslation() << std::endl; 
	    //std::cout << "new Rpc volume:dimensions" <<envX1<<","<<envX2<<","<<0.5*(Xcurr-currX)<<","<<envZ << std::endl; 
            //const Trk::TrackingVolume* rpcTrkVol = processRpc(rpcVol ,geoRpc,transfRpc);
            //trkVols.push_back(rpcTrkVol);  
            volSteps->push_back(Xcurr-currX);
            currX = Xcurr;
            openRpc = false;
          } else {
	    std::cout << "clash in Rpc definition!" << std::endl; 
          } 
        }   
        // close spacer if no further components
        if (openSpacer &&  compName[i].substr(0,1) != "C" && compName[i].substr(0,2) != "LB"){
          if (Xcurr-currX-(spacerlowXsize+spaceruppXsize)>= -tolerance ) {
            Trk::TrapezoidVolumeBounds* spacerBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,0.5*(Xcurr-currX));
            //HepTransform3D tr = HepRotateZ3D(90*deg)*HepRotateX3D(90*deg)*HepTranslate3D(currX+spacerBounds->halflengthZ(),0.,0.);
            HepTransform3D tr = HepTranslate3D(currX+spacerBounds->halflengthZ(),0.,0.)* HepRotateY3D(90*deg)*HepRotateZ3D(90*deg);
	    Trk::Volume* spacerVol =new Trk::Volume(new HepTransform3D(tr),spacerBounds);
            const Trk::TrackingVolume* spacerTrkVol = processSpacer(spacerVol ,geoSpacer, transfSpacer);
            trkVols.push_back(spacerTrkVol);  
            volSteps->push_back(Xcurr-currX);
            currX = Xcurr;
            openSpacer = false;
          } else {
	    std::cout << "currX,Xcurr,lowX,uppX" << currX <<"," << Xcurr<<","<<spacerlowXsize <<"," << spaceruppXsize<<std::endl;
	    std::cout << Xcurr-currX << "," << spacerlowXsize+spaceruppXsize << std::endl;  
	    std::cout << "clash in spacer definition!" << std::endl; 
          } 
        }   
        if (compName[i].substr(0,3) == "RPC" || compName[i].substr(0,3) == "Ded" ) {
	  std::cout << " RPC components for endcaps not coded " << std::endl;
        }
        if (compName[i].substr(0,1) == "C" || compName[i].substr(0,2) == "LB" ) {
          if (!openSpacer) {
	    // std::cout << "opening spacer,"<<currX<<std::endl;
            openSpacer = true;
            geoSpacer.clear();
	    geoSpacer.push_back(compGeo[i]);
            transfSpacer.clear();
	    transfSpacer.push_back(compTransf[i]);
            //establish temporary volume size
            spacerlowXsize = compVol[i]->center()[0]-currX;
            if (compCubBounds) {
              spaceruppXsize = compCubBounds->halflengthX();
              // check clash at low edge
              if (spacerlowXsize < compCubBounds->halflengthX()) 
                std::cout << "WARNING at spacer low edge - not enough space" << std::endl;
            }
            if (compTrdBounds) {
              spaceruppXsize = compTrdBounds->halflengthZ();
              // check clash at low edge
              if (spacerlowXsize < compTrdBounds->halflengthZ()) 
              std::cout << "WARNING at spacer low edge - not enough space" << std::endl;
            } 
          } else {
            geoSpacer.push_back(compGeo[i]);
	    transfSpacer.push_back(compTransf[i]);
            // check temporary volume size
            if (compCubBounds) {
               if ( compVol[i]->center()[0]-currX < compCubBounds->halflengthX()) 
                  std::cout << "WARNING at spacer low edge - not enough space" << std::endl;
               if ( compVol[i]->center()[0]+compCubBounds->halflengthX() > currX + spacerlowXsize+spaceruppXsize) 
	      spaceruppXsize += (compVol[i]->center()[0]+compCubBounds->halflengthX())-( currX + spacerlowXsize+spaceruppXsize);
            }
            if (compTrdBounds) {
               if ( compVol[i]->center()[0]-currX < compTrdBounds->halflengthZ()) 
                  std::cout << "WARNING at spacer low edge - not enough space" << std::endl;
               if ( compVol[i]->center()[0]+compTrdBounds->halflengthZ() > currX + spacerlowXsize+spaceruppXsize) 
	      spaceruppXsize += (compVol[i]->center()[0]+compTrdBounds->halflengthZ())-( currX + spacerlowXsize+spaceruppXsize);
            }
          }
          comp_processed = true;
        }
        if (compName[i].substr(0,3) == "MDT") {
          Trk::Volume* mdtVol;
	  Trk::TrapezoidVolumeBounds* mdtBounds; 
          if (lowX == currX) {
	    mdtBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,compTrdBounds->halflengthZ());
            mdtVol = new Trk::Volume(new HepTransform3D(compVol[i]->transform()),mdtBounds);
	  } else {
	    std::cout << "lowX,currX:" << lowX <<","<<currX << std::endl;
	    std::cout << "increase the basic Mdt volume" << std::endl;
          }
          const Trk::TrackingVolume* mdtTrkVol = processMdtTrd(mdtVol,compGeo[i],compTransf[i]);
          trkVols.push_back(mdtTrkVol); 
          volSteps->push_back(2*mdtBounds->halflengthZ());
          currX += 2.*mdtBounds->halflengthZ();
          comp_processed = true;
        }
        if ( !comp_processed ) std::cout << "unknown technology:" <<compName[i]<<std::endl;
      } // end loop over station children

      // there may be a spacer still open
      if (openSpacer) {
        if (maxX >= currX+spacerlowXsize+spaceruppXsize) {
          Trk::TrapezoidVolumeBounds* spacerBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,0.5*(maxX-currX));
          
          /*
	  Trk::Volume* spacerVol =new Trk::Volume(new HepTransform3D(HepRotateZ3D(90*deg)*HepRotateX3D(90*deg)
                                                           * HepTranslateZ3D(currX+spacerBounds->halflengthZ())),spacerBounds);
	  */
	  Trk::Volume* spacerVol =new Trk::Volume(new HepTransform3D(HepRotateY3D(90*deg)*HepRotateZ3D(90*deg)
                                                           * HepTranslateZ3D(currX+spacerBounds->halflengthZ())),spacerBounds);
          const Trk::TrackingVolume* spacerTrkVol = processSpacer(spacerVol ,geoSpacer, transfSpacer);
          trkVols.push_back(spacerTrkVol);  
          volSteps->push_back(maxX-currX);
          currX = maxX;
          openSpacer = false;
        } else {
	  std::cout <<"currX,maxX,lowX,uppX:"<< currX<<"," << maxX <<","<<spacerlowXsize<<"," <<spaceruppXsize<< std::endl; 
	  std::cout << "clash in spacer definition!(last volume)" << std::endl; 
        }          
      }
      // there may be an Rpc still open
      if (openRpc) {
	std::cout << "RPC not coded for trapezoid chambers"<<std::endl;
	std::cout << "maxX, other:"<< maxX <<"," << currX+rpclowXsize+rpcuppXsize << std::endl;
        if (maxX >= currX+rpclowXsize+rpcuppXsize) {
          Trk::CuboidVolumeBounds* rpcBounds = new Trk::CuboidVolumeBounds(0.5*(maxX-currX),envY,envZ);
	  Trk::Volume* rpcVol =new Trk::Volume(new HepTranslate3D(currX+rpcBounds->halflengthZ(),0.,0.),rpcBounds);
          const Trk::TrackingVolume* rpcTrkVol = processRpc(rpcVol ,geoRpc,transfRpc);
          trkVols.push_back(rpcTrkVol);  
          volSteps->push_back(maxX-currX);
          currX = maxX;
          openRpc = false;
        } else {
	  std::cout << "clash in Rpc definition!(last volume)" << std::endl; 
        }          
      }
      // create VolumeArray (1DX) 
      const Trk::TrackingVolumeArray* components = 0;
      const std::vector<const Trk::TrackingVolume*>* compVols = new std::vector<const Trk::TrackingVolume*>( trkVols );
      Trk::BinUtility* binUtility = new Trk::BinUtility1DX( -( envelope->halflengthZ() ), volSteps);
      if (m_trackingVolumeArrayCreator)  components = m_trackingVolumeArrayCreator->trapezoidVolumesArrayNav( *compVols, binUtility, false);
      // std::cout << "tracking volume array created" << std::endl;
     
   return components;  
}

// finalize
StatusCode Muon::MuonStationTypeBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    return StatusCode::SUCCESS;
}
//
const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processMdtBox(Trk::Volume*& vol,const GeoVPhysVol*& gv, HepTransform3D*& transf) const
{
  // tolerance
  double tol = 0.001;
  // std::cout << "processing MDT, number of children volumes:"<< gv->getNChildVols() <<std::endl; 
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::MaterialProperties*> x_mat;
  std::vector<double> x_thickness;
  double currX = -100000; 
  // here one could save time by not reading all tubes  
  for (unsigned int ich =0; ich< gv->getNChildVols(); ++ich) {
     const GeoVPhysVol* cv = &(*(gv->getChildVol(ich))); 
     const GeoLogVol* clv = cv->getLogVol();
     HepTransform3D transfc = gv->getXToChildVol(ich);        
     // std::cout << "MDT component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transfc.getTranslation()<<std::endl;
     // printChildren(cv);
     Trk::MaterialProperties* mdtMat=0;
       double xv = 0.;
       if ( (clv->getName()).substr(0,3)=="MDT") {
           xv = 13.0055;  // the half-thickness
           if (!m_mdtTubeMat ) {
	      const GeoTube* tube=dynamic_cast<const GeoTube*> (clv->getShape());
              double vol = 8*(tube->getRMax())*(tube->getZHalfLength())*xv;
	      // std::cout << " part of layer volume assigned to 1 tube:" << vol << std::endl;
	      // std::cout << "tube dimensions:" << tube->getRMax() << "," << tube->getRMin() << "," << tube->getZHalfLength() << std::endl;
              m_mdtTubeMat = getAveragedLayerMaterial(cv,vol,2*xv); 
           }        
           mdtMat = m_mdtTubeMat;
       }
       if ( (clv->getName())=="MultiLayerFoam") {
	 // std::cout << "processing MultiLayerFoam" << std::endl;  
          const GeoTrd* trd =  dynamic_cast<const GeoTrd*> (clv->getShape());
          if (trd){
            xv = trd->getXHalfLength1();
          } else {
	    std::cout << "MDT MultiFoam not trapezoid in MDT box chamber" << std::endl;
          }
          if (!m_mdtFoamMat && trd ) {
            double vol = 8*(trd->getXHalfLength1())*(trd->getYHalfLength1())*(trd->getZHalfLength());
            m_mdtFoamMat = getAveragedLayerMaterial(cv,vol,2*xv); 
          }
          if (!trd && m_mdtFoamMat) {
            xv = 0.5*m_mdtFoamMat->thickness();
          }
          mdtMat = m_mdtFoamMat;        
       }
       if (  transfc.getTranslation()[0] != currX ) {
         if (x_array.size() == 0 || transfc.getTranslation()[0] > x_array.back() ) {
           x_array.push_back(  transfc.getTranslation()[0] );
           x_mat.push_back(mdtMat);
           x_thickness.push_back( 2*xv );
           currX = transfc.getTranslation()[0];
	   // std::cout << "layer info included:" << clv->getName()<<"," << 2*xv <<","<< currX<< std::endl; 
         } else {
	  std::vector<double>::iterator xIter=x_array.begin();
	  std::vector<Trk::MaterialProperties*>::iterator mIter=x_mat.begin();
	  std::vector<double>::iterator tIter=x_thickness.begin();
          while ( transfc.getTranslation()[0] > *xIter ) {xIter++;mIter++;}
          x_array.insert(xIter,transfc.getTranslation()[0]);
          x_mat.insert(mIter,mdtMat);
          x_thickness.insert(tIter,2*xv);
          currX = transfc.getTranslation()[0];
        }
       }     
     
  }
  // create layers // 
  const Trk::PlaneLayer* layer;
  double thickness=0.;
  Trk::OverlapDescriptor* od=0;
  const Trk::CuboidVolumeBounds* volBounds = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(vol->volumeBounds()));
  if ( volBounds ) {
    double yv = volBounds->halflengthY();
    double zv = volBounds->halflengthZ();
    Trk::RectangleBounds* bounds=0;
    for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
      // x-y plane -> y-z plane
      bounds = new Trk::RectangleBounds(yv,zv); 
      thickness = x_thickness[iloop];
      Trk::MaterialProperties material=m_muonMaterial;
      if ( x_mat[iloop] ) material = *(x_mat[iloop]);      

      HepTransform3D* cTr = new HepTransform3D( (*transf) * HepTranslateX3D(x_array[iloop])
                                               * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      Trk::HomogenousLayerMaterial mdtMaterial(material, Trk::oppositePre);  
      layer = new Trk::PlaneLayer(cTr,
                                  bounds,
                                  mdtMaterial,
                                  thickness,
                                  od );
      layers.push_back(layer);
      // std::cout << "layer built ok"<<std::endl;
    }
  } 
  // create the BinnedArray
  std::cout << "number of Mdt layers:"<<layers.size()<<std::endl;
  std::vector<LayTr> layerOrder;
  std::vector<double> binSteps;
  // check if additional (navigation) layers needed
  double minX = - volBounds->halflengthX();
  if (layers.size()) {
     minX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = minX; 
     for (unsigned int i=0;i<layers.size();i++) { 
       const HepTransform3D* ltransf = new HepTransform3D(layers[i]->transform());
       layerOrder.push_back(LayTr(Trk::SharedObject<const Trk::Layer>(layers[i]), ltransf ));
       binSteps.push_back(ltransf->getTranslation()[0]+0.5*layers[i]->thickness()-currX);
       currX = ltransf->getTranslation()[0]+0.5*layers[i]->thickness();     
     }
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( minX, new std::vector<double>(binSteps));
  Trk::LayerArray* mdtLayerArray = 0;
  mdtLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     
  std::string name="MDT";
  const Trk::TrackingVolume* mdt= new Trk::TrackingVolume(*vol,
                                                          m_muonMaterial,
                                                          m_muonMagneticField,
                                                          mdtLayerArray,0,
                                                          name);         
  return mdt;
}
//
const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processMdtTrd(Trk::Volume*& vol,const GeoVPhysVol*& gv, HepTransform3D*& transf) const
{
  // std::cout << "processing MDT, number of children volumes:"<< gv->getNChildVols() <<std::endl; 
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::MaterialProperties*> x_mat;
  std::vector<double> x_thickness;
  double currX = -100000; 
  for (unsigned int ich =0; ich< gv->getNChildVols(); ++ich) {
     const GeoVPhysVol* cv = &(*(gv->getChildVol(ich))); 
     const GeoLogVol* clv = cv->getLogVol();
     HepTransform3D transfc = gv->getXToChildVol(ich);        
     //std::cout << "MDT component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transfc.getTranslation()<<std::endl;
       double xv = 0.;
       if (clv->getShape()->type()=="Trd"){
          const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
          double x1v = trd->getXHalfLength1();
          double x2v = trd->getXHalfLength2();
          double y1v = trd->getYHalfLength1();
          double y2v = trd->getYHalfLength2();
          double zv = trd->getZHalfLength();
          // std::cout << "dimensions:"<<x1v<<","<<x2v<<","<<y1v<<","<<y2v<<","<<zv<<std::endl; 
          if ( x1v==x2v ) xv = x1v;
       }
       Trk::MaterialProperties* mdtMat=0; 
       if ( (clv->getName()).substr(0,3)=="MDT") {
           xv = 13.0055;  // the half-thickness
           if (!m_mdtTubeMat ) {
	      const GeoTube* tube=dynamic_cast<const GeoTube*> (clv->getShape());
              double vol = 8*(tube->getRMax())*(tube->getZHalfLength())*xv;
	      // std::cout << " part of layer volume assigned to 1 tube:" << vol << std::endl;
	      // std::cout << "tube dimensions:" << tube->getRMax() << "," << tube->getRMin() << "," << tube->getZHalfLength() << std::endl;
              m_mdtTubeMat = getAveragedLayerMaterial(cv,vol,2*xv); 
           }        
           mdtMat = m_mdtTubeMat;
       }
       if ( (clv->getName())=="MultiLayerFoam") {
	 // std::cout << "processing MultiLayerFoam" << std::endl;  
          const GeoTrd* trd =  dynamic_cast<const GeoTrd*> (clv->getShape());
          if (trd){
            xv = trd->getXHalfLength1();
          } else {
	    std::cout << "MDT MultiFoam not trapezoid in MDT box chamber" << std::endl;
          }
          if (!m_mdtFoamMat && trd ) {
            double vol = 8*(trd->getXHalfLength1())*(trd->getYHalfLength1())*(trd->getZHalfLength());
            m_mdtFoamMat = getAveragedLayerMaterial(cv,vol,2*xv); 
          }
          if (!trd && m_mdtFoamMat) {
            xv = 0.5*m_mdtFoamMat->thickness();
          }
          mdtMat = m_mdtFoamMat;        
       }

       if (  transfc.getTranslation()[0] != currX ) {
         if (x_array.size() == 0 || transfc.getTranslation()[0] > x_array.back() ) {
           x_array.push_back(  transfc.getTranslation()[0] );
           x_mat.push_back( mdtMat );
           x_thickness.push_back( 2*xv );
           currX = transfc.getTranslation()[0];
         } else {
	  std::vector<double>::iterator xIter=x_array.begin();
	  std::vector<Trk::MaterialProperties*>::iterator mIter=x_mat.begin();
	  std::vector<double>::iterator tIter=x_thickness.begin();
          while ( transfc.getTranslation()[0] > *xIter ) {xIter++;mIter++;}
          x_array.insert(xIter,transfc.getTranslation()[0]);
          x_mat.insert(mIter,mdtMat);
          x_thickness.insert(tIter,2*xv);
          currX = transfc.getTranslation()[0];
        }
       }     
  }
  // create layers // 
  const Trk::PlaneLayer* layer;
  double thickness=0.;
  Trk::OverlapDescriptor* od=0;
  const Trk::TrapezoidVolumeBounds* volBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(vol->volumeBounds()));
  if ( volBounds ) {
    double x1v = volBounds->minHalflengthX();
    double x2v = volBounds->maxHalflengthX();
    double yv = volBounds->halflengthY();
    Trk::TrapezoidBounds* bounds=0;
    for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
      // std::cout << iloop << ","<< x_array[iloop] << "," << x_mat[iloop]<<std::endl;
      // x-y plane -> y-z plane
      bounds = new Trk::TrapezoidBounds(x1v,x2v,yv); 
      thickness = x_thickness[iloop];
      Trk::MaterialProperties material = m_muonMaterial;
      if (x_mat[iloop]) material = *(x_mat[iloop]);
      Trk::HomogenousLayerMaterial mdtMaterial(material, Trk::oppositePre);
      /*
      HepTransform3D* cTr = new HepTransform3D( (*transf) * HepRotateZ3D(-90*deg)*HepTranslateX3D(x_array[iloop]) * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      */
      HepTransform3D* cTr = new HepTransform3D( (*transf) * HepTranslateZ3D(x_array[iloop]) );
      layer = new Trk::PlaneLayer(cTr,
                                  bounds,
                                  mdtMaterial,
                                  thickness,
                                  od );
      layers.push_back(layer);
      // std::cout << "MDT layer position: X or Y ordering? "<<cTr->getTranslation()<<std::endl;
    }
  } 

  // create the BinnedArray
  std::cout << "number of Mdt layers:"<<layers.size()<<std::endl;
  std::vector<LayTr> layerOrder;
  std::vector<double> binSteps;
  // 
  double minX = - volBounds->halflengthZ();
  if (layers.size()) {
     minX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = minX; 
     for (unsigned int i=0;i<layers.size();i++) { 
       const HepTransform3D* ltransf = new HepTransform3D(layers[i]->transform());
       layerOrder.push_back(LayTr(Trk::SharedObject<const Trk::Layer>(layers[i]), ltransf ));
       binSteps.push_back(ltransf->getTranslation()[0]+0.5*layers[i]->thickness()-currX);
       currX = ltransf->getTranslation()[0]+0.5*layers[i]->thickness();     
     }
  }
  /*
  for (unsigned int i=0;i<layers.size();i++) { 
    const HepTransform3D* ltransf = new HepTransform3D(layers[i]->transform());
    layerOrder.push_back(LayTr(Trk::SharedObject<const Trk::Layer>(layers[i]), ltransf ));
    if (i>0) {
      binSteps.push_back(ltransf->getTranslation()[0] -currX);
    } else {
      minX = ltransf->getTranslation()[0];
    } 
    currX = ltransf->getTranslation()[0];
  }
  */

  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( minX, new std::vector<double>(binSteps));
  Trk::LayerArray* mdtLayerArray = 0;
  mdtLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     
  std::string name="MDT";
  const Trk::TrackingVolume* mdt= new Trk::TrackingVolume(*vol,
                                                          m_muonMaterial,
                                                          m_muonMagneticField,
                                                          mdtLayerArray,0,
                                                          name);         
  // std::cout << "Mdt processed with" << layers.size() << " layers" << std::endl;
  return mdt;
}
//
const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processRpc(Trk::Volume*& vol,std::vector<const GeoVPhysVol*> gv, std::vector<HepTransform3D*> transfc) const
{
  // layers correspond to DedModules and RpcModules; all substructures averaged in material properties
  std::vector<const Trk::Layer*> layers;
  for (unsigned int ic=0; ic<gv.size(); ++ic) {
    // std::cout << "processing Rpc component:"<< gv[ic]->getLogVol()->getName() <<std::endl;
    const GeoLogVol* glv = gv[ic]->getLogVol();
    if (glv->getShape()->type()=="Box") {
      const GeoBox* box = dynamic_cast<const GeoBox*> (glv->getShape());
      double xs = box->getXHalfLength();
      double ys = box->getYHalfLength();
      double zs = box->getZHalfLength();
      // std::cout << "dimensions:"<<box->getXHalfLength() << ","<<box->getYHalfLength() << ","<<box->getZHalfLength() << std::endl;
      // translating into layer; x dimension defines thickness
      const Trk::PlaneLayer* layer;
      double thickness=2*xs;
      Trk::OverlapDescriptor* od=0;
      Trk::RectangleBounds* bounds = new Trk::RectangleBounds(ys,zs); 
      HepTransform3D* cTr = new HepTransform3D((*transfc[ic]) * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      Trk::MaterialProperties rpcMat = m_muonMaterial;               // default
      if ( (glv->getName()).substr(0,3)=="Ded" ) {
        if (thickness == 50.0) {
          if (!m_rpcDed50) {
            double vol = 8*xs*ys*zs;
            m_rpcDed50 = getAveragedLayerMaterial(gv[ic],vol,2*xs);
          }
          rpcMat=*m_rpcDed50;  
        } else { std::cout << "Ded thickness different from 50:" << thickness << std::endl; }
      } else {
        if (thickness == 46.0) {
          if (!m_rpc46) {
            double vol = 8*xs*ys*zs;
            m_rpc46 = getAveragedLayerMaterial(gv[ic],vol,2*xs);
          }
          rpcMat=*m_rpc46;  
        } else { std::cout << "RPC module thickness different from 46:" << thickness << std::endl; }
      }
          
      Trk::HomogenousLayerMaterial rpcMaterial(rpcMat, Trk::oppositePre);
      layer = new Trk::PlaneLayer(cTr,
                                     bounds,
                                     rpcMaterial,
                                     thickness,
                                     od );
      layers.push_back(layer);
    }
    if (glv->getShape()->type()=="Trd") {
      const GeoTrd* trd = dynamic_cast<const GeoTrd*> (glv->getShape());
      double xs1 = trd->getXHalfLength1();
      double xs2 = trd->getXHalfLength2();
      double ys1 = trd->getYHalfLength1();
      double ys2 = trd->getYHalfLength2();
      double zs = trd->getZHalfLength();
      // translating into layer; x dimension defines thickness
      if (xs1==xs2 && ys1==ys2) {
        const Trk::PlaneLayer* layer;
        double thickness=2*xs1;
        Trk::OverlapDescriptor* od=0;
        Trk::RectangleBounds* bounds = new Trk::RectangleBounds(ys1,zs); 
        HepTransform3D* cTr = new HepTransform3D((*transfc[ic]) * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
        Trk::MaterialProperties rpcMat = m_muonMaterial;               // default
        if ( (glv->getName()).substr(0,3)=="Ded" ) {
          if (thickness == 50.0) {
            if (!m_rpcDed50) {
              double vol = 8*xs1*ys1*zs;
             m_rpcDed50 = getAveragedLayerMaterial(gv[ic],vol,2*xs1);
            }
            rpcMat=*m_rpcDed50;  
          } else { std::cout << "Ded thickness different from 50:" << thickness << std::endl; }
        } else {
          if (thickness == 46.0) {
            if (!m_rpc46) {
              double vol = 8*xs1*ys1*zs;
              m_rpc46 = getAveragedLayerMaterial(gv[ic],vol,2*xs1);
            }
            rpcMat=*m_rpc46;  
          } else { std::cout << "RPC module thickness different from 46:" << thickness << std::endl; }
        }
          
        Trk::HomogenousLayerMaterial rpcMaterial(rpcMat, Trk::oppositePre);
        layer = new Trk::PlaneLayer(cTr,
                                     bounds,
                                     rpcMaterial,
                                     thickness,
                                     od );
        layers.push_back(layer);
      } else {
        std::cout << "RPC true trapezoid layer, not coded yet" <<std::endl;
      }
      /*
      const Trk::PlaneLayer* layer;
      double thickness=2*xs;
      Trk::OverlapDescriptor* od=0;
      Trk::RectangleBounds* bounds = new Trk::RectangleBounds(ys,zs); 
      HepTransform3D* cTr = new HepTransform3D((*transfc[ic]) * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      const Trk::MaterialProperties material = getLayerMaterial(glv->getName(),thickness);
      Trk::LayerMaterialProperties rpcMaterial(material, Trk::oppositePre);
      layer = new Trk::PlaneLayer(cTr,
                                     bounds,
                                     rpcMaterial,
                                     thickness,
                                     od );
      layers.push_back(layer);
      */
    }
  } // end loop over Modules

  std::vector<const Trk::Layer*>* rpcLayers = new std::vector<const Trk::Layer*>(layers); 

  std::string name="RPC";
  const Trk::TrackingVolume* rpc= new Trk::TrackingVolume(*vol,
                                                          m_muonMaterial,
                                                          m_muonMagneticField,
                                                          rpcLayers,
                                                          name);         
  // std::cout << "Rpc processed with" << layers.size() << " layers" << std::endl;
  return rpc;
}
//

const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processSpacer(Trk::Volume*& vol,std::vector<const GeoVPhysVol*> gv, std::vector<HepTransform3D*> transf) const
{
  // spacers: one level below, assumed boxes
  std::vector<const Trk::Layer*> layers;
  for (unsigned int ic=0; ic<gv.size(); ++ic) {
    // std::cout << "processing spacer, number of children volumes:"<< gv[ic]->getNChildVols() <<std::endl; 
    for (unsigned int ich =0; ich< gv[ic]->getNChildVols(); ++ich) {
      const GeoVPhysVol* cv = &(*(gv[ic]->getChildVol(ich))); 
      const GeoLogVol* clv = cv->getLogVol();
      HepTransform3D transform = gv[ic]->getXToChildVol(ich);        
      // std::cout << "Spacer component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transf[ic]->getTranslation()<<","<<transform.getTranslation()<<std::endl;
      // std::cout << "combined transform:" << ((*transf[ic])*transform).getTranslation()<<std::endl; 
      if (clv->getShape()->type()=="Box") {
         const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
         double xs = box->getXHalfLength();
         double ys = box->getYHalfLength();
         double zs = box->getZHalfLength();
         // std::cout << "dimensions:"<<box->getXHalfLength() << ","<<box->getYHalfLength() << ","<<box->getZHalfLength() << std::endl;
         // translating into layer; find minimal size
	 const Trk::PlaneLayer* layer;
	 Trk::RectangleBounds* bounds;
         double thickness=0.;
         Trk::OverlapDescriptor* od=0;
         HepTransform3D* cTr = new HepTransform3D((*transf[ic])*transform);
         if (zs <= xs && zs <= ys ) { // x-y plane
	   bounds = new Trk::RectangleBounds(xs,ys); 
           thickness = zs;
         } else {
           if (xs <= ys && xs <= zs ) { // x-y plane -> y-z plane
	     bounds = new Trk::RectangleBounds(ys,zs); 
             thickness = xs;
             cTr = new HepTransform3D((*cTr) * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
           } else {  // x-y plane -> x-z plane
	     bounds = new Trk::RectangleBounds(xs,zs); 
             thickness = ys;
             cTr = new HepTransform3D((*cTr) * HepRotateX3D(90*deg));
           }
         } 
	 const Trk::MaterialProperties material = getLayerMaterial(clv->getMaterial()->getName(),thickness);
         Trk::HomogenousLayerMaterial spacerMaterial(material, Trk::oppositePre);

         layer = new Trk::PlaneLayer(cTr,
                                     bounds,
                                     spacerMaterial,
                                     thickness,
                                     od );
         layers.push_back(layer);
      }
    }
  }

  std::vector<const Trk::Layer*>* spacerLayers = new std::vector<const Trk::Layer*>(layers); 

  std::string name="Spacer";
  const Trk::TrackingVolume* spacer= new Trk::TrackingVolume(*vol,
                                                          m_muonMaterial,
                                                          m_muonMagneticField,
                                                          spacerLayers,
                                                          name);         
  std::cout << "spacer processed with" << layers.size() << " layers" << std::endl;
  return spacer;
}

const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processCscStation(const GeoVPhysVol* mv, std::string name) const
{

 // CSC stations have the particularity of displacement in Z between multilayer and the spacer - the envelope
 //   has to be derived from the component volume shape and component displacement
    bool isDiamond = false;
    double xMin=0; double xMed=0; double xMax=0; double y1=0; double y2=0; double z=0;
 //   std::cout << "processing CSC, number of children volumes:"<< mv->getNChildVols() <<std::endl;
 // find the shape and dimensions for the first component
    const GeoVPhysVol* cv = &(*(mv->getChildVol(0))); 
    const GeoLogVol* clv = cv->getLogVol();
    HepTransform3D transform = mv->getXToChildVol(0);        
  //   std::cout << "First CSC component:"<< clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transform.getTranslation()<<std::endl;
    if (clv->getShape()->type()=="Shift") {
      const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (clv->getShape());
      //std::cout << shift->getOp()->type()<<","<<(shift->getX()).getTranslation()<<std::endl;
      if (shift->getOp()->type()=="Union") {
        // that would be the union making the diamond/double trapezoid shape, let's retrieve the parameters 
        isDiamond = true;
   	const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion*> (shift->getOp());
        //std::cout << uni->getOpA()->type()<<","<< uni->getOpB()->type() << std::endl;
        if (uni->getOpA()->type()=="Trd") {
         const GeoTrd* trdA = dynamic_cast<const GeoTrd*> (uni->getOpA());
         
         double xs1 = trdA->getXHalfLength1();
         double xs2 = trdA->getXHalfLength2();
         double ys1 = trdA->getYHalfLength1();
         double ys2 = trdA->getYHalfLength2();
         double zs = trdA->getZHalfLength();
         //std::cout << "dimensionsA:"<<xs1<<","<<xs2<<","<<ys1<<","<< ys2 << ","<< zs << std::endl;
         
         xMin =  trdA->getYHalfLength1();
         xMed =  trdA->getYHalfLength2();
         y1   =  trdA->getZHalfLength();
         z    =  trdA->getXHalfLength1();
        }
        if (uni->getOpB()->type()=="Shift") {
          const GeoShapeShift* sh = dynamic_cast<const GeoShapeShift*> (uni->getOpB());
          std::cout << sh->getOp()->type()<<","<<(sh->getX()).getTranslation()<<std::endl;
          const GeoTrd* trdB = dynamic_cast<const GeoTrd*> (sh->getOp());
          
          double xs1 = trdB->getXHalfLength1();
          double xs2 = trdB->getXHalfLength2();
          double ys1 = trdB->getYHalfLength1();
          double ys2 = trdB->getYHalfLength2();
          double zs = trdB->getZHalfLength();
          //std::cout << "dimensionsB:"<<xs1<<","<<xs2<<","<<ys1<<","<< ys2 << ","<< zs << std::endl;
          
          if ( trdB->getYHalfLength1() != xMed ||  trdB->getXHalfLength1() != z )
                     std::cout <<"Something is wrong: dimensions of 2 trapezoids do not match"<<std::endl ;
          xMax =  trdB->getYHalfLength2();
          y2   =  trdB->getZHalfLength();
       }
      } //end Union
      if (shift->getOp()->type()=="Trd") {
        // that would be the trapezoid shape, let's retrieve the parameters 
         const GeoTrd* trd = dynamic_cast<const GeoTrd*> (shift->getOp());
         xMin =  trd->getYHalfLength1();
         xMed =  trd->getYHalfLength2();
         y1   =  trd->getZHalfLength();
         z    =  trd->getXHalfLength1();      
      } //end Trd
    } else {
      if (clv->getShape()->type()=="Trd"){
        // that would be the trapezoid shape, let's retrieve the parameters 
         const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
         xMin =  trd->getYHalfLength1();
         xMed =  trd->getYHalfLength2();
         y1   =  trd->getZHalfLength();
         z    =  trd->getXHalfLength1();      
      }
    }
// then loop over all components to get total Xsize & transforms 
   std::vector<HepTransform3D> compTransf;
   std::vector<std::string> compName;
   std::vector<const GeoVPhysVol*> compGeoVol;
   std::vector<double> xSizes;
   for (unsigned int ich =0; ich< mv->getNChildVols(); ++ich) {
     const GeoVPhysVol* cv = &(*(mv->getChildVol(ich))); 
     const GeoLogVol* clv = cv->getLogVol();
     HepTransform3D transform = mv->getXToChildVol(ich);        
     compTransf.push_back(transform);
     compName.push_back(clv->getName());
     compGeoVol.push_back(cv);
     //std::cout << "CSC component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transform.getTranslation()<<std::endl;
      if (clv->getShape()->type()=="Shift") {
       const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (clv->getShape());
       //std::cout << shift->getOp()->type()<<","<<(shift->getX()).getTranslation()<<std::endl;
       if (shift->getOp()->type()=="Union") {
	 // that would be the union making the diamond/double trapezoid shape, let's retrieve the parameters 
	 double xMin=0; double xMed=0; double xMax=0; double y1=0; double y2=0; double z=0;
   	 const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion*> (shift->getOp());
         //std::cout << uni->getOpA()->type()<<","<< uni->getOpB()->type() << std::endl;
         if (uni->getOpA()->type()=="Trd") {
          const GeoTrd* trdA = dynamic_cast<const GeoTrd*> (uni->getOpA());
          double xSize  =  trdA->getXHalfLength1();
          xSizes.push_back(xSize); 
        }
       } //end Union
      } // end Shift
      if (clv->getShape()->type()=="Trd") {
       const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
       double xSize  =  trd->getXHalfLength1();
       xSizes.push_back(xSize); 
      } // end Trd
      
      // printChildren(cv);  
   } 
   // this should be enough to build station envelope
   double xTotal = 0;
   for (unsigned int i=0;i<xSizes.size(); i++) xTotal += xSizes[i];
   //std::cout << "total station thickness:" << xTotal << std::endl;
   double zShift = 0;
   zShift = fabs(((compTransf.front()).getTranslation())[2])+  fabs(((compTransf.back()).getTranslation())[2]);
   //std::cout << "z displacement:" << zShift << std::endl;
   // calculate displacement with respect to GeoModel station volume
   // one way or the other, the station envelope is double trapezoid
   Trk::Volume* envelope;
   double envXMed = xMed;
   double envY1   = y1; 
   double envY2   = y2;
   std::vector<double>* volSteps=new std::vector<double>; 
   std::vector<const Trk::TrackingVolume*> components;
   if ( !isDiamond ) {
     Trk::TrapezoidVolumeBounds* cscBounds=0;
     Trk::TrapezoidVolumeBounds* compBounds=0;
     xMax = xMed;
     y2 = 0.5*zShift;
     //std::cout << "envelope dimensions:"<<xMin<<","<<envXMed<<","<<xMax<<","<<envY1<<","<<envY2<<","<<xTotal<<std::endl;
     // cscBounds = new Trk::DoubleTrapezoidVolumeBounds(xMin,xMed,xMax,y1,y2,xTotal); 
     cscBounds = new Trk::TrapezoidVolumeBounds(xMin,xMax,y1,xTotal); 
     //std::cout << "envelope bounds created"<<std::endl;
     // xy -> yz  rotation
     // the center of Volume is shifted by y1-y2 in y
     HepTransform3D* cTr = new HepTransform3D( HepRotateY3D(90*deg)*HepRotateZ3D(90*deg) );
     envelope = new Trk::Volume(cTr,cscBounds);
     //std::cout << "envelope created"<<std::endl;
     // components
     double xCurr = -xTotal;
     for (unsigned int ic = 0; ic< xSizes.size(); ic++) {
       // component volumes follow the envelope dimension
       xCurr += xSizes[ic];
       HepTransform3D* compTr = new HepTransform3D( HepRotateY3D(90*deg)*HepRotateZ3D(90*deg)*HepTranslateZ3D(xCurr));
       compBounds = new Trk::TrapezoidVolumeBounds(xMin,xMax,y1,xSizes[ic]);
       const Trk::LayerArray* cscLayerArray = processCSCTrdComponent(compGeoVol[ic],compBounds,compTr); 
       Trk::Volume* compVol = new Trk::Volume(compTr,compBounds);
       const Trk::TrackingVolume* compTV = new Trk::TrackingVolume( *compVol,
                                                                    m_muonMaterial,
                                                                    m_muonMagneticField,
                                                                    cscLayerArray,0,
                                                                    compName[ic]);                    
       components.push_back(compTV);
       volSteps->push_back(2*xSizes[ic]);  
       xCurr += xSizes[ic];
     } // end components
     
   } else {
     Trk::DoubleTrapezoidVolumeBounds* cscBounds=0;
     Trk::DoubleTrapezoidVolumeBounds* compBounds=0;
     if (xMed!=xMin && xMed!=xMax) {
       envXMed += zShift/(y1/(xMed-xMin)+y2/(xMed-xMax));
       envY1 = y1*(envXMed-xMin)/(xMed-xMin);
       envY2 = y2*(envXMed-xMax)/(xMed-xMax);
     } 
     //std::cout << "envelope dimensions:"<<xMin<<","<<envXMed<<","<<xMax<<","<<envY1<<","<<envY2<<","<<z+0.5*zShift<<std::endl;
     cscBounds = new Trk::DoubleTrapezoidVolumeBounds(xMin,envXMed,xMax,envY1,envY2,xTotal); 
     //std::cout << "envelope bounds created"<<std::endl;
     // xy -> yz  rotation
     // the center of DoubleTrapezoidVolume is shifted by (envY1-envY2) in y
     //HepTransform3D* cTr = new HepTransform3D(HepRotateZ3D(90*deg)*HepTranslateY3D(envY1-envY2));
     HepTransform3D* cTr = new HepTransform3D(HepRotateZ3D(90*deg)*HepTranslateZ3D(envY1-envY2));
     envelope = new Trk::Volume(cTr,cscBounds);
     //std::cout << "envelope created"<<std::endl;
     // components
     double xCurr = -xTotal;
     for (unsigned int ic = 0; ic< xSizes.size(); ic++) {
       // component volumes follow the envelope dimension
       xCurr += xSizes[ic];
       HepTransform3D* compTr = new HepTransform3D(HepRotateZ3D(90*deg)*HepTranslateZ3D(envY1-envY2)*HepTranslateY3D(xCurr));
       compBounds = new Trk::DoubleTrapezoidVolumeBounds(xMin,envXMed,xMax,envY1,envY2,xSizes[ic]); 
       const Trk::LayerArray* cscLayerArray = processCSCDiamondComponent(compGeoVol[ic],compBounds,compTr); 
       Trk::Volume* compVol = new Trk::Volume(compTr,compBounds);
       const Trk::TrackingVolume* compTV = new Trk::TrackingVolume( *compVol,
                                                                    m_muonMaterial,
                                                                    m_muonMagneticField,
                                                                    cscLayerArray,0,
                                                                    compName[ic]);                    
       components.push_back(compTV);
       volSteps->push_back(2*xSizes[ic]);  
       xCurr += xSizes[ic];
     } // end components
   }

 // convert component volumes into array 
   const Trk::BinnedArray<Trk::TrackingVolume>* compArray = 0; 
   if (components.size() && isDiamond) {
     Trk::BinUtility1DX* binUtil = new Trk::BinUtility1DX(-xTotal, volSteps);
     if (m_trackingVolumeArrayCreator) compArray = m_trackingVolumeArrayCreator->doubleTrapezoidVolumesArrayNav( components, binUtil, false);
   }
   if (components.size() && !isDiamond) {
     Trk::BinUtility1DX* binUtil = new Trk::BinUtility1DX(-xTotal, volSteps);
     if (m_trackingVolumeArrayCreator) compArray = m_trackingVolumeArrayCreator->trapezoidVolumesArrayNav( components, binUtil, false);
   }
 // ready to build the station prototype
 const Trk::TrackingVolume* csc_station = new Trk::TrackingVolume( *envelope,
                                                                 m_muonMaterial,
                                                                 m_muonMagneticField,
                                                                 0,compArray,
                                                                 name);                    
 return csc_station;    
}

std::vector<const Trk::TrackingVolume*> Muon::MuonStationTypeBuilder::processTgcStation(const GeoVPhysVol* mv) const
{
 // TGC stations 
  std::vector<const Trk::TrackingVolume*> tgc_stations;
 //  printChildren(mv);
  Trk::TrapezoidVolumeBounds* tgcBounds;
  Trk::Volume* envelope;
  for (unsigned int ich =0; ich< mv->getNChildVols(); ++ich) {
     const GeoVPhysVol* cv = &(*(mv->getChildVol(ich))); 
     const GeoLogVol* clv = cv->getLogVol();
     std::string tgc_name = clv->getName();
     HepTransform3D transform = mv->getXToChildVol(ich);        
     //std::cout << "TGC component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transform.getTranslation()<<std::endl;
     if (clv->getShape()->type()=="Trd") {
       const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
       double x1 = trd->getXHalfLength1();
       double x2 = trd->getXHalfLength2();
       double y1 = trd->getYHalfLength1();
       double y2 = trd->getYHalfLength2();
       double z = trd->getZHalfLength();
       //std::cout << "dimensions:"<<x1<<","<<x2<<","<<y1<<","<< y2 << ","<< z << std::endl;
       // define envelope
       tgcBounds = new Trk::TrapezoidVolumeBounds(y1,y2,z,x1); 
       // xy -> yz  rotation
       HepTransform3D* tTr = new HepTransform3D( transform * HepRotateY3D(90*deg)* HepRotateZ3D(90*deg) );
       envelope = new Trk::Volume(tTr,tgcBounds);
       const Trk::LayerArray* tgcLayerArray = processTGCComponent(cv,tgcBounds, tTr); 
       // ready to build the station prototype
       const Trk::TrackingVolume* tgc_station = new Trk::TrackingVolume( *envelope,
                                                                        m_muonMaterial,
                                                                        m_muonMagneticField,
									tgcLayerArray,0,
                                                                        tgc_name);                    

       if (tgc_station) tgc_stations.push_back(tgc_station);       
     } else {
       std::cout << "TGC component not trapezoid ?" << std::endl;
     }
   }
  
 return tgc_stations;    

}

const Trk::MaterialProperties Muon::MuonStationTypeBuilder::getLayerMaterial(std::string mat, double thickness) const
{
  const Trk::MaterialProperties* matprop;
  if (mat == "Aluminium") matprop=new Trk::MaterialProperties(thickness, 89.8689, 26.98154, 13., 2.7);         
  return *matprop;
}

const void Muon::MuonStationTypeBuilder::printChildren(const GeoVPhysVol* pv) const
{
  // subcomponents
  unsigned int nc = pv->getNChildVols();
  for (unsigned int ic=0; ic<nc; ic++) {
    HepTransform3D transf = pv->getXToChildVol(ic);
 
    /*
    std::cout << " dumping transform to subcomponent" << std::endl;
    std::cout << transf[0][0]<<"," <<transf[0][1]<<"," <<transf[0][2]<<","<<transf[0][3] << std::endl;
    std::cout << transf[1][0]<<"," <<transf[1][1]<<"," <<transf[1][2]<<","<<transf[1][3] << std::endl;
    std::cout << transf[2][0]<<"," <<transf[2][1]<<"," <<transf[2][2]<<","<<transf[2][3] << std::endl;
    */
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    //std::cout << "  ";
    //std::cout << "subcomponent:"<<ic<<":"<<clv->getName()<<", made of"<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<< ","<< transf.getTranslation()<<std::endl;
	 
          if ( clv->getShape()->type()=="Trd") {
	      const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
	      /*
	      std::cout << "dimensions:"<< trd->getXHalfLength1() <<","
                                        << trd->getXHalfLength2() <<","  
                                        << trd->getYHalfLength1() <<","  
                                        << trd->getYHalfLength2() <<","  
			                << trd->getZHalfLength() <<std::endl; 
	      */
          } 
          if ( clv->getShape()->type()=="Box") {
	      const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
	      /*
	        std::cout << "dimensions:"<< box->getXHalfLength() <<","
                                        << box->getYHalfLength() <<","  
			                << box->getZHalfLength() <<std::endl; 
	      */
          } 

	  printChildren(cv);
    }  
   

}

const double Muon::MuonStationTypeBuilder::get_x_size(const GeoVPhysVol* pv) const
{
  double xlow = 0;
  double xup  = 0; 
  // subcomponents
  unsigned int nc = pv->getNChildVols();
  for (unsigned int ic=0; ic<nc; ic++) {
    HepTransform3D transf = pv->getXToChildVol(ic);
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    double xh=0;
    std::string type =  clv->getShape()->type();
    if (type=="Trd") {
       const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
       xh = fmax(trd->getXHalfLength1(),trd->getXHalfLength2()); 
    } 
    if (type=="Box") {
       const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
       xh = box->getXHalfLength();
    } 
    if (type=="Tube") {
       const GeoTube* tube = dynamic_cast<const GeoTube*> (clv->getShape());
       xh = tube->getRMax();
    } 
    if ( type!="Trd" && type!="Box" && type!="Tube") std::cout <<"unknown component type? "<< type << std::endl;
    xlow = fmin(xlow,(transf.getTranslation())[0]-xh);
    xup  = fmax(xup ,(transf.getTranslation())[0]+xh);
  }  
  // std::cout << "x size:" << xlow <<","<< xup << std::endl;
   
  return fmax( -xlow, xup ); 
}

const Trk::MaterialProperties* Muon::MuonStationTypeBuilder::getAveragedLayerMaterial( const GeoVPhysVol* pv, double volume, double thickness) const
{  
  // loop through the whole hierarchy; collect material
  Trk::MaterialProperties* empty=new Trk::MaterialProperties(0.,10.e8,0.,0.,0.);
  Trk::MaterialProperties total = collectMaterial( pv, *empty, volume/thickness);
  //std::cout << " combined material thickness: "<< total.thickness() << std::endl; 
  //std::cout << " actual layer thickness: "<< thickness << std::endl; 
  // scaled material properties to the actual layer thickness
  if (total.thickness() > 0 ) { 
      total *= thickness/total.thickness();
      return new Trk::MaterialProperties(total); 
  }  
  return 0;
}

Trk::MaterialProperties Muon::MuonStationTypeBuilder::collectMaterial(const GeoVPhysVol* pv, Trk::MaterialProperties matProp,double sf) const
{
  // sf is surface of the new layer used to calculate the average 'thickness' of components
  // number of child volumes
     unsigned int nc = pv->getNChildVols();
  // add current volume
    const GeoLogVol* lv = pv->getLogVol();
    //  std::cout << "component:"<<lv->getName()<<", made of"<<lv->getMaterial()->getName()<<","<<lv->getShape()->volume()<<std::endl;
     
    if ( lv->getMaterial()->getName() != "Air" && (lv->getName()).substr(0,1)!="T") {
            
      //std::cout << "including this in material calculations"<<std::endl;
      // get material properties from GeoModel
      Trk::MaterialProperties newMP= m_materialConverter->convert( lv->getMaterial() );
      // current volume
      double vol = lv->getShape()->volume();
      // subtract children volumes
      for (unsigned int ic=0; ic<nc; ic++) {
        const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
        vol = vol - cv->getLogVol()->getShape()->volume(); 
      }
      double d = vol / sf;
      //std::cout << "corrected volume thickness:" << d << std::endl;

      //  std::cout << "old material properties:" << matProp.thickness() <<","<<matProp.x0()<<","<<matProp.zOverAtimesRho()<<","<< 
      //	matProp.averageZ()<<","<<matProp.dEdX() << std::endl;
      // std::cout << "new material properties:" << d <<","<<newMP.x0()<<","<<newMP.zOverAtimesRho()<<","<< 
      //	newMP.averageZ()<<","<<newMP.dEdX() << std::endl;
      // material properties with thickness
      // that would be an easy way if TrkGeometry update
      // new.setThickness(d);
      // and this is as it works in the meantime
      Trk::MaterialProperties newUpdate(d,newMP.x0(),newMP.zOverAtimesRho(),newMP.averageZ(),newMP.dEdX());  
      // combine
      matProp.addMaterial(newUpdate);     
      //std::cout << "combined material properties:" << matProp.thickness() <<","<<matProp.x0()<<","<<matProp.zOverAtimesRho()<<","<< 
      //	matProp.averageZ()<<","<<matProp.dEdX() << std::endl;
   } 

  // subcomponents
  // skip children volume if we deal with G10 ( not correctly described )
  //if ( lv->getName() != "G10" ) { 
    for (unsigned int ic=0; ic<nc; ic++) {
      const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
      matProp = collectMaterial( cv, matProp, sf);
    }
    //}  
    //std::cout << " current layer thickness:" << matProp.thickness() << std::endl;    
  return matProp;
}

const Trk::LayerArray* Muon::MuonStationTypeBuilder::processCSCTrdComponent(const GeoVPhysVol*& pv, Trk::TrapezoidVolumeBounds*& compBounds, HepTransform3D*& transf) const {

  // tolerance
  double tol = 0.001;
  std::string name = pv->getLogVol()->getName();
  // std::cout << "processing CSC component, number of children volumes:"<< pv->getLogVol()->getName() << "," << pv->getNChildVols() <<std::endl; 
  // printChildren(pv);
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::MaterialProperties> x_mat;
  std::vector<double> x_thickness;
  double currX = -100000;
  // while waiting for better suggestion, define a single material layer
  Trk::MaterialProperties matCSC;
  double thickness =2* compBounds->halflengthY();
  double minX = compBounds->minHalflengthX();
  double maxX = compBounds->maxHalflengthX();
  double halfY = compBounds->halflengthY();
  if (name.substr( name.size()-5,5 ) == "CSC01" ) {
    if (!m_matCSC01 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matCSC01 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matCSC = *m_matCSC01; 
  } 
  if (name == "CSCspacer" ) {
    if (!m_matCSCspacer1 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matCSCspacer1 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matCSC = *m_matCSCspacer1; 
  } 
  x_array.push_back(0.);
  x_mat.push_back(matCSC);
  x_thickness.push_back(thickness );
  // create layers
  const Trk::PlaneLayer* layer;
  Trk::OverlapDescriptor* od=0;
  for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
    Trk::TrapezoidBounds* bounds= new Trk::TrapezoidBounds(minX,maxX,halfY); ;
    HepTransform3D* cTr = new HepTransform3D( (*transf) * HepTranslateZ3D(x_array[iloop]) ); // this won't work for multiple layers !!! //
    Trk::HomogenousLayerMaterial cscMaterial(x_mat[iloop], Trk::oppositePre);  
    layer = new Trk::PlaneLayer(cTr,
                                bounds,
                                cscMaterial,
                                x_thickness[iloop],
                                od );
    layers.push_back(layer);
    //std::cout << "CSC layer built ok"<<std::endl;
  }
  // create the BinnedArray
  std::cout << "number of Csc layers:"<<layers.size()<<std::endl;
  std::vector<LayTr> layerOrder;
  std::vector<double> binSteps;
  double lowX = - compBounds->halflengthZ() ;

  if (layers.size()) {
     lowX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = lowX; 
     for (unsigned int i=0;i<layers.size();i++) { 
       const HepTransform3D* ltransf = new HepTransform3D(layers[i]->transform());
       layerOrder.push_back(LayTr(Trk::SharedObject<const Trk::Layer>(layers[i]), ltransf ));
       binSteps.push_back(ltransf->getTranslation()[0]+0.5*layers[i]->thickness()-currX);
       currX = ltransf->getTranslation()[0]+0.5*layers[i]->thickness();     
     }
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( lowX, new std::vector<double>(binSteps));
  Trk::LayerArray* cscLayerArray = 0;
  cscLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     

  return cscLayerArray;

} 

const Trk::LayerArray* Muon::MuonStationTypeBuilder::processCSCDiamondComponent(const GeoVPhysVol*& pv, Trk::DoubleTrapezoidVolumeBounds*& compBounds, HepTransform3D*& transf) const {

  // tolerance
  double tol = 0.001;
  std::string name = pv->getLogVol()->getName();
  //std::cout << "processing CSC component, number of children volumes:"<< pv->getLogVol()->getName() << "," << pv->getNChildVols() <<std::endl; 
  //printChildren(pv);
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::MaterialProperties> x_mat;
  std::vector<double> x_thickness;
  double currX = -100000;
  // while waiting for better suggestion, define a single material layer
  Trk::MaterialProperties matCSC;
  double thickness = 2*compBounds->halflengthZ();
  double minX = compBounds->minHalflengthX();
  double medX = compBounds->medHalflengthX();
  double maxX = compBounds->maxHalflengthX();
  double halfY1 = compBounds->halflengthY1();
  double halfY2 = compBounds->halflengthY2();
  double halfZ = compBounds->halflengthZ();
  if (name.substr( name.size()-5,5 ) == "CSC02" ) {
    if (!m_matCSC02 ) { 
      double vol = ( (minX + medX)*2*halfY1+(medX+maxX)*2*halfY2 ) * thickness;
      m_matCSC02 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matCSC = *m_matCSC02; 
  } 
  if (name == "CSCspacer" ) {
    if (!m_matCSCspacer2 ) { 
      double vol = ( (minX + medX)*2*halfY1+(medX+maxX)*2*halfY2 ) * thickness;
      m_matCSCspacer2 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matCSC = *m_matCSCspacer2; 
  } 
  x_array.push_back(0.);
  x_mat.push_back(matCSC);
  x_thickness.push_back(thickness );
  // create layers
  const Trk::PlaneLayer* layer;
  Trk::OverlapDescriptor* od=0;
  for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
    Trk::DiamondBounds* bounds= new Trk::DiamondBounds(minX,medX,maxX,halfY1,halfY2); ;
    HepTransform3D* cTr = new HepTransform3D( (*transf)* HepRotateX3D(90*deg)  * HepTranslateZ3D(x_array[iloop]) ); // this won't work for multiple layers !!! //
    Trk::HomogenousLayerMaterial cscMaterial(x_mat[iloop], Trk::oppositePre);  
    layer = new Trk::PlaneLayer(cTr,
                                bounds,
                                cscMaterial,
                                x_thickness[iloop],
                                od );
    layers.push_back(layer);
    //std::cout << "CSC layer built ok"<<std::endl;
  }
  // create the BinnedArray
  std::cout << "number of Csc layers:"<<layers.size()<<std::endl;
  std::vector<LayTr> layerOrder;
  std::vector<double> binSteps;
  double lowX = - compBounds->halflengthZ() ;
  currX = lowX;
   
  for (unsigned int i=0;i<layers.size();i++) { 
    const HepTransform3D* ltransf = new HepTransform3D(layers[i]->transform());
    layerOrder.push_back(LayTr(Trk::SharedObject<const Trk::Layer>(layers[i]), ltransf ));
    if (i>0) {
      binSteps.push_back(ltransf->getTranslation()[0] -currX);
    } else {
      binSteps.push_back(2*(ltransf->getTranslation()[0] -currX));
    } 
    currX = ltransf->getTranslation()[0];
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( lowX, new std::vector<double>(binSteps));
  Trk::LayerArray* cscLayerArray = 0;
  cscLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     

  return cscLayerArray;

} 

const Trk::LayerArray* Muon::MuonStationTypeBuilder::processTGCComponent(const GeoVPhysVol*& pv, Trk::TrapezoidVolumeBounds*& tgcBounds, HepTransform3D*& transf) const {

  // tolerance
  double tol = 0.001;
  std::string name = pv->getLogVol()->getName();
  // std::cout << "processing TGC component, number of children volumes:"<< pv->getLogVol()->getName() << "," << pv->getNChildVols() <<std::endl; 
  // printChildren(pv);
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::MaterialProperties> x_mat;
  std::vector<double> x_thickness;
  double currX = -100000;
  // while waiting for better suggestion, define a single material layer
  Trk::MaterialProperties matTGC;
  double thickness =2*tgcBounds->halflengthZ();
  double minX = tgcBounds->minHalflengthX();
  double maxX = tgcBounds->maxHalflengthX();
  double halfY = tgcBounds->halflengthY();
  //std::cout << "tgc bounds half y:" << tgcBounds->halflengthY() << std::endl; 
  if ( tgcBounds->halflengthZ() == 35.00 ) {
    if (!m_matTGC01 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matTGC01 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matTGC = *m_matTGC01; 
  } 
  if ( tgcBounds->halflengthZ() == 21.85 ) {
    if (!m_matTGC06 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matTGC06 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matTGC = *m_matTGC06; 
  } 
  x_array.push_back(0.);
  x_mat.push_back(matTGC);
  x_thickness.push_back(thickness );
  // create layers
  const Trk::PlaneLayer* layer;
  Trk::OverlapDescriptor* od=0;
  for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
    Trk::TrapezoidBounds* bounds= new Trk::TrapezoidBounds(minX,maxX,halfY); ;
    HepTransform3D* cTr = new HepTransform3D( (*transf) * HepTranslateZ3D(x_array[iloop]) ); // this won't work for multiple layers !!! //
    Trk::HomogenousLayerMaterial tgcMaterial(x_mat[iloop], Trk::oppositePre);  
    layer = new Trk::PlaneLayer(cTr,
                                bounds,
                                tgcMaterial,
                                x_thickness[iloop],
                                od );
    layers.push_back(layer);
    //std::cout << "TGC layer built ok"<<std::endl;
  }
  // create the BinnedArray
  //std::cout << "number of Tgc layers:"<<layers.size()<<std::endl;
  std::vector<LayTr> layerOrder;
  std::vector<double> binSteps;
  // 
  double lowX = - tgcBounds->halflengthZ();
  if (layers.size()) {
     lowX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = lowX; 
     for (unsigned int i=0;i<layers.size();i++) { 
       const HepTransform3D* ltransf = new HepTransform3D(layers[i]->transform());
       layerOrder.push_back(LayTr(Trk::SharedObject<const Trk::Layer>(layers[i]), ltransf ));
       binSteps.push_back(ltransf->getTranslation()[0]+0.5*layers[i]->thickness()-currX);
       currX = ltransf->getTranslation()[0]+0.5*layers[i]->thickness();     
     }
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( lowX, new std::vector<double>(binSteps));
  Trk::LayerArray* tgcLayerArray = 0;
  tgcLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     

  return tgcLayerArray;

} 
