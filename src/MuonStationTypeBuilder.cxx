///////////////////////////////////////////////////////////////////
// MuonStationTypeBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonStationTypeBuilder.h"
//MuonSpectrometer include
#include "MuonGeoModel/MuonDetectorManager.h"
#include "MuonGeoModel/MuonStation.h"
#include "MuonGeoModel/MdtReadoutElement.h"
#include "MuonIdHelpers/MuonIdHelper.h"
// Trk
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkDetDescrUtils/BinUtility1DX.h"
#include "TrkDetDescrUtils/BinUtility1DY.h"
#include "TrkDetDescrUtils/BinUtility1DZ.h"
#include "TrkDetDescrUtils/BinnedArray.h"
#include "TrkDetDescrUtils/NavBinnedArray1D.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/DoubleTrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkVolumes/CombinedVolumeBounds.h"
#include "TrkVolumes/VolumeExcluder.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkSurfaces/TrapezoidBounds.h"
#include "TrkSurfaces/DiamondBounds.h"
#include "TrkGeometrySurfaces/SubtractedPlaneSurface.h"
#include "TrkGeometry/CylinderLayer.h"
#include "TrkGeometry/DiscLayer.h"
#include "TrkGeometry/PlaneLayer.h"
#include "TrkGeometry/SubtractedPlaneLayer.h"
#include "TrkGeometry/ExtendedMaterialProperties.h"
#include "TrkGeometry/LayerMaterialProperties.h"
#include "TrkGeometry/HomogenousLayerMaterial.h"
#include "TrkGeometry/OverlapDescriptor.h"
#include "TrkMagFieldInterfaces/IMagneticFieldTool.h"
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

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Geometry/Transform3D.h"

// STD
#include <map>

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
  m_trackingVolumeArrayCreator("Trk::TrackingVolumeArrayCreator/TrackingVolumeArrayCreator"),
  m_magFieldTool("Trk::MagneticFieldTool/AtlasMagneticFieldTool"),
  m_mdtTubeMat(0),
  m_mdtFoamMat(),
  m_rpc46(0),
  m_rpcDed(),
  m_rpcLayer(0),
  m_rpcExtPanel(0),
  m_rpcMidPanel(0),
  m_matCSC01(0),
  m_matCSCspacer1(0),
  m_matCSC02(0),
  m_matCSCspacer2(0),
  m_matTGC01(0),
  m_matTGC06(0),
  m_log(0)
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
    
    m_log = new MsgStream(msgSvc(), name());

    StatusCode s = AlgTool::initialize();
    if (s.isFailure()) *m_log<< MSG::INFO << "failing to initialize?" << endreq;

    // Get DetectorStore service
    //
    StoreGateSvc* m_detStore=0;
    StatusCode ds = service("DetectorStore",m_detStore);
    if (ds.isFailure()) {
        *m_log << MSG::FATAL << "DetectorStore service not found !" << endreq;
    }
    // get Muon Spectrometer Description Manager
    // s = m_detStore->retrieve(m_muonMgr, m_muonMgrLocation);

    ds = m_detStore->retrieve(m_muonMgr);

    if (ds.isFailure()) {
        *m_log << MSG::ERROR << "Could not get MuonDetectorManager, no layers for muons will be built. " << endreq;
    }

    *m_log << MSG::DEBUG << m_muonMgr->geometryVersion() << endreq; 
    
    // Retrieve the magnetic field tool   ----------------------------------------------------    
    if (m_magFieldTool.retrieve().isFailure())
    {
      *m_log << MSG::FATAL << "Failed to retrieve tool " << m_magFieldTool << endreq;
      return StatusCode::FAILURE;
    } else
      *m_log << MSG::INFO << "Retrieved tool " << m_magFieldTool << endreq;

    // Retrieve the tracking volume array creator   -------------------------------------------    
    if (m_trackingVolumeArrayCreator.retrieve().isFailure()) {
      *m_log << MSG::FATAL << "Failed to retrieve tool " << m_trackingVolumeArrayCreator << endreq;
      return StatusCode::FAILURE;
    } else
      *m_log << MSG::INFO << "Retrieved tool " << m_trackingVolumeArrayCreator << endreq;

    // default (trivial) muon material properties 
    m_muonMaterial = new Trk::MaterialProperties(0.,10e8,0.);

    Trk::MagneticFieldProperties muonMagneticFieldProperties(&(*m_magFieldTool), Trk::RealisticField);    
    m_muonMagneticField = muonMagneticFieldProperties;

    m_materialConverter= new Trk::GeoMaterialConverter();
     
    *m_log << MSG::INFO  << name() <<" initialize() successful" << endreq;    
    
    return StatusCode::SUCCESS;
}

const Trk::TrackingVolumeArray* Muon::MuonStationTypeBuilder::processBoxStationComponents(const GeoVPhysVol* mv, Trk::CuboidVolumeBounds* envelope) const
{

    *m_log << MSG::DEBUG  << name() <<" processing station components for " <<mv->getLogVol()->getName()<< endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////
   
   double tolerance = 0.0001;   

   // loop over children volumes; ( make sure they do not exceed enveloping volume boundaries ?)
   // split into connected subvolumes ( assume ordering along X unless otherwise )
      std::vector<Trk::Volume*> compVol;
      std::vector<std::string> compName;
      std::vector<const GeoVPhysVol*> compGeo;
      std::vector<HepTransform3D> compTransf;
      double halfZ=0.;   
      double halfX1=0.;
      double halfX2=0.;
      double halfY1=0.;
      double halfY2=0.;
      for (unsigned int ich =0; ich< mv->getNChildVols(); ++ich) 
      {
        const GeoVPhysVol* cv = &(*(mv->getChildVol(ich))); 
        const GeoLogVol* clv = cv->getLogVol();
        HepTransform3D transf = mv->getXToChildVol(ich);        
      // TEMPORARY CORRECTION 
        if ( (mv->getLogVol()->getName()).substr(0,3)=="BMF" && (clv->getName()).substr(0,2)=="LB" ) {
	  *m_log<< MSG::DEBUG << "TEMPORARY MANUAL CORRECTION OF BMF SPACER LONG BEAM POSITION" << endreq;
            transf = transf * HepTranslate3D(-37.5,0.,0.);
        } 
        // retrieve volumes for components
	Trk::VolumeBounds* volBounds=0; 
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
 	  //std::cout << "WARNING:component shape not Box nor Trapezoid, determining the x size from subcomponents" << std::endl; 
          double xSize = get_x_size(cv);
	  //std::cout << "estimated x size:" << xSize << std::endl;
          //printChildren(cv);
          volBounds = new Trk::CuboidVolumeBounds(xSize,envelope->halflengthY(),envelope->halflengthZ());
        }
	vol = new Trk::Volume(new HepTransform3D(transf),volBounds);
	//std::cout <<"volume center:"<< vol->center() << ","<< ich << std::endl;
	std::string cname = clv->getName();
	std::string vname = mv->getLogVol()->getName();
        int nameSize = vname.size()-8;
        if (cname.substr(0,nameSize) == vname.substr(0,nameSize)) cname = cname.substr(nameSize,cname.size()-nameSize);  
        // order in X
        if (compVol.size()==0 || vol->center()[0]>=compVol.back()->center()[0]){
          compVol.push_back(vol);
          compName.push_back(cname);
          compGeo.push_back(cv);
          compTransf.push_back(transf);
        } else {
	  std::vector<Trk::Volume*>::iterator volIter=compVol.begin();
	  std::vector<std::string>::iterator  nameIter=compName.begin();
	  std::vector<const GeoVPhysVol*>::iterator  geoIter=compGeo.begin();
	  std::vector<HepTransform3D>::iterator  transfIter=compTransf.begin();
          while ( vol->center()[0]>= (*volIter)->center()[0]) {volIter++;nameIter++;geoIter++;transfIter++;}
          compVol.insert(volIter,vol);
          compName.insert(nameIter,cname);
          compGeo.insert(geoIter,cv);
          compTransf.insert(transfIter,transf);
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
      std::vector<HepTransform3D> transfSpacer;
      std::vector<HepTransform3D> transfRpc;
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
        // close Rpc if no further components
        if (openRpc  && compName[i].substr(0,3) != "RPC" && compName[i].substr(0,3) != "Ded"){
          // low edge of current volume
          double Xcurr = compVol[i]->center()[0]-compBounds->halflengthX();
          if (Xcurr >= currX+rpclowXsize+rpcuppXsize) {
            Trk::CuboidVolumeBounds* rpcBounds = new Trk::CuboidVolumeBounds(0.5*(Xcurr-currX),envY,envZ);
	    Trk::Volume* rpcVol =new Trk::Volume(new HepTranslate3D(currX+rpcBounds->halflengthX(),0.,0.),rpcBounds);
            const Trk::TrackingVolume* rpcTrkVol = processRpc(rpcVol ,geoRpc,transfRpc);
            delete rpcVol;
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
            const Trk::TrackingVolume* spacerTrkVol = processSpacer(spacerVol ,geoSpacer, transfSpacer);
            delete spacerVol;
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
	  Trk::CuboidVolumeBounds* mdtBounds=0; 
          // remove z shift in transform !! bugfix !!
          double zShift = compVol[i]->transform().getTranslation()[2];
	  if ( fabs(zShift)>0 ) {
	    *m_log << MSG::DEBUG << "unusual z shift for subvolume:" << zShift << endreq; 
          }  
          //                                 (HepTranslateZ3D(-zShift)*(*compTransf[i])).getTranslation() <<std::endl;
          if (lowX == currX) {
	    mdtBounds = new Trk::CuboidVolumeBounds(compBounds->halflengthX(),envY,envZ);
            mdtVol = new Trk::Volume(new HepTransform3D(HepTranslateZ3D(-zShift)*compVol[i]->transform()),mdtBounds);
	  } else {
            if ( fabs(lowX-currX)>0.002 ) { 
	      *m_log << MSG::WARNING  << "Mdt volume size does not match the envelope:lowX,currX:" << lowX <<","<<currX << std::endl;
	      *m_log << MSG::WARNING  << "adjusting Mdt volume " << std::endl;
	    }
	    mdtBounds = new Trk::CuboidVolumeBounds(compBounds->halflengthX()+0.5*(lowX-currX),envY,envZ);
            mdtVol = new Trk::Volume(new HepTransform3D(HepTranslateX3D(0.5*(currX-lowX))
							*HepTranslateZ3D(-zShift)*compVol[i]->transform()),mdtBounds);
          }
          double shiftSign = 1.; 
          if (fabs(zShift) > 0.) {
	    std::string stName = mv->getLogVol()->getName();
	    if ( stName.substr(0,4)=="BIR3" || stName.substr(0,4)=="BIR5" || stName.substr(0,4)=="BIR7" || stName.substr(0,5)=="BIR10" ) shiftSign = -1.;
          }
          const Trk::TrackingVolume* mdtTrkVol = processMdtBox(mdtVol,compGeo[i],
					  new HepTransform3D(HepTranslateZ3D(-zShift)*compTransf[i]),
							       shiftSign*fabs(zShift));
          trkVols.push_back(mdtTrkVol); 
          volSteps->push_back(2.*mdtBounds->halflengthX()); 
          currX += 2.*mdtBounds->halflengthX();
          delete mdtVol;
          comp_processed = true;
          zShift = 0.;
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
          delete spacerVol;
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
          delete rpcVol;
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
      Trk::BinUtility* binUtility = new Trk::BinUtility1DX( -( envelope->halflengthX() ), volSteps);
      if (m_trackingVolumeArrayCreator)  components = m_trackingVolumeArrayCreator->cuboidVolumesArrayNav( trkVols, binUtility, false);
      // std::cout << "tracking volume array created" << std::endl;

      for (size_t i=0; i<compVol.size(); i++) delete compVol[i];

     
   return components;  
}

const Trk::TrackingVolumeArray* Muon::MuonStationTypeBuilder::processTrdStationComponents(const GeoVPhysVol* mv, Trk::TrapezoidVolumeBounds* envelope ) const
{
    *m_log << MSG::DEBUG  << name() <<" processing station components for " <<mv->getLogVol()->getName() << endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////

    double tolerance = 0.0001;   
   
   // loop over children volumes; ( make sure they do not exceed enveloping volume boundaries ?)
   // split into connected subvolumes ( assume ordering along X unless otherwise )
      std::vector<Trk::Volume*> compVol;
      std::vector<std::string> compName;
      std::vector<const GeoVPhysVol*> compGeo;
      std::vector<HepTransform3D> compTransf;
      double halfZ=0.;   
      double halfX1=0.;
      double halfX2=0.;
      double halfY1=0.;
      double halfY2=0.;
      for (unsigned int ich =0; ich< mv->getNChildVols(); ++ich) 
      {
	// std::cout << "next component:"<< ich << std::endl;
        const GeoVPhysVol* cv = &(*(mv->getChildVol(ich))); 
        const GeoLogVol* clv = cv->getLogVol();
        HepTransform3D transf = mv->getXToChildVol(ich);        
        // std::cout << "component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","  <<transf.getTranslation()<<std::endl;
        // retrieve volumes for components
	Trk::VolumeBounds* volBounds=0; 
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
 	  //std::cout<<"WARNING:component shape not Box nor Trapezoid, determining the x size from subcomponents"<<std::endl; 
	  double xSize = get_x_size(cv);
          // printChildren(cv);
	  if (clv->getName().substr(0,1)!="C" && clv->getName().substr(0,2)!="LB")
	    transf = transf*HepRotateY3D(90*deg)*HepRotateZ3D(90*deg);
          volBounds = new Trk::TrapezoidVolumeBounds(envelope->minHalflengthX(),envelope->maxHalflengthX(),envelope->halflengthY(),xSize);
        }
	vol = new Trk::Volume(new HepTransform3D(transf),volBounds);
	// std::cout <<"volume center:"<< vol->center() << ","<< ich << std::endl;
	std::string cname = clv->getName();
	std::string vname = mv->getLogVol()->getName();
        int nameSize = vname.size()-8;
        if (cname.substr(0,nameSize) == vname.substr(0,nameSize)) cname = cname.substr(nameSize,cname.size()-nameSize);  
        // order in X
        if (compVol.size()==0 || vol->center()[0]>=compVol.back()->center()[0]){
          compVol.push_back(vol);
          compName.push_back(cname);
          compGeo.push_back(cv);
          compTransf.push_back(transf);
        } else {
	  std::vector<Trk::Volume*>::iterator volIter=compVol.begin();
	  std::vector<std::string>::iterator  nameIter=compName.begin();
	  std::vector<const GeoVPhysVol*>::iterator  geoIter=compGeo.begin();
	  std::vector<HepTransform3D>::iterator  transfIter=compTransf.begin();
          while ( vol->center()[0]>= (*volIter)->center()[0]) {volIter++;nameIter++;geoIter++;transfIter++;}
          compVol.insert(volIter,vol);
          compName.insert(nameIter,cname);
          compGeo.insert(geoIter,cv);
          compTransf.insert(transfIter,transf);
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
      std::vector<HepTransform3D> transfSpacer;
      std::vector<HepTransform3D> transfRpc;
      double spacerlowXsize=0; 
      double spaceruppXsize=0; 
      double rpclowXsize=0; 
      double rpcuppXsize=0; 
      double Xcurr=0;
      double lowX=0.;
      double uppX=0.;
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
        }   
        // close spacer if no further components
        if (openSpacer &&  compName[i].substr(0,1) != "C" && compName[i].substr(0,2) != "LB"){
          if (Xcurr-currX-(spacerlowXsize+spaceruppXsize)>= -tolerance ) {
            Trk::TrapezoidVolumeBounds* spacerBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,0.5*(Xcurr-currX));
            //HepTransform3D tr = HepRotateZ3D(90*deg)*HepRotateX3D(90*deg)*HepTranslate3D(currX+spacerBounds->halflengthZ(),0.,0.);
            HepTransform3D tr = HepTranslate3D(currX+spacerBounds->halflengthZ(),0.,0.)* HepRotateY3D(90*deg)*HepRotateZ3D(90*deg);
	    Trk::Volume* spacerVol =new Trk::Volume(new HepTransform3D(tr),spacerBounds);
            const Trk::TrackingVolume* spacerTrkVol = processSpacer(spacerVol ,geoSpacer, transfSpacer);
            delete spacerVol;   
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
          Trk::Volume* mdtVol=0;
	  Trk::TrapezoidVolumeBounds* mdtBounds=0; 
	  //std::cout <<"processMdtTrd:"<<envX1 <<","<<envX2 <<","<<envY<<","<<compTrdBounds->halflengthZ()<<std::endl;
          if (lowX == currX) {
	    mdtBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,compTrdBounds->halflengthZ());
            mdtVol = new Trk::Volume(new HepTransform3D(compVol[i]->transform()),mdtBounds);
	  } else {
            if (fabs(lowX-currX)>0.002 ) {
	      *m_log << MSG::WARNING  << "Mdt volume size does not match the envelope:lowX,currX:" << lowX <<","<<currX << std::endl;
	      *m_log << MSG::WARNING  << "adjusting Mdt volume " << std::endl;
	    }
	    mdtBounds = new Trk::TrapezoidVolumeBounds(envX1,envX2,envY,compTrdBounds->halflengthZ()+0.5*(lowX-currX));
            mdtVol = new Trk::Volume(new HepTransform3D(HepTranslateZ3D(0.5*(currX-lowX))*compVol[i]->transform()),mdtBounds);
          }
          const Trk::TrackingVolume* mdtTrkVol = processMdtTrd(mdtVol,compGeo[i],&compTransf[i]);
          trkVols.push_back(mdtTrkVol); 
          volSteps->push_back(2*mdtBounds->halflengthZ());
          currX += 2.*mdtBounds->halflengthZ();
          delete mdtVol;
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
      Trk::BinUtility* binUtility = new Trk::BinUtility1DX( -( envelope->halflengthZ() ), volSteps);
      if (m_trackingVolumeArrayCreator)  components = m_trackingVolumeArrayCreator->trapezoidVolumesArrayNav( trkVols, binUtility, false);
      // std::cout << "tracking volume array created" << std::endl;

      for (size_t i=0; i<compVol.size(); i++) delete compVol[i];
     
   return components;  
}

// finalize
StatusCode Muon::MuonStationTypeBuilder::finalize()
{
    
    delete m_matCSCspacer1;
    delete m_matCSCspacer2;
    delete m_materialConverter;
    delete m_matCSC01;
    delete m_matCSC02;
    delete m_matTGC01;
    delete m_matTGC06;
    delete m_mdtTubeMat;
    delete m_rpcLayer;
    delete m_rpcMidPanel;
    delete m_rpcExtPanel;
    delete m_muonMaterial;
    for (unsigned int i=0;i<m_rpcDed.size();i++) delete m_rpcDed[i];
    for (unsigned int i=0;i<m_mdtFoamMat.size();i++) delete m_mdtFoamMat[i];

    *m_log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    delete m_log; 
    m_log = 0;
 
    return StatusCode::SUCCESS;
}
//
const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processMdtBox(Trk::Volume*& vol,const GeoVPhysVol*& gv, HepTransform3D* transf, double zShift) const
{
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<double> x_ref;
  std::vector<Trk::ExtendedMaterialProperties> x_mat;
  std::vector<double> x_thickness;
  std::vector<int> x_active;
  double currX = -100000; 
  // here one could save time by not reading all tubes  
  for (unsigned int ich =0; ich< gv->getNChildVols(); ++ich) {
    const GeoVPhysVol* cv = &(*(gv->getChildVol(ich))); 
    const GeoLogVol* clv = cv->getLogVol();
    HepTransform3D transfc = gv->getXToChildVol(ich);        
    //std::cout << "MDT component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transfc.getTranslation()<<std::endl;
    // printChildren(cv);
    Trk::ExtendedMaterialProperties mdtMat(0.,10e8,10e8,0.,0.,0.);
    double xv = 0.;
    int active = 0;
    if ( (clv->getName()).substr(0,3)=="MDT") {
      xv = 13.0055;  // the half-thickness
      if (!m_mdtTubeMat ) {
	const GeoTube* tube=dynamic_cast<const GeoTube*> (clv->getShape());
	double volume = 8*(tube->getRMax())*(tube->getZHalfLength())*xv;
	//std::cout << " part of layer volume assigned to 1 tube:" << volume << std::endl;
	//std::cout << "tube dimensions:" << tube->getRMax() << "," << tube->getRMin() << "," << tube->getZHalfLength() << std::endl;
	m_mdtTubeMat = getAveragedLayerMaterial(cv,volume,2*xv); 
      }        
      mdtMat = *m_mdtTubeMat;
      active = 1; 
    }
    if ( (clv->getName())=="MultiLayerFoam") {
      //std::cout << "processing MultiLayerFoam" << std::endl;  
      xv = decodeX(clv->getShape());
      for (unsigned int i=0;i<m_mdtFoamMat.size();i++) {
        if ( fabs(xv-0.5*m_mdtFoamMat[i]->thickness())<0.001 ) {
	  mdtMat = *(m_mdtFoamMat[i]);
          break;
	}
      }
      if (mdtMat.thickness()<0.001) {
        //std::cout << "adding new box MDT:foam thickness:"<< xv << std::endl;
        const Trk::CuboidVolumeBounds* cub =  dynamic_cast<const Trk::CuboidVolumeBounds*> (&(vol->volumeBounds()));
	double volume = 8*(cub->halflengthY())*(cub->halflengthZ())*xv;
	m_mdtFoamMat.push_back(getAveragedLayerMaterial(cv,volume,2*xv)); 
        mdtMat = *(m_mdtFoamMat.back());
      }
    }
    if (  transfc.getTranslation()[0] != currX ) {
      if (x_array.size() == 0 || transfc.getTranslation()[0] > x_array.back() ) {
	x_array.push_back(  transfc.getTranslation()[0] );
	x_mat.push_back(mdtMat);
	x_thickness.push_back( 2*xv );
        x_active.push_back(active);
	currX = transfc.getTranslation()[0];
	if ( fabs(transfc.getTranslation()[1])>0.001) {
	  // code 2.corrdinate shift
	  double ref =  transfc.getTranslation()[2]+1e5;
	  ref += int(1000*transfc.getTranslation()[1])*10e6; 
	  x_ref.push_back( ref ) ;
	} else {
	  x_ref.push_back( transfc.getTranslation()[2] ) ;
	}
	// std::cout << "layer info included:" << clv->getName()<<"," << 2*xv <<","<< currX<< std::endl; 
      } else {
	std::vector<double>::iterator xIter=x_array.begin();
	std::vector<Trk::ExtendedMaterialProperties>::iterator mIter=x_mat.begin();
	std::vector<double>::iterator tIter=x_thickness.begin();
	std::vector<double>::iterator rIter=x_ref.begin();
	std::vector<int>::iterator aIter=x_active.begin();
	while ( transfc.getTranslation()[0] > *xIter ) {xIter++;mIter++;rIter++;}
	x_array.insert(xIter,transfc.getTranslation()[0]);
	x_mat.insert(mIter,mdtMat);
	x_thickness.insert(tIter,2*xv);
        x_active.insert(aIter,active);
	if ( fabs(transfc.getTranslation()[1])>0.001) {
	  // code 2.corrdinate shift
	  double sign = (transfc.getTranslation()[1]>0.) ? 1. : -1.;
	  double ref =  transfc.getTranslation()[2]+sign*1e5;            
	  ref += int(1000*transfc.getTranslation()[1])*10e6; 
	  x_ref.insert( rIter,ref ) ;
	} else {
	  x_ref.insert(rIter,transfc.getTranslation()[2] ) ;
	}
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
    Trk::RectangleBounds* rbounds = new Trk::RectangleBounds(yv,zv); 
    Trk::SharedObject<const Trk::SurfaceBounds> bounds(rbounds);
    for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
      // x-y plane -> y-z plane
      thickness = x_thickness[iloop];
      HepTransform3D* cTr = new HepTransform3D( (*transf) * HepTranslateX3D(x_array[iloop])
						* HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      Trk::HomogenousLayerMaterial mdtMaterial(x_mat[iloop]);  
      layer = new Trk::PlaneLayer(cTr,
                                  bounds,
                                  mdtMaterial,
                                  thickness,
                                  od );
      layer->setRef(x_ref[iloop]-zShift);
      //std::cout << " reference value set for layer:"<<iloop<<","<<layer->getRef()<<std::endl;
      //std::cout << transf->getTranslation() << std::endl;
      layers.push_back(layer);
      
      // make preliminary identification of active layers
      layer->setLayerType(x_active[iloop]);
    }
  } 
  // create the BinnedArray
  //std::cout << "number of Mdt layers:"<<layers.size()<<std::endl;
  std::vector<Trk::SharedObject<const Trk::Layer> > layerOrder;
  std::vector<double> binSteps;
  // check if additional (navigation) layers needed
  // fix lower and upper bound of step vector to volume boundary
  // double minX = - volBounds->halflengthX();
  double minX = transf->getTranslation()[0] - volBounds->halflengthX();
  if (layers.size()) {
    //minX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
    currX = minX; 
    for (unsigned int i=0;i<layers.size();i++) { 
      const HepTransform3D ltransf = layers[i]->transform();
      layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers[i]));
      if (i<layers.size()-1) { 
	binSteps.push_back(ltransf.getTranslation()[0]+0.5*layers[i]->thickness()-currX);
	currX = ltransf.getTranslation()[0]+0.5*layers[i]->thickness();     
      }
    }
    binSteps.push_back(transf->getTranslation()[0]+volBounds->halflengthX()-currX);
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( minX, new std::vector<double>(binSteps));
  Trk::LayerArray* mdtLayerArray = 0;
  mdtLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     
  std::string name="MDT";
  const Trk::TrackingVolume* mdt= new Trk::TrackingVolume(*vol,
                                                          *m_muonMaterial,
                                                          m_muonMagneticField,
                                                          mdtLayerArray,0,
                                                          name);         
  delete transf;
  return mdt;
}
//
const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processMdtTrd(Trk::Volume*& vol,const GeoVPhysVol*& gv, HepTransform3D* transf) const
{
   // std::cout << "processing MDT, number of children volumes:"<< gv->getNChildVols() <<std::endl; 
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::ExtendedMaterialProperties> x_mat;
  std::vector<double> x_thickness;
  std::vector<double> x_ref;
  std::vector<int> x_active;
  double currX = -100000; 
  for (unsigned int ich =0; ich< gv->getNChildVols(); ++ich) {
    const GeoVPhysVol* cv = &(*(gv->getChildVol(ich))); 
    const GeoLogVol* clv = cv->getLogVol();
    HepTransform3D transfc = gv->getXToChildVol(ich);        
    //std::cout << "MDT component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transfc.getTranslation()<<std::endl;
    double xv = 0.;
    int active = 0;
    if (clv->getShape()->type()=="Trd"){
      const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
      double x1v = trd->getXHalfLength1();
      double x2v = trd->getXHalfLength2();
      if ( x1v==x2v ) xv = x1v;
    }
    Trk::ExtendedMaterialProperties mdtMat(0.,10e8,10e8,0.,0.,0.); 
    if ( (clv->getName()).substr(0,3)=="MDT") {
      xv = 13.0055;  // the half-thickness
      if (!m_mdtTubeMat ) {
	const GeoTube* tube=dynamic_cast<const GeoTube*> (clv->getShape());
	double volume = 8*(tube->getRMax())*(tube->getZHalfLength())*xv;
	// std::cout << " part of layer volume assigned to 1 tube:" << vol << std::endl;
	// std::cout << "tube dimensions:" << tube->getRMax() << "," << tube->getRMin() << "," << tube->getZHalfLength() << std::endl;
	m_mdtTubeMat = getAveragedLayerMaterial(cv,volume,2*xv); 
      }        
      mdtMat = *m_mdtTubeMat;
      active = 1;
    }
    if ( (clv->getName())=="MultiLayerFoam") {
      //std::cout << "processing MultiLayerFoam" << std::endl;  
      xv = decodeX(clv->getShape());
      for (unsigned int i=0;i<m_mdtFoamMat.size();i++) {
        if ( fabs(xv-0.5*m_mdtFoamMat[i]->thickness())<0.001 ) {
	  mdtMat = *(m_mdtFoamMat[i]);
          break;
	}
      }
      if (mdtMat.thickness()<0.001) {
        //std::cout << "adding new trd MDT:foam thickness:"<< xv << std::endl;
        const Trk::TrapezoidVolumeBounds* trd =  dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(vol->volumeBounds()));
	double volume = 4*(trd->minHalflengthX()+trd->maxHalflengthX())*(trd->halflengthY())*xv;
	m_mdtFoamMat.push_back(getAveragedLayerMaterial(cv,volume,2*xv)); 
        mdtMat = *(m_mdtFoamMat.back());
      }
    }
    
    if (  transfc.getTranslation()[0] != currX ) {
      if (x_array.size() == 0 || transfc.getTranslation()[0] > x_array.back() ) {
	x_array.push_back(  transfc.getTranslation()[0] );
	x_mat.push_back( mdtMat );
	x_thickness.push_back( 2*xv );
	x_ref.push_back(  transfc.getTranslation()[2] );
	currX = transfc.getTranslation()[0];
	x_active.push_back(active);
      } else {
	std::vector<double>::iterator xIter=x_array.begin();
	std::vector<Trk::ExtendedMaterialProperties>::iterator mIter=x_mat.begin();
	std::vector<double>::iterator tIter=x_thickness.begin();
	std::vector<double>::iterator rIter=x_ref.begin();
	std::vector<int>::iterator aIter=x_active.begin();
	while ( transfc.getTranslation()[0] > *xIter ) {xIter++;mIter++;rIter++;}
	x_array.insert(xIter,transfc.getTranslation()[0]);
	x_mat.insert(mIter,mdtMat);
	x_thickness.insert(tIter,2*xv);
	x_ref.insert(rIter,transfc.getTranslation()[2] );
	x_active.insert(aIter,active);
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
    // x-y plane -> y-z plane
    Trk::TrapezoidBounds* tbounds= new Trk::TrapezoidBounds(x1v,x2v,yv);
    Trk::SharedObject<const Trk::SurfaceBounds> bounds(tbounds);
    for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {      
      thickness = x_thickness[iloop];
      Trk::HomogenousLayerMaterial mdtMaterial(x_mat[iloop]);
      /*
	HepTransform3D* cTr = new HepTransform3D( (*transf) * HepRotateZ3D(-90*deg)*HepTranslateX3D(x_array[iloop]) * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      */
      HepTransform3D* cTr = new HepTransform3D( (*transf) * HepTranslateZ3D(x_array[iloop]) );
      layer = new Trk::PlaneLayer(cTr,
                                  bounds,
                                  mdtMaterial,
                                  thickness,
                                  od );
      layer->setRef(x_ref[iloop]);
      //std::cout << " reference value set for layer:"<<iloop<<","<<layer->getRef()<<std::endl;
      layers.push_back(layer);
      // make preliminary identification of active layers
      layer->setLayerType(x_active[iloop]);
    } 
    
    // create the BinnedArray
    //std::cout << "number of Mdt layers:"<<layers.size()<<std::endl;
    std::vector<Trk::SharedObject<const Trk::Layer> > layerOrder;
    std::vector<double> binSteps;
    // 
    // double minX = - volBounds->halflengthZ();
    double minX = transf->getTranslation()[0] - volBounds->halflengthZ();
    if (layers.size()) {
      //minX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
      currX = minX; 
      for (unsigned int i=0;i<layers.size();i++) { 
	const HepTransform3D ltransf = layers[i]->transform();
	layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers[i]));
	if ( i < layers.size()-1 ) { 
	  binSteps.push_back(ltransf.getTranslation()[0]+0.5*layers[i]->thickness()-currX);
	  currX = ltransf.getTranslation()[0]+0.5*layers[i]->thickness();
	}     
      }
      binSteps.push_back( transf->getTranslation()[0] + volBounds->halflengthZ() - currX );
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
							    *m_muonMaterial,
							    m_muonMagneticField,
							    mdtLayerArray,0,
							    name);         
    // std::cout << "Mdt processed with" << layers.size() << " layers" << std::endl;
    return mdt;
  }
  return 0;
}
//
const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processRpc(Trk::Volume*& vol,std::vector<const GeoVPhysVol*> gv, std::vector<HepTransform3D> transfc) const
{
  // layers correspond to DedModules and RpcModules; all substructures averaged in material properties
  std::vector<const Trk::Layer*> layers;
  for (unsigned int ic=0; ic<gv.size(); ++ic) {
    //std::cout << "processing Rpc component:"<< gv[ic]->getLogVol()->getName() <<std::endl;
    const GeoLogVol* glv = gv[ic]->getLogVol();
    const GeoShape* shape = glv->getShape();
    if (shape->type()!="Box" && shape->type()!="Trd") {
      const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (shape);
      const GeoShape* subt = 0;
      while (sub) {
        subt = sub->getOpA();
        sub = dynamic_cast<const GeoShapeSubtraction*> (subt);
      }
      shape = subt; 
    }
    if (shape && shape->type()=="Box") {
      const GeoBox* box = dynamic_cast<const GeoBox*> (shape);
      double xs = box->getXHalfLength();
      double ys = box->getYHalfLength();
      double zs = box->getZHalfLength();
      // std::cout << "dimensions:"<<box->getXHalfLength() << ","<<box->getYHalfLength() << ","<<box->getZHalfLength() << std::endl;
      // translating into layer; x dimension defines thickness
      const Trk::PlaneLayer* layer;
      double thickness=2*xs;
      Trk::OverlapDescriptor* od=0;
      Trk::RectangleBounds* rbounds = new Trk::RectangleBounds(ys,zs);
      Trk::SharedObject<const Trk::SurfaceBounds> bounds(rbounds); 
      HepTransform3D* cTr = new HepTransform3D( transfc[ic] * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      Trk::ExtendedMaterialProperties rpcMat(0.,10e8,10e8,0.,0.,0.);               // default
      if ( (glv->getName()).substr(0,3)=="Ded" ) {
        // find if material exists already
        bool found = false;
        for (unsigned int i=0;i<m_rpcDed.size();i++) {
	  if (fabs(thickness-m_rpcDed[i]->thickness())<0.001) { rpcMat = *(m_rpcDed[i]); found = true; break; }
        } 
	if (!found) {
	  double vol = 8*xs*ys*zs;
	  Trk::ExtendedMaterialProperties* rpcDed = getAveragedLayerMaterial(gv[ic],vol,2*xs);
	  m_rpcDed.push_back(rpcDed);
	  rpcMat = *rpcDed;
	}
      } else {
        //printChildren(gv[ic]);
        if (fabs(thickness-46.0)<0.001) {
          if (!m_rpc46) {
            double vol = 8*xs*ys*zs;
            m_rpc46 = getAveragedLayerMaterial(gv[ic],vol,2*xs);
          }
          rpcMat=*m_rpc46;  
        } else { *m_log << MSG::WARNING << name() << "RPC module thickness different from 46:" << thickness << endreq; }
      }
          
      Trk::HomogenousLayerMaterial rpcMaterial(rpcMat);
      layer = new Trk::PlaneLayer(cTr,
				  bounds,
				  rpcMaterial,
				  thickness,
				  od );
      layers.push_back(layer);
      // make preliminary identification of active layers
      if ((glv->getName()).substr(0,3)!="Ded" ) {
        layer->setLayerType(1);
      } else {
        layer->setLayerType(0);
      }
    } else if (shape && shape->type()=="Trd") {
      const GeoTrd* trd = dynamic_cast<const GeoTrd*> (shape);
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
        Trk::RectangleBounds* rbounds = new Trk::RectangleBounds(ys1,zs); 
        Trk::SharedObject<const Trk::SurfaceBounds> bounds(rbounds);
        HepTransform3D* cTr = new HepTransform3D( transfc[ic] * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
        Trk::ExtendedMaterialProperties rpcMat(0.,10e8,10e8,0.,0.,0.);               // default
        if ( (glv->getName()).substr(0,3)=="Ded" ) {
	  // find if material exists already
	  bool found = false;
	  for (unsigned int i=0;i<m_rpcDed.size();i++) {
	    if (fabs(thickness-m_rpcDed[i]->thickness())<0.001) { rpcMat = *(m_rpcDed[i]); found = true; break; }
	  } 
	  if (!found) {
	    double vol = 8*xs1*ys1*zs;
	    Trk::ExtendedMaterialProperties* rpcDed = getAveragedLayerMaterial(gv[ic],vol,2*xs1);
	    m_rpcDed.push_back(rpcDed);
	    rpcMat = *rpcDed;
	  }
          // create Ded layer
	  Trk::HomogenousLayerMaterial rpcMaterial(rpcMat);
	  layer = new Trk::PlaneLayer(cTr, bounds, rpcMaterial, thickness, od );
          layers.push_back(layer);
	  layer->setLayerType(0);
        } else {
          // RPC layer; step one level below to resolve strip planes
          //printChildren(gv[ic]);
          unsigned int ngc = gv[ic]->getNChildVols();
          for (unsigned int igc=0; igc<ngc; igc++) {
	    HepTransform3D trgc;
             if (transfc[ic].getTranslation()[0]>vol->center()[0]) {
               trgc = HepRotateZ3D(180*deg)*(gv[ic]->getXToChildVol(igc)); 
             } else {
               trgc = gv[ic]->getXToChildVol(igc);
             }
             const GeoVPhysVol* gcv = &(*(gv[ic]->getChildVol(igc)));
             const GeoLogVol* gclv = gcv->getLogVol();
             const GeoShape* lshape = gclv->getShape();
             while (lshape->type()=="Subtraction") {
	       const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (lshape); 
	       lshape = sub->getOpA(); 
	     } 
	     const GeoTrd* gtrd = dynamic_cast<const GeoTrd*> (lshape);
	     double gx = gtrd->getXHalfLength1();
	     double gy = gtrd->getYHalfLength1();
	     double gz = gtrd->getZHalfLength();
            
             if ( (gclv->getName()).substr(0,6)=="RPC_AL" ) {
	       if (fabs(gx-5.0) < 0.001) {
		 if (!m_rpcExtPanel) {
		   double vol = 8*gx*gy*gz;
                   m_rpcExtPanel = getAveragedLayerMaterial(gcv,vol,2*gx);
                 }
                 rpcMat=*m_rpcExtPanel;
	       } else if (fabs(gx - 4.3) < 0.001) {
		 if (!m_rpcMidPanel) {
		   double vol = 8*gx*gy*gz;
                   m_rpcMidPanel = getAveragedLayerMaterial(gcv,vol,2*gx);
                 }
                 rpcMat=*m_rpcMidPanel;
	       } else {
		 *m_log << MSG::WARNING << name() << "unknown RPC panel:" << gx << endreq;               
               }
	       // create Rpc panel layers 
               thickness = 2*gx;
	       Trk::HomogenousLayerMaterial rpcMaterial(rpcMat);
	       layer = new Trk::PlaneLayer(new HepTransform3D(HepTranslate3D(trgc.getTranslation())*(*cTr)),
					   bounds, rpcMaterial, thickness, od );
	       layers.push_back(layer);
	       layer->setLayerType(0);
	     } else if  ( (gclv->getName())=="Rpclayer" ) {
               if ( fabs(gx-6.85)>0.001 )  *m_log << MSG::WARNING << name() << " unusual thickness of RPC layer:" << 2*gx << endreq;
	       if (!m_rpcLayer) {                   
		 double vol = 8*gx*gy*gz;
                 // material allocated to two strip planes ( gas volume suppressed )
		 m_rpcLayer = getAveragedLayerMaterial(gcv,vol,2*gx);
	       }
	       rpcMat=*m_rpcLayer;
               // define 1 layer for 2 strip planes
               thickness = 2*gx;
	       Trk::HomogenousLayerMaterial rpcMaterial(rpcMat);
	       layer = new Trk::PlaneLayer(new HepTransform3D(HepTranslate3D(trgc.getTranslation())*(*cTr)),
					   bounds, rpcMaterial, thickness, od );
	       layers.push_back(layer);
	       layer->setLayerType(1);               
             } else {
	       *m_log << MSG::WARNING << name()  << "unknown RPC component? " << gclv->getName() << endreq;
             }
	  }
          delete cTr;
	}
      } else {
         *m_log << MSG::WARNING << name() << "RPC true trapezoid layer, not coded yet" <<endreq;
      }
    } else {
      *m_log << MSG::WARNING << name() << "RPC layer shape not recognized" <<endreq;
    }
  } // end loop over Modules

  std::vector<const Trk::Layer*>* rpcLayers = new std::vector<const Trk::Layer*>(layers); 
  std::string name="RPC";
  const Trk::TrackingVolume* rpc= new Trk::TrackingVolume(*vol,
                                                          *m_muonMaterial,
                                                          m_muonMagneticField,
                                                          rpcLayers,
                                                          name);         
  *m_log << MSG::DEBUG << " Rpc component volume processed with" << layers.size() << " layers"  << endreq;
  return rpc;
}
//

const Trk::TrackingVolume* Muon::MuonStationTypeBuilder::processSpacer(Trk::Volume*& vol,std::vector<const GeoVPhysVol*> gv, std::vector<HepTransform3D> transf) const
{
  // spacers: one level below, assumed boxes
  std::vector<const Trk::Layer*> layers;
  // resolve child volumes
  std::vector<const GeoVPhysVol*>::iterator vIter = gv.begin();
  std::vector<HepTransform3D>::iterator tIter = transf.begin();
  while (vIter != gv.end() ) {
    if ((*vIter)->getNChildVols()) {
      for (unsigned int ich=0; ich<(*vIter)->getNChildVols();++ich) {
	gv.push_back(&(*((*vIter)->getChildVol(ich))));
        transf.push_back(HepTransform3D( (*tIter)*(*vIter)->getXToChildVol(ich)));  
        gv.erase(vIter); transf.erase(tIter);
      }
    } else { vIter++; tIter++; }
  } 
  // translate into layers
  for (unsigned int ic=0; ic<gv.size(); ++ic) {
    const GeoLogVol* clv = gv[ic]->getLogVol();        
    // std::cout << "Spacer component:"<<ic<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transf[ic].getTranslation()<<","<<transform.getTranslation()<<std::endl; 
    Trk::ExtendedMaterialProperties cmat = m_materialConverter->convertExtended(clv->getMaterial());
    Trk::OverlapDescriptor* od=0;
    if (clv->getShape()->type()=="Box") {
      const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
      double xs = box->getXHalfLength();
      double ys = box->getYHalfLength();
      double zs = box->getZHalfLength();
      // translating into layer; find minimal size
      const Trk::PlaneLayer* layer;
      Trk::RectangleBounds* rbounds=0;
      double thickness=0.;
      HepTransform3D* cTr = 0;
      if (zs <= xs && zs <= ys ) { // x-y plane
	rbounds = new Trk::RectangleBounds(xs,ys); 
	thickness = zs;
        cTr = new HepTransform3D(transf[ic]);
      } else if (xs <= ys && xs <= zs ) { // x-y plane -> y-z plane
	rbounds = new Trk::RectangleBounds(ys,zs); 
	thickness = xs; 
	cTr = new HepTransform3D( transf[ic] * HepRotateY3D(90*deg) * HepRotateZ3D(90*deg));
      } else {  // x-y plane -> x-z plane
	rbounds = new Trk::RectangleBounds(xs,zs); 
	thickness = ys; 
	cTr = new HepTransform3D( transf[ic] * HepRotateX3D(90*deg));
      } 
      Trk::ExtendedMaterialProperties material(thickness,cmat.x0(),cmat.l0(),cmat.averageA(),
					       cmat.averageZ(),cmat.averageRho());  
      Trk::HomogenousLayerMaterial spacerMaterial(material);
      Trk::SharedObject<const Trk::SurfaceBounds> bounds(rbounds);
      
      layer = new Trk::PlaneLayer(cTr,
				  bounds,
				  spacerMaterial,
				  thickness,
				  od,0 );
      layers.push_back(layer);
    } else if ( clv->getShape()->type()=="Subtraction") {
      const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (clv->getShape());
      if (sub && sub->getOpA()->type()=="Box" && sub->getOpB()->type()=="Box") {
        // LB
        const GeoBox* boxA = dynamic_cast<const GeoBox*> (sub->getOpA());
        const GeoBox* boxB = dynamic_cast<const GeoBox*> (sub->getOpB());
	Trk::RectangleBounds* rbounds = new Trk::RectangleBounds(boxA->getYHalfLength(),boxA->getZHalfLength());
        double thickness = 0.5*(boxA->getXHalfLength()-boxB->getXHalfLength());
        double shift     = 0.5*(boxA->getXHalfLength()+boxB->getXHalfLength());
	Trk::ExtendedMaterialProperties material;
	Trk::HomogenousLayerMaterial spacerMaterial;
	Trk::SharedObject<const Trk::SurfaceBounds> bounds;
        if (thickness>0.) {
	  material = Trk::ExtendedMaterialProperties(thickness,cmat.x0(),cmat.l0(),cmat.averageA(),
						   cmat.averageZ(),cmat.averageRho());  
	  spacerMaterial = Trk::HomogenousLayerMaterial(material);
	  bounds = Trk::SharedObject<const Trk::SurfaceBounds>(rbounds);
	  Trk::PlaneLayer* layx = new Trk::PlaneLayer(new HepTransform3D( transf[ic]*HepTranslateX3D(shift)*HepRotateY3D(90*deg)*HepRotateZ3D(90*deg) ), bounds, spacerMaterial, thickness, od, 0 );
	  Trk::PlaneLayer* layxx = new Trk::PlaneLayer(*layx,HepTranslateX3D(-2*shift));
	  layers.push_back(layx)  ; 
	  layers.push_back(layxx) ; layxx->setLayerType(0);
	}
	rbounds = new Trk::RectangleBounds(boxB->getXHalfLength(),boxA->getZHalfLength());
        thickness = 0.5*(boxA->getYHalfLength()-boxB->getYHalfLength());
        if (thickness>0.) {
	  material = Trk::ExtendedMaterialProperties(thickness,cmat.x0(),cmat.l0(),cmat.averageA(),
						     cmat.averageZ(),cmat.averageRho());  
	  spacerMaterial = Trk::HomogenousLayerMaterial(material);
	  shift     = 0.5*(boxA->getYHalfLength()+boxB->getYHalfLength());
	  bounds = Trk::SharedObject<const Trk::SurfaceBounds>(rbounds);
	  Trk::PlaneLayer* lay = new Trk::PlaneLayer(new HepTransform3D( transf[ic]*HepTranslateY3D(shift)*HepRotateX3D(90*deg)), bounds, spacerMaterial, thickness, od, 0 );
	  Trk::PlaneLayer* layy = new Trk::PlaneLayer(*lay,HepTranslateY3D(-2*shift));
	  layers.push_back(lay)  ; 
	  layers.push_back(layy) ; layy->setLayerType(0);
	}
	rbounds = new Trk::RectangleBounds(boxB->getXHalfLength(),boxB->getYHalfLength());
        thickness = 0.5*(boxA->getZHalfLength()-boxB->getZHalfLength());
        if (thickness>0.) {
	  material = Trk::ExtendedMaterialProperties(thickness,cmat.x0(),cmat.l0(),cmat.averageA(),cmat.averageZ(),cmat.averageRho());  
	  spacerMaterial = Trk::HomogenousLayerMaterial(material);
	  shift     = 0.5*(boxA->getZHalfLength()+boxB->getZHalfLength());
	  bounds = Trk::SharedObject<const Trk::SurfaceBounds>(rbounds);
	  Trk::PlaneLayer* layz = new Trk::PlaneLayer(new HepTransform3D( transf[ic]*HepTranslateZ3D(shift)), bounds, spacerMaterial, thickness, od, 0 );
	  Trk::PlaneLayer* layzz = new Trk::PlaneLayer(*layz,HepTranslateZ3D(-2*shift));
	  layers.push_back(layz)  ;
	  layers.push_back(layzz) ; layzz->setLayerType(0);
	}
      } else if (sub) {
	std::vector<std::pair<const GeoShape*,HepTransform3D> > subVs;
	const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (sub->getOpB());
	if (shift) subVs.push_back(std::pair<const GeoShape*,HepTransform3D>(shift->getOp(),shift->getX()));
        const GeoShape* shape = sub->getOpA();
        while (shape->type()=="Subtraction") {
          const GeoShapeSubtraction* subtr = dynamic_cast<const GeoShapeSubtraction*> (shape);
          const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (subtr->getOpB());
          if (shift) subVs.push_back(std::pair<const GeoShape*,HepTransform3D>(shift->getOp(),shift->getX()));
          shape = subtr->getOpA();  
	}
        const GeoBox* box = dynamic_cast<const GeoBox*> (shape);
        if (box && subVs.size()==4) {
	  Trk::Volume* v1 = 0; Trk::Volume* v2 = 0; Trk::VolumeExcluder* volExcl=0; 
          const GeoBox* sb1 = dynamic_cast<const GeoBox*> (subVs[0].first);
          if (sb1) v1 = new Trk::Volume(new HepTransform3D(subVs[0].second),
					new Trk::CuboidVolumeBounds(sb1->getXHalfLength(),
								    sb1->getYHalfLength(),
								    sb1->getZHalfLength()));
          const GeoBox* sb2 = dynamic_cast<const GeoBox*> (subVs[1].first);
          if (sb2) v2 = new Trk::Volume(new HepTransform3D(subVs[1].second),
					new Trk::CuboidVolumeBounds(sb2->getXHalfLength(),
								    sb2->getYHalfLength(),
								    sb2->getZHalfLength()));
           
          const GeoBox* boxB = dynamic_cast<const GeoBox*> (subVs[2].first);
          if (boxB && v1 && v2) {
	    Trk::RectangleBounds* rbounds = new Trk::RectangleBounds(box->getYHalfLength(),box->getZHalfLength());
	    double thickness = 0.5*(box->getXHalfLength()-boxB->getXHalfLength());
	    double shift     = 0.5*(box->getXHalfLength()+boxB->getXHalfLength());
	    Trk::SharedObject<const Trk::SurfaceBounds> bounds(rbounds);
	    Trk::Volume* cVol = new Trk::Volume(new HepTranslateX3D(-shift),
						new Trk::CombinedVolumeBounds(v1,v2,false));
            volExcl=new Trk::VolumeExcluder(cVol);
	    Trk::SubtractedPlaneSurface* subPlane=new Trk::SubtractedPlaneSurface(Trk::PlaneSurface(new HepTransform3D(transf[ic]*HepTranslateX3D(shift)*HepRotateY3D(90*deg)*HepRotateZ3D(90*deg)), bounds), volExcl, false);
	    Trk::ExtendedMaterialProperties material(thickness,cmat.x0(),cmat.l0(),cmat.averageA(),
						     cmat.averageZ(),cmat.averageRho());  
	    Trk::HomogenousLayerMaterial spacerMaterial(material);
	    Trk::SubtractedPlaneLayer* layx = new Trk::SubtractedPlaneLayer(subPlane,spacerMaterial, thickness, od, 0 );
            delete subPlane;
	    Trk::SubtractedPlaneLayer* layxx = new Trk::SubtractedPlaneLayer(*layx,HepTranslateX3D(-2*shift));
	    layers.push_back(layx) ; 
	    layers.push_back(layxx) ; layxx->setLayerType(0);

	    rbounds = new Trk::RectangleBounds(boxB->getXHalfLength(),box->getZHalfLength());
	    thickness = subVs[2].second.getTranslation().mag();
	    Trk::VolumeExcluder* volEx=new Trk::VolumeExcluder(new Trk::Volume(*cVol, 
									       HepTranslateX3D(2*shift)));
	    bounds = Trk::SharedObject<const Trk::SurfaceBounds>(rbounds);
	    subPlane = new Trk::SubtractedPlaneSurface(Trk::PlaneSurface(new HepTransform3D(transf[ic]*HepRotateX3D(90*deg)), bounds), volEx, false);
	    material = Trk::ExtendedMaterialProperties(thickness,cmat.x0(),cmat.l0(),cmat.averageA(),
							   cmat.averageZ(),cmat.averageRho());  
	    spacerMaterial = Trk::HomogenousLayerMaterial(material);
	    Trk::SubtractedPlaneLayer* lay = new Trk::SubtractedPlaneLayer(subPlane,spacerMaterial, thickness, od, 0 );
            delete subPlane;
	    layers.push_back(lay) ; 
	  }          
	}
      } else {
	std::cout << "unresolved spacer component "<< clv->getName() << std::endl; 
      }    
    } else {
      std::cout << "unresolved spacer component "<< clv->getName() << std::endl; 
    }
  }

  std::vector<const Trk::Layer*>::iterator lIt = layers.begin();
  for ( ; lIt!=layers.end();lIt++) if ((*lIt)->thickness()<0.) layers.erase(lIt); 

  std::vector<const Trk::Layer*>* spacerLayers = new std::vector<const Trk::Layer*>(layers); 

  std::string name="Spacer";
  const Trk::TrackingVolume* spacer= new Trk::TrackingVolume(*vol,
							     *m_muonMaterial,
							     m_muonMagneticField,
							     spacerLayers,
							     name);         
  //std::cout << "spacer processed with" << layers.size() << " layers" << std::endl;
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
    //std::cout << "First CSC component:"<< clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transform.getTranslation()<<std::endl;
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
         xMin =  trdA->getYHalfLength1();
         xMed =  trdA->getYHalfLength2();
         y1   =  trdA->getZHalfLength();
         z    =  trdA->getXHalfLength1();
        }
        if (uni->getOpB()->type()=="Shift") {
          const GeoShapeShift* sh = dynamic_cast<const GeoShapeShift*> (uni->getOpB());
          //std::cout << sh->getOp()->type()<<","<<(sh->getX()).getTranslation()<<std::endl;
          const GeoTrd* trdB = dynamic_cast<const GeoTrd*> (sh->getOp());          
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
                                                                    *m_muonMaterial,
                                                                    m_muonMagneticField,
                                                                    cscLayerArray,0,
                                                                    compName[ic]);                    
       delete compVol;
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
     HepTransform3D* cTr = new HepTransform3D(HepRotateY3D(90*deg)*HepRotateZ3D(90*deg)*HepTranslateY3D(envY1-envY2));
     envelope = new Trk::Volume(cTr,cscBounds);
     //std::cout << "envelope created"<<std::endl;
     // components
     double xCurr = -xTotal;
     for (unsigned int ic = 0; ic< xSizes.size(); ic++) {
       // component volumes follow the envelope dimension
       xCurr += xSizes[ic];
       //HepTransform3D* compTr = new HepTransform3D(HepRotateZ3D(90*deg)*HepTranslateZ3D(envY1-envY2)*HepTranslateY3D(xCurr));
       HepTransform3D* compTr = new HepTransform3D(HepRotateY3D(90*deg)*HepRotateZ3D(90*deg)*HepTranslateY3D(envY1-envY2)*HepTranslateZ3D(xCurr));
       compBounds = new Trk::DoubleTrapezoidVolumeBounds(xMin,envXMed,xMax,envY1,envY2,xSizes[ic]); 
       const Trk::LayerArray* cscLayerArray = processCSCDiamondComponent(compGeoVol[ic],compBounds,compTr); 
       Trk::Volume* compVol = new Trk::Volume(compTr,compBounds);
       const Trk::TrackingVolume* compTV = new Trk::TrackingVolume( *compVol,
                                                                    *m_muonMaterial,
                                                                    m_muonMagneticField,
                                                                    cscLayerArray,0,
                                                                    compName[ic]);                    
       delete compVol;
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
								   *m_muonMaterial,
								   m_muonMagneticField,
								   0,compArray,
								   name);                    
 
 delete envelope;
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
     //std::cout << "tgc name:" << tgc_name << std::endl; 
     HepTransform3D transform = mv->getXToChildVol(ich);        
     //std::cout << "TGC component:"<<ich<<":" << clv->getName() <<", made of "<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<","<<transform.getTranslation()<<std::endl;
     const GeoShape* baseShape = clv->getShape();
     if (baseShape->type()=="Subtraction") {
       const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (baseShape);
       if (sub) baseShape = sub->getOpA();
     }
     
     if (baseShape->type()=="Trd") {
       const GeoTrd* trd = dynamic_cast<const GeoTrd*> (baseShape);
       double x1 = trd->getXHalfLength1();
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
                                                                        *m_muonMaterial,
                                                                        m_muonMagneticField,
									tgcLayerArray,0,
                                                                        tgc_name);                    

       delete envelope;
       if (tgc_station) tgc_stations.push_back(tgc_station);       
      
     } else {
       std::cout << "TGC component not trapezoid ?" << std::endl;
     }
  }
 return tgc_stations;    

}

void Muon::MuonStationTypeBuilder::printChildren(const GeoVPhysVol* pv) const
{
  // subcomponents
  unsigned int nc = pv->getNChildVols();
  for (unsigned int ic=0; ic<nc; ic++) {
    HepTransform3D transf = pv->getXToChildVol(ic);
 
    //
    /*
    std::cout << " dumping transform to subcomponent" << std::endl;
    std::cout << transf[0][0]<<"," <<transf[0][1]<<"," <<transf[0][2]<<","<<transf[0][3] << std::endl;
    std::cout << transf[1][0]<<"," <<transf[1][1]<<"," <<transf[1][2]<<","<<transf[1][3] << std::endl;
    std::cout << transf[2][0]<<"," <<transf[2][1]<<"," <<transf[2][2]<<","<<transf[2][3] << std::endl;
    */
    //
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    std::cout << "  ";
    std::cout << "subcomponent:"<<ic<<":"<<clv->getName()<<", made of"<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<< ","<< transf.getTranslation()<<std::endl;
	 
          if ( clv->getShape()->type()=="Trd") {
	      const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
	      //
	      std::cout << "dimensions:"<< trd->getXHalfLength1() <<","
                                        << trd->getXHalfLength2() <<","  
                                        << trd->getYHalfLength1() <<","  
                                        << trd->getYHalfLength2() <<","  
			                << trd->getZHalfLength() <<std::endl; 
	      //
          } 
          if ( clv->getShape()->type()=="Box") {
	      const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
	      //
	        std::cout << "dimensions:"<< box->getXHalfLength() <<","
                                        << box->getYHalfLength() <<","  
			                << box->getZHalfLength() <<std::endl; 
	      //
          } 

	  printChildren(cv);
    }  
   

}

double Muon::MuonStationTypeBuilder::get_x_size(const GeoVPhysVol* pv) const
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
    if (type=="Subtraction") {
       xh = decodeX(clv->getShape());
    }
 
    xlow = fmin(xlow,(transf.getTranslation())[0]-xh);
    xup  = fmax(xup ,(transf.getTranslation())[0]+xh);
  }  
  // std::cout << "x size:" << xlow <<","<< xup << std::endl;
   
  return fmax( -xlow, xup ); 
}

Trk::ExtendedMaterialProperties* Muon::MuonStationTypeBuilder::getAveragedLayerMaterial( const GeoVPhysVol* pv, double volume, double thickness) const
{  

  *m_log << MSG::DEBUG << name() << "::getAveragedLayerMaterial:processing "<< std::endl; 
  // loop through the whole hierarchy; collect material
  Trk::ExtendedMaterialProperties* total=new Trk::ExtendedMaterialProperties(0.,10.e8,10.e8,0.,0.,0.);
  collectMaterial( pv, total, volume/thickness);
  *m_log << MSG::DEBUG << name() << " combined material thickness: "<< total->thickness() << std::endl; 
  *m_log << MSG::DEBUG << name() << " actual layer thickness: "<< thickness << std::endl; 
  // scaled material properties to the actual layer thickness
  if (total->thickness() > 0 ) { 
    double scale = thickness/total->thickness();
    //std::cout << "scaling factor:"<<scale<<std::endl;
    Trk::ExtendedMaterialProperties* scaled =
      new Trk::ExtendedMaterialProperties(scale*total->thickness(),
					  scale*total->x0(),
					  scale*total->l0(),
					  total->averageA(),
					  total->averageZ(),
					  total->averageRho()/scale);
    delete total;    
    *m_log<< MSG::DEBUG << "averaged material:d,x0,dInX0:"<<scaled->thickness()<<","<<scaled->x0()<<","<<scaled->thicknessInX0()<<endreq;

    return scaled;
  }  
  *m_log << MSG::DEBUG << name() << "::getAverageLayerMaterial returns 0 " << std::endl; 
  return 0;
}

 void Muon::MuonStationTypeBuilder::collectMaterial(const GeoVPhysVol* pv, Trk::ExtendedMaterialProperties*& matProp,double sf) const
{
  // sf is surface of the new layer used to calculate the average 'thickness' of components
  // number of child volumes
  unsigned int nc = pv->getNChildVols();
  // add current volume 
  const GeoLogVol* lv = pv->getLogVol(); 
  //std::cout << "component:"<<lv->getName()<<", made of"<<lv->getMaterial()->getName()<<","<<lv->getShape()->type()<<std::endl;
     
  if ( lv->getMaterial()->getName() != "Air" && (lv->getName()).substr(0,1)!="T") {
    // get material properties from GeoModel
    Trk::ExtendedMaterialProperties newMP= m_materialConverter->convertExtended( lv->getMaterial() ); 
    // current volume
    double vol = getVolume(lv->getShape());
    // subtract children volumes
    for (unsigned int ic=0; ic<nc; ic++) {
      const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
      if ( getVolume(cv->getLogVol()->getShape()) > vol ) {
	//std::cout << "WARNING:collect material : child volume bigger than mother volume" << std::endl; 
      } else {
	vol = vol - getVolume(cv->getLogVol()->getShape()); 
      }
    }
    double d = vol / sf;
    double totalThickness = d + matProp->thickness();
    double fold = matProp->thickness()/totalThickness;
    double fnew = d/totalThickness; 
  
    //std::cout << "old material:f,d,x0:"<< fold<<","<<matProp->thickness()<<","<<matProp->x0()<<std::endl;
    //std::cout << "new material:f,d,x0:"<< fnew<<","<<newMP.thickness()<<","<<newMP.x0()<<std::endl;
 
    Trk::ExtendedMaterialProperties* newUpdate=new Trk::ExtendedMaterialProperties(totalThickness,
                                              1./(fnew/newMP.x0()+fold/matProp->x0()),
                                              1./(fnew/newMP.l0()+fold/matProp->l0()),
                                              fnew*newMP.averageA()+fold*matProp->averageA(),
                                              fnew*newMP.averageZ()+fold*matProp->averageZ(),
                                              fnew*newMP.averageRho()+fold*matProp->averageRho());  

    delete matProp;
    matProp = newUpdate;
    //std::cout << "sum material:d,x0:"<<matProp->thickness()<<","<<matProp->x0()<<std::endl;
  } 
  // subcomponents
  // skip children volume if we deal with G10 ( not correctly described )
  //if ( lv->getName() != "G10" ) { 
  for (unsigned int ic=0; ic<nc; ic++) {
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    collectMaterial( cv, matProp, sf);
  }
}

double Muon::MuonStationTypeBuilder::getVolume( const GeoShape* shape) const {
  //
  double volume = 0.;
  
  if (shape->type()=="Shift" ) {
    const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (shape);
    volume = getVolume(shift->getOp());
  } else if (shape->type()=="Subtraction" ) {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction* > (shape);
    double volA = getVolume(sub->getOpA());
    double volB = getVolume(sub->getOpB());
    // protection against subtraction of large volumes 
    if (volA > volB){
      volume = volA - volB;
    } else {
      volume = volA;
    }
  } else if (shape->type()=="Union") {  
    const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion* > (shape);
    double volA = getVolume(uni->getOpA());
    double volB = getVolume(uni->getOpB());
    volume = volA+volB; 
  } else {
    volume = shape->volume();
  }
  return volume;
} 

const Trk::LayerArray* Muon::MuonStationTypeBuilder::processCSCTrdComponent(const GeoVPhysVol*& pv, Trk::TrapezoidVolumeBounds*& compBounds, HepTransform3D*& transf) const {

  // tolerance
  // double tol = 0.001;
  std::string name = pv->getLogVol()->getName();
  //std::cout << "processing CSC component, number of children volumes:"<< name << "," << pv->getNChildVols() <<std::endl; 
  //printChildren(pv);
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::ExtendedMaterialProperties> x_mat;
  std::vector<double> x_thickness;
  std::vector<int> x_active;
  double currX = -100000;
  // while waiting for better suggestion, define a single material layer
  Trk::ExtendedMaterialProperties matCSC(0.,10e8,10e8,0.,0.,0.);
  double thickness =2* compBounds->halflengthZ();
  double minX = compBounds->minHalflengthX();
  double maxX = compBounds->maxHalflengthX();
  double halfY = compBounds->halflengthY();
  double halfZ = compBounds->halflengthZ();
  if (name.substr( name.size()-5,5 ) == "CSC01" ) {
    if (!m_matCSC01 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matCSC01 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matCSC = *m_matCSC01; 
    // retrieve number of gas gaps and their position -> turn them into active layers
    // step 1 level below
    const GeoVPhysVol* cv1 = &(*(pv->getChildVol(0)));
    for (unsigned int ic=0; ic < cv1->getNChildVols(); ic++) {
      HepTransform3D transfc = cv1->getXToChildVol(ic);
      const GeoVPhysVol* cv = &(*(cv1->getChildVol(ic)));
      const GeoLogVol* clv = cv->getLogVol();
      if ( clv->getName() == "CscArCO2" ) {
        double xl = transfc.getTranslation()[0];
        if (x_array.size() == 0 || xl >= x_array.back()) { 
          x_array.push_back(xl);
        } else {
          unsigned int ix = 0; 
          while ( ix < x_array.size() && x_array[ix] < xl  ) {ix++;}
          x_array.insert(x_array.begin()+ix,xl);
        } 
      }
    } 
    //
    //std::cout << " Csc multilayer has " << x_array.size() << " active layers" << std::endl; 
    // TO DO: in the following, the material should be scaled to the actual layer thickness 
    if (x_array.size()==0) {
      x_array.push_back(0.);
      x_mat.push_back(matCSC);
      x_thickness.push_back(thickness );
      x_active.push_back(1);
    } else if (x_array.size()==1) {
      double xthick = 2*fmin(x_array[0]+halfZ,halfZ-x_array[0]);
      double scale = xthick/thickness;
      Trk::ExtendedMaterialProperties xmatCSC(scale*matCSC.thickness(),scale*matCSC.x0(),scale*matCSC.l0(),
					      matCSC.averageA(),matCSC.averageZ(),matCSC.averageRho()/scale); 
      x_mat.push_back(xmatCSC);
      x_thickness.push_back(xthick);
      x_active.push_back(1); 
    } else {
      double currX = -halfZ;
      for (unsigned int il=0; il < x_array.size(); il++) {
        double xthick;
	if (il<x_array.size()-1) {
          xthick = 2*fmin(x_array[il]-currX,0.5*(x_array[il+1]-x_array[il])) ;
	} else {
          xthick = 2*fmin(x_array[il]-currX,halfZ-x_array[il]);
        }
	x_thickness.push_back(xthick);
	double scale = xthick/thickness;
	Trk::ExtendedMaterialProperties xmatCSC(scale*matCSC.thickness(),scale*matCSC.x0(),scale*matCSC.l0(),
					      matCSC.averageA(),matCSC.averageZ(),matCSC.averageRho()/scale); 
        x_mat.push_back(xmatCSC);
        currX = x_array[il]+0.5*x_thickness.back();
        x_active.push_back(1);
      }
    }
    //for (unsigned int il=0; il < x_array.size(); il++) {
    //  std::cout << "csc active layers:" << il << "," << x_array[il]<< ","<<x_thickness[il] << std::endl;
    //}    
  } 
  if (name == "CSCspacer" ) {
    if (!m_matCSCspacer1 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matCSCspacer1 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    x_array.push_back(0.);
    x_mat.push_back(*m_matCSCspacer1);
    x_thickness.push_back(thickness );
    x_active.push_back(0);
  } 
  // create layers
  const Trk::PlaneLayer* layer;
  Trk::OverlapDescriptor* od=0;
  Trk::TrapezoidBounds* tbounds= new Trk::TrapezoidBounds(minX,maxX,halfY); 
  Trk::SharedObject<const Trk::SurfaceBounds> bounds(tbounds);
  for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
    HepTransform3D* cTr = new HepTransform3D(  HepTranslateX3D(x_array[iloop]) * (*transf) ); // this won't work for multiple layers !!! //
    Trk::HomogenousLayerMaterial cscMaterial(x_mat[iloop]);  
    layer = new Trk::PlaneLayer(cTr,
                                bounds,
                                cscMaterial,
                                x_thickness[iloop],
                                od );
    layers.push_back(layer);
    // make preliminary identification of active layers
    layer->setLayerType(x_active[iloop]);

    //std::cout << "CSC layer built ok"<<std::endl;

  }

  // create the BinnedArray
  //std::cout << "number of Csc layers:"<<layers.size()<<std::endl;
  std::vector<Trk::SharedObject<const Trk::Layer> > layerOrder;
  std::vector<double> binSteps;
  double xShift = transf->getTranslation()[0];
  double lowX = - compBounds->halflengthZ()+xShift ;

  if (layers.size()) {
    // lowX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = lowX; 
     for (unsigned int i=0;i<layers.size()-1;i++) { 
       const HepTransform3D ltransf = HepTranslateX3D(x_array[i]);
       layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers[i]));
       binSteps.push_back(ltransf.getTranslation()[0]+0.5*layers[i]->thickness()-currX+xShift);
       currX = ltransf.getTranslation()[0]+0.5*layers[i]->thickness()+xShift;     
       //std::cout << "binstep:"<<i <<","<<binSteps[i]<<std::endl;
     }
     const HepTransform3D ltransf = HepTranslateX3D(x_array.back());
     layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers.back()));
     binSteps.push_back(compBounds->halflengthZ()-currX+xShift);
     //std::cout << "binstep:"<<layers.size() <<","<<binSteps.back()<<std::endl;
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( lowX, new std::vector<double>(binSteps));
  Trk::LayerArray* cscLayerArray = 0;
  cscLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     

  return cscLayerArray;

} 

const Trk::LayerArray* Muon::MuonStationTypeBuilder::processCSCDiamondComponent(const GeoVPhysVol*& pv, Trk::DoubleTrapezoidVolumeBounds*& compBounds, HepTransform3D*& transf) const {

  // tolerance
  // double tol = 0.001;
  std::string name = pv->getLogVol()->getName();
  //std::cout << "processing CSC component, number of children volumes:"<< name << "," << pv->getNChildVols() <<std::endl; 
  //printChildren(pv);
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::ExtendedMaterialProperties> x_mat;
  std::vector<double> x_thickness;
  std::vector<int> x_active;
  double currX = -100000;
  // while waiting for better suggestion, define a single material layer
  Trk::ExtendedMaterialProperties matCSC(0.,10e8,10e8,0.,0.,0.);
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
    // retrieve number of gas gaps and their position -> turn them into active layers
    // step 1 level below
    const GeoVPhysVol* cv1 = &(*(pv->getChildVol(0)));
    for (unsigned int ic=0; ic < cv1->getNChildVols(); ic++) {
      HepTransform3D transfc = cv1->getXToChildVol(ic);
      const GeoVPhysVol* cv = &(*(cv1->getChildVol(ic)));
      const GeoLogVol* clv = cv->getLogVol();
      if ( clv->getName() == "CscArCO2" ) {
        double xl = transfc.getTranslation()[0];
        if (x_array.size() == 0 || xl >= x_array.back()) { 
          x_array.push_back(xl);
        } else {
          unsigned int ix = 0; 
          while ( ix < x_array.size() && x_array[ix] < xl  ) {ix++;}
          x_array.insert(x_array.begin()+ix,xl);
        } 
      }
    } 
    //
    //std::cout << " Csc multilayer has " << x_array.size() << " active layers" << std::endl; 
    // TO DO: in the following, the material should be scaled to the actual layer thickness 
    if (x_array.size()==0) {
      x_array.push_back(0.);
      x_mat.push_back(matCSC);
      x_thickness.push_back(thickness );
      x_active.push_back(1);
    } else if (x_array.size()==1) {
      x_mat.push_back(matCSC);
      x_thickness.push_back(2*fmin(x_array[0]+halfZ,halfZ-x_array[0]));
      x_active.push_back(1);
    } else {
      double currX = -halfZ;
      for (unsigned int il=0; il < x_array.size(); il++) {
        double xthick = 0.;        
	if (il<x_array.size()-1) {
          xthick = 2*fmin(x_array[il]-currX,0.5*(x_array[il+1]-x_array[il])); 
	  x_thickness.push_back(xthick);
	} else {
          xthick = 2*fmin(x_array[il]-currX,halfZ-x_array[il]);
	  x_thickness.push_back(xthick);
        }
        double scale = xthick/thickness;
        Trk::ExtendedMaterialProperties xmatCSC(scale*matCSC.thickness(),scale*matCSC.x0(),
						scale*matCSC.l0(),matCSC.averageA(),
						matCSC.averageZ(),matCSC.averageRho()/scale); 
        x_mat.push_back(xmatCSC);
        currX = x_array[il]+0.5*x_thickness.back();
        x_active.push_back(1); 
      }
    }
    //for (unsigned int il=0; il < x_array.size(); il++) {
    //  std::cout << "csc active layers:" << il << "," << x_array[il]<< ","<<x_thickness[il] << std::endl;
    //}    
  } 
  if (name == "CSCspacer" ) {
    if (!m_matCSCspacer2 ) { 
      double vol = ( (minX + medX)*2*halfY1+(medX+maxX)*2*halfY2 ) * thickness;
      m_matCSCspacer2 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matCSC = *m_matCSCspacer2; 
    x_array.push_back(0.);
    x_mat.push_back(matCSC);
    x_thickness.push_back(thickness );
    x_active.push_back(0);
  } 
  // create layers
  const Trk::PlaneLayer* layer;
  Trk::OverlapDescriptor* od=0;
  Trk::DiamondBounds* dbounds= new Trk::DiamondBounds(minX,medX,maxX,halfY1,halfY2); ;
  Trk::SharedObject<const Trk::SurfaceBounds> bounds(dbounds);
  for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
    HepTransform3D* cTr = new HepTransform3D( HepTranslateX3D(x_array[iloop]) * (*transf) ); // this won't work for multiple layers !!! //
    Trk::HomogenousLayerMaterial cscMaterial(x_mat[iloop]);  
    layer = new Trk::PlaneLayer(cTr,
                                bounds,
                                cscMaterial,
                                x_thickness[iloop],
                                od );
    layers.push_back(layer);
    // make preliminary identification of active layers
    layer->setLayerType(x_active[iloop]);
    //std::cout << "CSC layer built ok"<<std::endl;
  }
   
  // create the BinnedArray
  //std::cout << "number of Csc layers:"<<layers.size()<<std::endl;
  std::vector<Trk::SharedObject<const Trk::Layer> > layerOrder;
  std::vector<double> binSteps;
  double xShift = transf->getTranslation()[0];
  double lowX = - compBounds->halflengthZ()+xShift ;

  if (layers.size()) {
    // lowX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = lowX; 
     for (unsigned int i=0;i<layers.size()-1;i++) { 
       const HepTransform3D ltransf = HepTranslateX3D(x_array[i]);
       layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers[i]));
       binSteps.push_back(ltransf.getTranslation()[0]+0.5*layers[i]->thickness()-currX+xShift);
       currX = ltransf.getTranslation()[0]+0.5*layers[i]->thickness()+xShift;     
       //std::cout << "binstep:"<<i <<","<<binSteps[i]<<std::endl;
     }
     const HepTransform3D ltransf = HepTranslateX3D(x_array.back());
     layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers.back()));
     binSteps.push_back(compBounds->halflengthZ()-currX+xShift);
     //std::cout << "binstep:"<<layers.size() <<","<<binSteps.back()<<std::endl;
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
  //std::cout << "processing TGC component, number of children volumes:"<< pv->getLogVol()->getName() << "," << pv->getNChildVols() <<std::endl; 
  //printChildren(pv);
  std::vector<const Trk::PlaneLayer*> layers;
  std::vector<double> x_array;
  std::vector<Trk::ExtendedMaterialProperties> x_mat;
  std::vector<double> x_thickness;
  double currX = -100000;
  // while waiting for better suggestion, define a single material layer
  Trk::ExtendedMaterialProperties matTGC(0.,10e8,10e8,0.,0.,0.);
  double minX = tgcBounds->minHalflengthX();
  double maxX = tgcBounds->maxHalflengthX();
  double halfY = tgcBounds->halflengthY();
  double halfZ = tgcBounds->halflengthZ();
  double thickness =2*halfZ;
  //std::cout << "tgc bounds half y:" << tgcBounds->halflengthY() << std::endl; 
  if ( fabs( tgcBounds->halflengthZ() - 35.00) < tol ) {
    if (!m_matTGC01 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matTGC01 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matTGC = *m_matTGC01; 
  } else if ( fabs( tgcBounds->halflengthZ() - 21.85) < tol ) {
    if (!m_matTGC06 ) { 
      double vol = (minX + maxX)*2*halfY*thickness;
      m_matTGC06 = getAveragedLayerMaterial(pv,vol,thickness); 
    }
    matTGC = *m_matTGC06; 
  } else {
    std::cout << "unknown TGC material:" << tgcBounds->halflengthZ()  << std::endl;
  }

  for (unsigned int ic=0; ic < pv->getNChildVols(); ic++) {
    HepTransform3D transfc = pv->getXToChildVol(ic);
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    if ( clv->getName() == "muo::TGCGas" ) {
      double xl = transfc.getTranslation()[0];
      if (x_array.size() == 0 || xl >= x_array.back()) { 
	x_array.push_back(xl);
      } else {
	unsigned int ix = 0; 
	while ( ix < x_array.size() && x_array[ix] < xl  ) {ix++;}
	x_array.insert(x_array.begin()+ix,xl);
      } 
    }
  } 
  //
  // std::cout << " TGC multilayer has " << x_array.size() << " active layers" << std::endl; 
  //
  double activeThick=0.;
  if (x_array.size()==0) {
    x_array.push_back(0.);
    x_thickness.push_back(thickness );
    activeThick = thickness;
  } else if (x_array.size()==1) {
    x_thickness.push_back(2*fmin(x_array[0]+halfZ,halfZ-x_array[0]));
    activeThick += x_thickness.back();
  } else {
    double currX = -halfZ;
    for (unsigned int il=0; il < x_array.size(); il++) {
      if (il<x_array.size()-1) {
	x_thickness.push_back(2*fmin(x_array[il]-currX,0.5*(x_array[il+1]-x_array[il])));
      } else {
	x_thickness.push_back(2*fmin(x_array[il]-currX,halfZ-x_array[il]));
      }
      currX = x_array[il]+0.5*x_thickness.back();
      activeThick += x_thickness.back();
    }
  }
  // rescale material to match the combined thickness of active layers
  //std::cout << "thickness:active,env:" << activeThick <<"," << thickness << std::endl;
  double scale = activeThick/thickness;
  matTGC = Trk::ExtendedMaterialProperties(scale*matTGC.thickness(),scale*matTGC.x0(),scale*matTGC.l0(),
					   matTGC.averageA(),matTGC.averageZ(),matTGC.averageRho()/scale);

  for (unsigned int il=0; il < x_array.size(); il++) {
    scale = x_thickness[il]/activeThick; 
    Trk::ExtendedMaterialProperties xmatTGC(scale*matTGC.thickness(),scale*matTGC.x0(),scale*matTGC.l0(),
					    matTGC.averageA(),matTGC.averageZ(),matTGC.averageRho()/scale);
    x_mat.push_back(xmatTGC);
    //std::cout << "tgc active layers:" << il << "," << x_array[il]<< ","<<x_thickness[il] << std::endl;
  }    

  // create layers
  const Trk::PlaneLayer* layer;
  Trk::OverlapDescriptor* od=0;
  Trk::TrapezoidBounds* tbounds= new Trk::TrapezoidBounds(minX,maxX,halfY); ;
  Trk::SharedObject<const Trk::SurfaceBounds> bounds(tbounds);
  for (unsigned int iloop=0; iloop<x_array.size(); iloop++) {
    HepTransform3D* cTr = new HepTransform3D( HepTranslateX3D(x_array[iloop])*(*transf) ); // this won't work for multiple layers !!! //
    Trk::HomogenousLayerMaterial tgcMaterial(x_mat[iloop]);  
    layer = new Trk::PlaneLayer(cTr,
                                bounds,
                                tgcMaterial,
                                x_thickness[iloop],
                                od );
    layers.push_back(layer);
    // make preliminary identification of active layers
    layer->setLayerType(1);
    //std::cout << "TGC layer built ok"<<std::endl;
  }
  // create the BinnedArray
  //std::cout << "number of Tgc layers:"<<layers.size()<<std::endl;
  std::vector<Trk::SharedObject<const Trk::Layer> > layerOrder;
  std::vector<double> binSteps;
  // 
  double xShift = transf->getTranslation()[0]; 
  double lowX = - halfZ+xShift;
  if (layers.size()) {
     //lowX = layers[0]->transform().getTranslation()[0]-0.5*layers[0]->thickness();
     currX = lowX; 
     for (unsigned int i=0;i<layers.size()-1;i++) { 
       const HepTransform3D ltransf = HepTranslateX3D(x_array[i]);
       layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers[i]));
       binSteps.push_back(ltransf.getTranslation()[0]+0.5*layers[i]->thickness()-currX+xShift);
       currX = ltransf.getTranslation()[0]+0.5*layers[i]->thickness()+xShift;     
     }
     const HepTransform3D ltransf = HepTranslateX3D(x_array.back());
     layerOrder.push_back(Trk::SharedObject<const Trk::Layer>(layers.back()));
     binSteps.push_back( halfZ-currX+xShift);
  }
  Trk::BinUtility* binUtility = new Trk::BinUtility1DX( lowX, new std::vector<double>(binSteps));
  Trk::LayerArray* tgcLayerArray = 0;
  tgcLayerArray = new Trk::NavBinnedArray1D<Trk::Layer>(layerOrder, binUtility, new HepTransform3D());     

  return tgcLayerArray;

} 

double Muon::MuonStationTypeBuilder::decodeX(const GeoShape* sh) const 
{
  double xHalf = 0;

  const GeoTrd*  trd = dynamic_cast<const GeoTrd*> (sh);
  const GeoBox*  box = dynamic_cast<const GeoBox*> (sh);
  const GeoTube* tub = dynamic_cast<const GeoTube*> (sh);
  const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (sh);
  const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion*> (sh);
  const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);

  if (!trd && !box && !tub && !shift && !uni && !sub ) 
    std::cout << "MuonStationTypeBuilder::decodeX : unknown shape type ?" <<sh->type() << std::endl;  

  //  if (trd ) std::cout << "trapezoid dimensions:" << trd->getXHalfLength1()<<"," << trd->getXHalfLength2()<<
  //	      "," << trd->getYHalfLength1()<<"," << trd->getYHalfLength2()<<"," << trd->getZHalfLength() << std::endl; 
  if (trd) xHalf = fmax( trd->getXHalfLength1(), trd->getXHalfLength2() );
  if (box) xHalf = box->getXHalfLength();
  // if (box ) std::cout << "box dimensions:" << box->getXHalfLength()<<"," << box->getYHalfLength()<<"," << box->getZHalfLength() << std::endl; 
  if (tub) xHalf = tub->getRMax();

  if (sub) {
    // be careful to handle properly GeoModel habit of subtracting large volumes from smaller ones
    // std::cout << " decoding subtraction:" << sub->getOpA()->type() << "," << sub->getOpB()->type() << std::endl;
    double xA = decodeX( sub->getOpA() );
    // double xB = decodeX( sub->getOpB() );
    xHalf = xA;
  }
  if (uni) {
    // std::cout << " decoding union:" << uni->getOpA()->type() << "," << uni->getOpB()->type() << std::endl;
    double xA = decodeX( uni->getOpA() );
    double xB = decodeX( uni->getOpB() );
    xHalf = fmax(xA,xB);
  }
  if (shift) {
    // std::cout << " decoding shift:" << shift->getOp()->type() <<"," << shift->getX().getTranslation() << std::endl;
    double xA = decodeX( shift->getOp() );
    double xB = shift->getX().getTranslation()[0]; 
    xHalf = xA + fabs(xB);
  }

  // std::cout << "MuonStationTypeBuilder::decodeX : returns " << xHalf << std::endl;  
  return xHalf;
}
//

std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>*> Muon::MuonStationTypeBuilder::createLayerRepresentation(const Trk::TrackingVolume* trVol) const
{

  const Trk::Layer* layRepr = 0;   
  if (!trVol) return std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>*> (layRepr,0);

  std::vector<const Trk::Layer*>* multi = new std::vector<const Trk::Layer*>;

  // retrieve volume envelope

  const Trk::CuboidVolumeBounds* cubBounds= dynamic_cast<const Trk::CuboidVolumeBounds*> (&(trVol->volumeBounds()));
  const Trk::TrapezoidVolumeBounds* trdBounds= dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(trVol->volumeBounds()));
  const Trk::DoubleTrapezoidVolumeBounds* dtrdBounds= dynamic_cast<const Trk::DoubleTrapezoidVolumeBounds*> (&(trVol->volumeBounds()));

  const Trk::PlaneLayer* layer = 0;

  if (cubBounds) { 
    double thickness = 2*cubBounds->halflengthX();
    double sf        = 4*cubBounds->halflengthZ()*cubBounds->halflengthY();
    //const std::vector<const Trk::Surface*>* surfs = cubBounds->decomposeToSurfaces(HepTransform3D());
    //const Trk::RectangleBounds* rbounds = dynamic_cast<const Trk::RectangleBounds*> (&(*(surfs))[0]->bounds());
    Trk::RectangleBounds* rbounds=new Trk::RectangleBounds(cubBounds->halflengthY(),cubBounds->halflengthZ());
    Trk::SharedObject<const Trk::SurfaceBounds> bounds(rbounds);
    Trk::OverlapDescriptor* od=0;
    Trk::ExtendedMaterialProperties matProp = collectStationMaterial(trVol,sf);
    if (matProp.thickness() > thickness) {
      *m_log << MSG::DEBUG << " thickness of combined station material exceeds station size:" << trVol->volumeName()<<endreq;
    } else if (matProp.thickness()<thickness && matProp.thickness()>0.) {
      //if (matProp.thickness()> 0.)  matProp *= thickness/matProp.thickness();
      double scale = thickness/matProp.thickness();
      matProp = Trk::ExtendedMaterialProperties(scale*matProp.thickness(),scale*matProp.x0(),scale*matProp.l0(),
						matProp.averageA(),matProp.averageZ(),matProp.averageRho()/scale);
    }
    Trk::HomogenousLayerMaterial mat(matProp);  
    layer = new Trk::PlaneLayer(new HepTransform3D(trVol->transform()*HepRotateY3D(90*deg) * HepRotateZ3D(90*deg) ),
                                bounds, mat, thickness, od, 1 );
    //for (size_t i=0; i<surfs->size(); i++) delete (*surfs)[i];
    //delete surfs;
    // multilayers
    if ( trVol->confinedVolumes()) {
      const std::vector<const Trk::TrackingVolume*> vols = trVol->confinedVolumes()->arrayObjects();
      if (vols.size()>1) {
      for (unsigned int i=0;i<vols.size();i++) {
	Trk::ExtendedMaterialProperties matMulti = collectStationMaterial(vols[i],sf);
        multi->push_back(new Trk::PlaneLayer(new HepTransform3D(vols[i]->transform()*HepRotateY3D(90*deg) * HepRotateZ3D(90*deg) ),
                                bounds, Trk::HomogenousLayerMaterial(matMulti), matMulti.thickness(), od, 1 ));
	//std::cout <<"multilayer added:"<< i<<","<<matMulti.thickness()<< std::endl;
      }}      
    }
  }
  if (trdBounds) {
    double thickness = 2*trdBounds->halflengthZ();
    double sf        = 2*(trdBounds->minHalflengthX()+trdBounds->maxHalflengthX())*trdBounds->halflengthY();
    const std::vector<const Trk::Surface*>* surfs = trdBounds->decomposeToSurfaces(HepTransform3D());
    const Trk::TrapezoidBounds* tbounds = dynamic_cast<const Trk::TrapezoidBounds*> (&(*(surfs))[0]->bounds());
    Trk::SharedObject<const Trk::SurfaceBounds> bounds(new Trk::TrapezoidBounds(*tbounds));
    Trk::OverlapDescriptor* od=0;
    Trk::ExtendedMaterialProperties matProp = collectStationMaterial(trVol,sf);
    if (matProp.thickness() > thickness) {
      *m_log << MSG::DEBUG << " thickness of combined station material exceeds station size:" << trVol->volumeName()<<endreq;
    } else if (matProp.thickness()<thickness && matProp.thickness()>0.) {
      //if (matProp.thickness()> 0.)  matProp *= thickness/matProp.thickness();
      double scale = thickness/matProp.thickness();
      matProp = Trk::ExtendedMaterialProperties(scale*matProp.thickness(),scale*matProp.x0(),scale*matProp.l0(),
						matProp.averageA(),matProp.averageZ(),matProp.averageRho()/scale);
    }
    Trk::HomogenousLayerMaterial mat(matProp);  
    layer = new Trk::PlaneLayer(new HepTransform3D(trVol->transform()),
				bounds, mat, thickness, od, 1 );
    for (size_t i=0; i<surfs->size(); i++) delete (*surfs)[i];
    delete surfs;
    // multilayers
    if ( trVol->confinedVolumes() ) {
      const std::vector<const Trk::TrackingVolume*> vols = trVol->confinedVolumes()->arrayObjects();
      if (vols.size()>1) {
      for (unsigned int i=0;i<vols.size();i++) {
	Trk::ExtendedMaterialProperties matMulti = collectStationMaterial(vols[i],sf);
        multi->push_back(new Trk::PlaneLayer(new HepTransform3D(vols[i]->transform()),
					     bounds, Trk::HomogenousLayerMaterial(matMulti), matMulti.thickness(), od, 1 ));
	//std::cout <<"multilayer added:"<< i<<","<<matMulti.thickness()<< std::endl;
      }}      
    }
  }
  if (dtrdBounds) {
    double thickness = 2*dtrdBounds->halflengthZ();
    double sf        = 2*(dtrdBounds->minHalflengthX()+dtrdBounds->medHalflengthX())*dtrdBounds->halflengthY1()
                      +2*(dtrdBounds->medHalflengthX()+dtrdBounds->maxHalflengthX())*dtrdBounds->halflengthY2();
    const std::vector<const Trk::Surface*>* surfs = dtrdBounds->decomposeToSurfaces(HepTransform3D());
    const Trk::DiamondBounds* dbounds = dynamic_cast<const Trk::DiamondBounds*> (& (*(surfs))[0]->bounds());
    Trk::SharedObject<const Trk::SurfaceBounds> bounds(new Trk::DiamondBounds(*dbounds));
    Trk::OverlapDescriptor* od=0;
    Trk::ExtendedMaterialProperties matProp = collectStationMaterial(trVol,sf);
    if (matProp.thickness() > thickness) {
      *m_log << MSG::DEBUG << " thickness of combined station material exceeds station size:" << trVol->volumeName()<<endreq;
    } else if (matProp.thickness()<thickness && matProp.thickness()>0.) {
      //if (matProp.thickness()> 0.)  matProp *= thickness/matProp.thickness();
      double scale = thickness/matProp.thickness();
      matProp = Trk::ExtendedMaterialProperties(scale*matProp.thickness(),scale*matProp.x0(),scale*matProp.l0(),
						matProp.averageA(),matProp.averageZ(),matProp.averageRho()/scale);
    }
    Trk::HomogenousLayerMaterial mat(matProp);  
    layer = new Trk::PlaneLayer(new HepTransform3D(trVol->transform()),
                                bounds, mat, thickness, od, 1 );
    for (size_t i=0; i<surfs->size(); i++) delete (*surfs)[i];
    delete surfs;
    // multilayers
    if ( trVol->confinedVolumes()) {
      const std::vector<const Trk::TrackingVolume*> vols = trVol->confinedVolumes()->arrayObjects();
      if (vols.size()>1) {
      for (unsigned int i=0;i<vols.size();i++) {
	Trk::ExtendedMaterialProperties matMulti = collectStationMaterial(vols[i],sf);
        multi->push_back(new Trk::PlaneLayer(new HepTransform3D(vols[i]->transform()),
					     bounds, Trk::HomogenousLayerMaterial(matMulti), matMulti.thickness(), od, 1 ));
	//std::cout <<"multilayer added:"<< i<<","<<matMulti.thickness()<< std::endl;
      }}      
    }
  }

  //std::cout << "station:"<<trVol->volumeName()<<",thickness:"<<layer->thickness()<<std::endl;
  //std::cout << "station:"<<trVol->volumeName()<<",dInX0:"<<layer->layerMaterialProperties()->fullMaterial(layer->center())->thicknessInX0()<<std::endl;
  //if (trVol->confinedVolumes()) std::cout << "subvolumes:"<<trVol->confinedVolumes()->arrayObjects().size()<< std::endl;
 
  layRepr = layer;

  if (!multi->size()) {delete multi; multi=0;} 
  return std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>*> (layRepr,multi);
}

Trk::ExtendedMaterialProperties Muon::MuonStationTypeBuilder::collectStationMaterial(const Trk::TrackingVolume* vol, double sf) const
{
  Trk::ExtendedMaterialProperties mat(0.,10.e8,10e8,0.,0.,0.);
  // sf is surface of the new layer used to calculate the average 'thickness' of components
  // layers
  if (vol->confinedLayers()){
    const std::vector<const Trk::Layer*> lays = vol->confinedLayers()->arrayObjects();
    for (unsigned il=0; il<lays.size(); il++) {
      const Trk::MaterialProperties* mProp = lays[il]->layerMaterialProperties()->fullMaterial(lays[il]->surfaceRepresentation().center());
      const Trk::ExtendedMaterialProperties* mLay = dynamic_cast<const Trk::ExtendedMaterialProperties*> (mProp);
      if (mLay) {
	double totalD = mat.thickness()+mLay->thickness();
	double f1 = mat.thickness()/totalD; double f2 = mLay->thickness()/totalD;  
	mat = Trk::ExtendedMaterialProperties(totalD, 1./(f1/mat.x0()+f2/mLay->x0()),
					      1./(f1/mat.l0()+f2/mLay->l0()),
					      f1*mat.averageA()+f2*mLay->averageA(),
					      f1*mat.averageZ()+f2*mLay->averageZ(),
					      f1*mat.averageRho()+f2*mLay->averageRho());
      }
    }
  } 
  if (vol->confinedArbitraryLayers()){
    const std::vector<const Trk::Layer*> lays = *(vol->confinedArbitraryLayers());
    for (unsigned il=0; il<lays.size(); il++) {
      const Trk::MaterialProperties* mProp = lays[il]->layerMaterialProperties()->fullMaterial(lays[il]->surfaceRepresentation().center());
      const Trk::ExtendedMaterialProperties* mLay = dynamic_cast<const Trk::ExtendedMaterialProperties*> (mProp);
      // scaling factor
      const Trk::RectangleBounds* rect = dynamic_cast<const Trk::RectangleBounds*> (&(lays[il]->surfaceRepresentation().bounds()));
      if (rect && mLay) {
	double scale = 4*rect->halflengthX()*rect->halflengthY()/sf;
	double totalD = mat.thickness()+scale*mLay->thickness();
	double f1 = mat.thickness()/totalD; double f2 = scale*mLay->thickness()/totalD;  
	mat = Trk::ExtendedMaterialProperties(totalD,1./(f1/mat.x0()+f2/mLay->x0()),
					      1./(f1/mat.l0()+f2/mLay->l0()),
					      f1*mat.averageA()+f2*mLay->averageA(),
					      f1*mat.averageZ()+f2*mLay->averageZ(),
					      f1*mat.averageRho()+f2*mLay->averageRho());
      }
    }
  } 
  // subvolumes
  if (vol->confinedVolumes()){
    const std::vector<const Trk::TrackingVolume*> subVols = vol->confinedVolumes()->arrayObjects();
    for (unsigned iv=0; iv<subVols.size(); iv++) {
      if (subVols[iv]->confinedLayers()){
	const std::vector<const Trk::Layer*> lays = subVols[iv]->confinedLayers()->arrayObjects();
	for (unsigned il=0; il<lays.size(); il++) {
	  const Trk::MaterialProperties* mProp = lays[il]->layerMaterialProperties()->fullMaterial(lays[il]->surfaceRepresentation().center());
	  const Trk::ExtendedMaterialProperties* mLay = dynamic_cast<const Trk::ExtendedMaterialProperties*> (mProp);
          if (mLay) { 
	    double totalD = mat.thickness()+mLay->thickness();
	    double f1 = mat.thickness()/totalD; double f2 = mLay->thickness()/totalD;  
	    mat = Trk::ExtendedMaterialProperties(totalD, f1*mat.x0()+f2*mLay->x0(), f1*mat.l0()+f2*mLay->l0(),
						  f1*mat.averageA()+f2*mLay->averageA(),f1*mat.averageZ()+f2*mLay->averageZ(),
						  f1*mat.averageRho()+f2*mLay->averageRho());
	  }
        }
      } 
      if (subVols[iv]->confinedArbitraryLayers()){
	const std::vector<const Trk::Layer*> lays = *(subVols[iv]->confinedArbitraryLayers());
	for (unsigned il=0; il<lays.size(); il++) {
	  const Trk::MaterialProperties* mProp = lays[il]->layerMaterialProperties()->fullMaterial(lays[il]->surfaceRepresentation().center());
	  const Trk::ExtendedMaterialProperties* mLay = dynamic_cast<const Trk::ExtendedMaterialProperties*> (mProp);
	  // scaling factor
	  const Trk::RectangleBounds* rect = dynamic_cast<const Trk::RectangleBounds*> (&(lays[il]->surfaceRepresentation().bounds()));
	  if (rect && mLay) {
	    double scale = 4*rect->halflengthX()*rect->halflengthY()/sf;
	    double totalD = mat.thickness()+scale*mLay->thickness();
	    double f1 = mat.thickness()/totalD; double f2 = scale*mLay->thickness()/totalD;  
	    mat = Trk::ExtendedMaterialProperties(totalD, f1*mat.x0()+f2*mLay->x0(), f1*mat.l0()+f2*mLay->l0(),
						  f1*mat.averageA()+f2*mLay->averageA(),f1*mat.averageZ()+f2*mLay->averageZ(),
						  f1*mat.averageRho()+f2*mLay->averageRho());
	  }
        }
      }     
    }
  }
  return mat;
}
