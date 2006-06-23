///////////////////////////////////////////////////////////////////
// MuonInertMaterialBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// Muon
#include "MuonTrackingGeometry/MuonInertMaterialBuilder.h"
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
#include "GeoModelKernel/GeoShapeUnion.h"
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoTube.h"

// constructor
Muon::MuonInertMaterialBuilder::MuonInertMaterialBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgTool(t,n,p),
  m_muonMgrLocation("MuonMgr"),
  m_magFieldTool(0),
  m_magFieldToolName("Trk::MagneticFieldTool"),
  m_magFieldToolInstanceName("ATLAS_TrackingMagFieldTool")
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
}

// destructor
Muon::MuonInertMaterialBuilder::~MuonInertMaterialBuilder()
{}

// Athena standard methods
// initialize
StatusCode Muon::MuonInertMaterialBuilder::initialize()
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

    ds = m_detStore->retrieve(m_muonMgr);

    if (ds.isFailure()) {
        log << MSG::ERROR << "Could not get MuonDetectorManager, no layers for muons will be built. " << endreq;
    }

    log << MSG::INFO << m_muonMgr->geometryVersion() << endreq; 
    
    s = toolSvc()->retrieveTool(m_magFieldToolName, m_magFieldToolInstanceName, m_magFieldTool);
    if (s.isFailure())
    {
      log << MSG::ERROR << "Could not retrieve " << m_magFieldToolName << " from ToolSvc. MagneticField will be 0. " << endreq;
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

const std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumes()
 const
{
  std::vector<const Trk::DetachedTrackingVolume*> mInert;

  if (m_muonMgr) { 
    // retrieve muon objects from GeoModel
    const std::vector<const Trk::TrackingVolume*>* msTypes = buildDetachedTrackingVolumeTypes();
    std::cout << "got " << msTypes->size() << "prototypes" << std::endl;
    if (msTypes->size()) {
      GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++) 
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild))); 
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
        if (vname.substr(0,10)=="BTBevelled")     // 
	{
            HepTransform3D transf = top->getXToChildVol(ichild);
	    std::cout << "position:" << transf.getTranslation() << std::endl;
            std::vector<const Trk::TrackingVolume*>::const_iterator msTypeIter = msTypes->begin();
            for (; msTypeIter != msTypes->end(); ++msTypeIter) { 
              std::string msTypeName = (*msTypeIter)->volumeName();
              if (msTypeName == vname) {
                const Trk::TrackingVolume* msTV = *msTypeIter;
                const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(msTypeName,msTV);
                const Trk::DetachedTrackingVolume* newStat = typeStat->clone(vname,transf);
                mInert.push_back(newStat);
             }
            } // msType
	} //BT tube
      }  // top child loop
    } //msTypes   
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonObjects=new std::vector<const Trk::DetachedTrackingVolume*>(mInert);
  std::cout << "muonInertMaterialBuilder returns:" << (*muonObjects).size() << std::endl;
  return muonObjects; 
}

const std::vector<const Trk::TrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumeTypes() const 
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building muon object types" << endreq;    
///////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<const Trk::TrackingVolume*> objs;

    if (m_muonMgr){
      GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++) 
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild))); 
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
        if (vname.substr(0,10)=="BTBevelled")     // 
	{
          std::cout << " muon object found:" << vname << std::endl;
          // does it exist already ?
          bool done = false;
	  for (unsigned int i=0; i< objs.size(); i++) {
            if (vname == objs[i]->volumeName()) done=true;
          }             
          if (!done) {
             std::cout << ", made of"<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<<std::endl;
             Trk::CylinderVolumeBounds* cyl=decodeBaseTube(clv->getShape());
             if (cyl) {
        	// ready to build the prototype
         	const Trk::Volume* env= new Trk::Volume(new HepTransform3D(),cyl);
	        const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *env,
                                                                       m_muonMaterial,
                                                                       m_muonMagneticField,
                                                                       0,0,
                                                                       vname);         
           
                objs.push_back(newType); 
             }
          }
	}
      }
   }
   
///////////////////////////////////////////////////////////////////////////////////////
   const std::vector<const Trk::TrackingVolume*>* mObjects = new std::vector<const Trk::TrackingVolume*>(objs); 
   return mObjects;  
}

// finalize
StatusCode Muon::MuonInertMaterialBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    return StatusCode::SUCCESS;
}
//
const void Muon::MuonInertMaterialBuilder::decodeShape(const GeoShape* sh) const
{
  std::cout << "  decoding shape " << sh->type() << std::endl;
  if ( sh->type() == "Trd" ) {
    const GeoTrd* trd = dynamic_cast<const GeoTrd*> (sh);
    double halfX1 = trd->getXHalfLength1();
    double halfX2 = trd->getXHalfLength2();
    double halfY1 = trd->getYHalfLength1();
    double halfY2 = trd->getYHalfLength2();
    double halfZ  = trd->getZHalfLength();
    std::cout <<"dimensions:" <<halfX1 << "," <<halfX2 << "," <<halfY1 << "," <<halfY2 << "," <<halfZ << std::endl;
  }
  if ( sh->type() == "Tube" ) {
    const GeoTube* tub=dynamic_cast<const GeoTube*> (sh);
    double rMax = tub->getRMax(); 
    double rMin = tub->getRMin(); 
    double hz   = tub->getZHalfLength(); 
    std::cout <<"dimensions:" << rMax <<"," << rMin << "," << hz <<std::endl;

  }
  if ( sh->type() == "Box" ) {
     const GeoBox* box = dynamic_cast<const GeoBox*> (sh);
     double halfX = box->getXHalfLength();
     double halfY = box->getYHalfLength();
     double halfZ = box->getZHalfLength();
     std::cout <<"dimensions:" << halfX << "," << halfY << "," << halfZ << std::endl;
     
  }
  if ( sh->type() == "Union" ) {
    const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion*> (sh);
    const GeoShape* opA = uni->getOpA();
    const GeoShape* opB = uni->getOpB();
    std::cout << "union of " << opA->type() << " and " << opB->type() << std::endl;
    decodeShape(opA);
    decodeShape(opB);
  }
  if ( sh->type() == "Subtraction" ) {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
    const GeoShape* opA = sub->getOpA();
    const GeoShape* opB = sub->getOpB();
    std::cout << "subtraction of " << opA->type() << " and " << opB->type() << std::endl;
    decodeShape(opA);
    decodeShape(opB);
  }
  if ( sh->type() == "Shift" ) {
    const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (sh);
    const GeoShape* op = shift->getOp();
    const HepTransform3D& transf = shift->getX();
    std::cout << "shift of " << op->type() << " with transform " << transf.getTranslation() << std::endl;
    std::cout << " dumping shift transform" << std::endl;
    std::cout << transf[0][0]<<"," <<transf[0][1]<<"," <<transf[0][2]<<","<<transf[0][3] << std::endl;
    std::cout << transf[1][0]<<"," <<transf[1][1]<<"," <<transf[1][2]<<","<<transf[1][3] << std::endl;
    std::cout << transf[2][0]<<"," <<transf[2][1]<<"," <<transf[2][2]<<","<<transf[2][3] << std::endl;
    decodeShape(op);
  }
}
const Trk::CylinderVolumeBounds* Muon::MuonInertMaterialBuilder::decodeBaseTube(const GeoShape* sh) const
{
  std::cout << "  decoding shape " << sh->type() << std::endl;
  const Trk::CylinderVolumeBounds* cyl=0; 
  while ( sh->type() == "Subtraction" ) {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
    std::cout << sub->getOpA()->type()<<","<<sub->getOpB()->type() << std::endl;
    sh = sub->getOpA();
  }
  std::cout << "  shape A " << sh->type() << std::endl;
  if ( sh->type() == "Tube" ) {
    const GeoTube* tub=dynamic_cast<const GeoTube*> (sh);
    double rMax = tub->getRMax(); 
    double rMin = tub->getRMin(); 
    double hz   = tub->getZHalfLength(); 
    
    std::cout <<"dimensions:" << rMax <<"," << rMin << "," << hz <<std::endl;
    return new Trk::CylinderVolumeBounds(rMax,hz);
  }
  return cyl;
}

const void Muon::MuonInertMaterialBuilder::printChildren(const GeoVPhysVol* pv) const
{
  // subcomponents
  unsigned int nc = pv->getNChildVols();
  std::cout << "Number of child volumes:" << nc << std::endl;
  for (unsigned int ic=0; ic<nc; ic++) {
    HepTransform3D transf = pv->getXToChildVol(ic);
   
    std::cout << " dumping transform to subcomponent" << std::endl;
    std::cout << transf[0][0]<<"," <<transf[0][1]<<"," <<transf[0][2]<<","<<transf[0][3] << std::endl;
    std::cout << transf[1][0]<<"," <<transf[1][1]<<"," <<transf[1][2]<<","<<transf[1][3] << std::endl;
    std::cout << transf[2][0]<<"," <<transf[2][1]<<"," <<transf[2][2]<<","<<transf[2][3] << std::endl;
   
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    std::cout << "  ";
    std::cout << "subcomponent:"<<ic<<":"<<clv->getName()<<", made of"<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<< ","<< transf.getTranslation()<<std::endl;
	 
          if ( clv->getShape()->type()=="Trd") {
	      const GeoTrd* trd = dynamic_cast<const GeoTrd*> (clv->getShape());
	      std::cout << "dimensions:"<< trd->getXHalfLength1() <<","
                                        << trd->getXHalfLength2() <<","  
                                        << trd->getYHalfLength1() <<","  
                                        << trd->getYHalfLength2() <<","  
			                << trd->getZHalfLength() <<std::endl; 
          } 
          if ( clv->getShape()->type()=="Box") {
	      const GeoBox* box = dynamic_cast<const GeoBox*> (clv->getShape());
	        std::cout << "dimensions:"<< box->getXHalfLength() <<","
                                        << box->getYHalfLength() <<","  
			                << box->getZHalfLength() <<std::endl; 
          } 

	  printChildren(cv);
    }  
   

}

