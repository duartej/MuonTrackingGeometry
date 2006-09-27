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
//#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
//#include "TrkGeometry/CylinderLayer.h"
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

// StoreGate
#include "StoreGate/StoreGateSvc.h"

// BField
#include "BFieldAth/MagFieldAthena.h"

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

// STD
#include <map>

#include <iostream>
#include <fstream>
// Gaudi
#include "GaudiKernel/MsgStream.h"

#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoShapeShift.h"
#include "GeoModelKernel/GeoTube.h"
//mw
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoShapeUnion.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoPgon.h"

//#include "TrkVolumes/BevelledCylinderVolumeBounds.h"
//#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkSurfaces/EllipseBounds.h"
#include "TrkEventPrimitives/GlobalPosition.h"



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

// mw
    m_materialConverter= new Trk::GeoMaterialConverter();


    log << MSG::INFO  << name() <<" initialize() successful" << endreq;
    
  return StatusCode::SUCCESS;
}

const std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumes()
 const
{
  std::vector<const Trk::DetachedTrackingVolume*> mInert;

  if (m_muonMgr) {
    // retrieve muon station prototypes from GeoModel
    const std::vector<const Trk::TrackingVolume*>* msTypes = buildDetachedTrackingVolumeTypes();
    std::cout << " obtained " << msTypes->size() << " prototypes" << std::endl;

    // retrieve muon objects from GeoModel
    if (msTypes->size()) {
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++)
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild)));
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
//        if (vname.substr(0,14)=="BTBevelledLong" ||
        if (vname.substr(0,10)=="BTBevelled" ||
	       (vname.substr(0,6)=="BTCold" && vname.substr(0,9)!="BTColdRib")
	       || vname.substr(0,3)=="ECT")     //  ECT
//        if (vname.substr(0,17)=="BTColdLongSegment" || vname.substr(0,18)=="BTColdShortSegment")
//           if (vname.substr(0,6)=="BTCold" && vname.substr(0,9)!="BTColdRib" && vname.substr(0,17)!="BTColdLongSegment" && vname.substr(0,18)!="BTColdShortSegment")
	{
	  HepTransform3D transf = top->getXToChildVol(ichild);
//	  HepTransform3D transf =HepTransform3D();
	  std::cout <<"MW volume: "<<vname<<" position:" << transf.getTranslation() << std::endl;
	  std::vector<const Trk::TrackingVolume*>::const_iterator msTypeIter = msTypes->begin();

	  for (; msTypeIter != msTypes->end(); ++msTypeIter) {
		std::string msTypeName = (*msTypeIter)->volumeName();
		if (msTypeName == vname) {
		  const Trk::TrackingVolume* msTV = *msTypeIter;
		  const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(msTypeName,msTV);
		  const Trk::DetachedTrackingVolume* newStat = typeStat->clone(vname,transf);

		  mInert.push_back(newStat);
		  std::cout<<"MW building inert element: "<<vname<<" position "<<transf.getTranslation()<<std::endl;
		}
	  } // msType
	}
      }
    }
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonObjects=new std::vector<const Trk::DetachedTrackingVolume*>(mInert);
  std::cout << "muonInertMaterialBuilder returns:" << (*muonObjects).size() << std::endl;
///////////////////////////////////////////////////////////////////////////////////////
   checkObject(muonObjects);
///////////////////////////////////////////////////////////////////////////////////////
  return muonObjects;

}

const std::vector<const Trk::TrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumeTypes() const
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building muon object types" << endreq;
///////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<const Trk::TrackingVolume*> objs;


    if (m_muonMgr){

      // link to top tree
      std::cout << "number of top tree:" << m_muonMgr->getNumTreeTops() << std::endl;
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));


      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++)
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild)));
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
        if ((vname.size()>7 && vname.substr(vname.size()-7,7) !="Station") ||  vname.size()<8 )     // Inert element
	{
	  std::cout << " INERT muon object found:" << vname <<" "<<ichild<< std::endl;



	  if (vname.substr(0,10)=="BTBevelled"){
			std::cout << vname<<","<<clv->getShape()->type() <<","<<std::endl;

	          // is it in objs vector?
		bool inObjs = false;
		for (unsigned int i=0; i< objs.size(); i++) {
		  if (vname == objs[i]->volumeName()) inObjs=true;
		}

                if ( !inObjs) {
                  Trk::BevelledCylinderVolumeBounds* envBounds = decodeBevelledCylinder(clv->getShape());

		    if (envBounds) {
			const Trk::Volume* envelope= new Trk::Volume(new HepTransform3D(),envBounds);
			Trk::MaterialProperties newMaterialProp= m_materialConverter->convert( clv->getMaterial() );


			const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *envelope,
										newMaterialProp,
										m_muonMagneticField,
										0,0,  // second: confinedVolumes
										vname);


			objs.push_back(newType);
                  }
	        }
	  }

// Cold segments
	  if (vname.substr(0,6)=="BTCold" && vname.substr(0,9)!="BTColdRib"){
			std::cout << vname<<","<<clv->getShape()->type() <<","<<std::endl;

	          // is it in objs vector?
		bool inObjs = false;
		for (unsigned int i=0; i< objs.size(); i++) {
		  if (vname == objs[i]->volumeName()) inObjs=true;
		}

                if ( !inObjs) {
                  Trk::TrapezoidVolumeBounds* envBounds = decodeColdSegment(clv->getShape());
//                  Trk::CuboidVolumeBounds* envBounds = decodeColdSegment(clv->getShape());

		    if (envBounds) {
		        HepTransform3D* tTr = new HepTransform3D( HepRotateX3D(90*deg) );
			const Trk::Volume* envelope= new Trk::Volume( tTr, envBounds);
			Trk::MaterialProperties newMaterialProp= m_materialConverter->convert( clv->getMaterial() );


			const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *envelope,
										newMaterialProp,
										m_muonMagneticField,
										0,0,  // second: confinedVolumes
										vname);


			objs.push_back(newType);
                  }
	        }
	  }


// ECT (EndcapToroid) segments
	  if (vname.substr(0,3)=="ECT" ){
//	  if (vname.substr(0,7)=="ECTWall" ){
			std::cout << vname<<","<<clv->getShape()->type() <<","<<std::endl;

	          // is it in objs vector?
		bool inObjs = false;
		for (unsigned int i=0; i< objs.size(); i++) {
		  if (vname == objs[i]->volumeName()) inObjs=true;
		}

                if ( !inObjs) {
                  Trk::VolumeBounds* envBounds = decodeECTSegment(clv->getShape());
//Trk::CylinderVolumeBounds* envBounds = new Trk::CylinderVolumeBounds(300.,600.,10000.);

		    if (envBounds) {
		    std::cout<<"MW///////////// "<<vname<<std::endl;
                    envBounds->dump(std::cout);
		    std::cout<<std::endl;
		    std::cout<<"MW/////////////"<<std::endl;
//		        HepTransform3D* tTr = new HepTransform3D( HepRotateX3D(90*deg) );
                        HepTransform3D* tTr;
                        if (vname.substr(0,11)=="ECTEndplate" ||vname.substr(0,7)=="ECTWall" ){
			 const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (envBounds);
                         tTr = new HepTransform3D( HepTranslateY3D(
			 trd->halflengthY()*(trd->maxHalflengthX()+trd->minHalflengthX())
			 /(trd->maxHalflengthX()-trd->minHalflengthX())
			 )  *HepTranslateY3D(0.));
			}
			else
			 tTr = new HepTransform3D(HepTranslateY3D(0.));
			const Trk::Volume* envelope= new Trk::Volume( tTr, envBounds);
			Trk::MaterialProperties newMaterialProp= m_materialConverter->convert( clv->getMaterial() );


			const Trk::TrackingVolume* newType= new Trk::TrackingVolume(
									*envelope,
									newMaterialProp,
									m_muonMagneticField,
									0,0,  // second: confinedVolumes
									vname);


			objs.push_back(newType);
                    }

	        }
	  }



	}


      }
   }

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





Trk::BevelledCylinderVolumeBounds* Muon::MuonInertMaterialBuilder::decodeBevelledCylinder(const GeoShape* sh) const
{

  Trk::BevelledCylinderVolumeBounds* bcyl=0;
  //std::cout << "  decoding shape " << sh->type() << std::endl;
  while ( sh->type() == "Subtraction" ) {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
    sh = sub->getOpA();
  }
  // std::cout << "  found the last shape " << sh->type() << std::endl;
  if ( sh->type() == "Tube" ) {
    const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);
    //std::cout <<"Found tube:" << tube->getRMax() <<"," << tube->getRMin() << "," << tube->getZHalfLength() <<std::endl;
//  MW hard-coded angle, to be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    double theta = 0.3927;
    double theta = 22.5*deg;
    return new Trk::BevelledCylinderVolumeBounds(tube->getRMin(),tube->getRMax(),
              tube->getZHalfLength(),theta,theta);
  }
  return bcyl;
}



///////////////////////////////////////////////////////

Trk::VolumeBounds* Muon::MuonInertMaterialBuilder::decodeECTSegment(const GeoShape* sh) const
{
  std::cout << "MW ECT  decoding shape " << sh->type() << std::endl;

  while ( sh->type() == "Subtraction" || sh->type() == "Union" ) {
    if (sh->type() == "Subtraction") {
      const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
      	std::cout << "MW ECT "<<sub->type()<<","<<sub->getOpA()->type()<<","<<sub->getOpB()->type() << std::endl;
	if (sub->getOpA()->type()=="Pgon") {
		const GeoPgon* box = dynamic_cast<const GeoPgon*>(sub->getOpA());
		std::cout <<"MW ECT Pgon "<<box->getSPhi() <<","<<box->getDPhi()<<","<<box->getNSides()<<","<< box->getNPlanes()<<std::endl;
		unsigned int i = box->getNPlanes();
		for (unsigned int j=0; j<i; j++){
		std::cout <<box->getZPlane(j)<<","<<box->getRMinPlane(j)<<","<<box->getRMaxPlane(j)<<","<<std::endl;
		}
//		return new Trk::CylinderVolumeBounds(40.,50.,500.);
	//          return new Trk::CylinderVolumeBounds(box->getRMinPlane(0),box->getRMaxPlane(0),fabs(box->getZPlane(1)-box->getZPlane(0)));
	}
	if (sub->getOpB()->type()=="Shift"){
	const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*>(sub->getOpB());
	std::cout <<"MW ECT Shift "<<shift->getOp()->type()<<","<<
	(shift->getX()).getTranslation()<<(shift->getX()).getRotation() << std::endl;
	if (shift->getOp()->type()=="Box") {
		const GeoBox* box = dynamic_cast<const GeoBox*>(shift->getOp());
		std::cout <<"MW ECT Box "<<box->getXHalfLength()<<","<<box->getYHalfLength()<<","<<box->getZHalfLength()<<","<< std::endl;
	}
	if (shift->getOp()->type()=="Pgon") {
		const GeoPgon* box = dynamic_cast<const GeoPgon*>(shift->getOp());
		std::cout <<"MW ECT Pgon "<<box->getSPhi() <<","<<box->getDPhi()<<","<<box->getNSides()<<","<< box->getNPlanes()<<std::endl;
		unsigned int i = box->getNPlanes();
		for (unsigned int j=0; j<i; j++){
		std::cout <<box->getZPlane(j)<<","<<box->getRMinPlane(j)<<","<<box->getRMaxPlane(j)<<","<<std::endl;
		}
//		return new Trk::CylinderVolumeBounds(40.,50.,500.);
	//          return new Trk::CylinderVolumeBounds(box->getRMinPlane(0),box->getRMaxPlane(0),fabs(box->getZPlane(1)-box->getZPlane(0)));
	}
	}
	sh = sub->getOpA();
    }
    else {
      const GeoShapeUnion* sub = dynamic_cast<const GeoShapeUnion*> (sh);
      	std::cout << "MW ECT "<<sub->type()<<","<<sub->getOpA()->type()<<","<<sub->getOpB()->type() << std::endl;
	if (sub->getOpA()->type()=="Pgon") {
		const GeoPgon* box = dynamic_cast<const GeoPgon*>(sub->getOpA());
		std::cout <<"MW ECT Pgon "<<box->getSPhi() <<","<<box->getDPhi()<<","<<box->getNSides()<<","<< box->getNPlanes()<<std::endl;
		unsigned int i = box->getNPlanes();
		for (unsigned int j=0; j<i; j++){
		std::cout <<box->getZPlane(j)<<","<<box->getRMinPlane(j)<<","<<box->getRMaxPlane(j)<<","<<std::endl;
		}
//		return new Trk::CylinderVolumeBounds(40.,50.,500.);
	//          return new Trk::CylinderVolumeBounds(box->getRMinPlane(0),box->getRMaxPlane(0),fabs(box->getZPlane(1)-box->getZPlane(0)));
	}
	if (sub->getOpB()->type()=="Shift"){
	const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*>(sub->getOpB());
	std::cout <<"MW ECT Shift "<<shift->getOp()->type()<<","<<
	(shift->getX()).getTranslation()<<(shift->getX()).getRotation() << std::endl;
	if (shift->getOp()->type()=="Box") {
		const GeoBox* box = dynamic_cast<const GeoBox*>(shift->getOp());
		std::cout <<"MW ECT Box "<<box->getXHalfLength()<<","<<box->getYHalfLength()<<","<<box->getZHalfLength()<<","<< std::endl;
	}
	if (shift->getOp()->type()=="Pgon") {
		const GeoPgon* box = dynamic_cast<const GeoPgon*>(shift->getOp());
		std::cout <<"MW ECT Pgon "<<box->getSPhi() <<","<<box->getDPhi()<<","<<box->getNSides()<<","<< box->getNPlanes()<<std::endl;
		unsigned int i = box->getNPlanes();
		for (unsigned int j=0; j<i; j++){
		std::cout <<box->getZPlane(j)<<","<<box->getRMinPlane(j)<<","<<box->getRMaxPlane(j)<<","<<std::endl;
		}
//		return new Trk::CylinderVolumeBounds(40.,50.,500.);
	//          return new Trk::CylinderVolumeBounds(box->getRMinPlane(0),box->getRMaxPlane(0),fabs(box->getZPlane(1)-box->getZPlane(0)));
	}
	}
	sh = sub->getOpA();
    }

  }
  std::cout << "ECT  found the last shape " << sh->type() << std::endl;

  if ( sh->type() == "Trd" ) {
    const GeoTrd* trd=dynamic_cast<const GeoTrd*> (sh);
    std::cout <<"Found Trd:" <<trd->getXHalfLength1()<<" "<< trd->getXHalfLength2()<<" " <<trd->getYHalfLength1()<<" "<< trd->getYHalfLength2()<<" "<< trd->getZHalfLength()<<std::endl;
        double halfX1 = trd->getXHalfLength1();
        double halfX2 = trd->getXHalfLength2();
        double halfY1 = trd->getYHalfLength1();
        //double halfY2 = trd->getYHalfLength2();
        double halfZ  = trd->getZHalfLength();
//        Trk::CuboidVolumeBounds* volBounds = new Trk::CuboidVolumeBounds(halfX1,halfY1,halfZ);
        return new Trk::TrapezoidVolumeBounds(fmin(halfX1,halfX2),fmax(halfX1,halfX2),halfZ,halfY1);
  }
  if ( sh->type() == "Box" ) {
    const GeoBox* box=dynamic_cast<const GeoBox*> (sh);
    std::cout <<"Found Box:" <<box->getXHalfLength()<<" "<<box->getYHalfLength()<<" "<< box->getZHalfLength()<<std::endl;
        double halfX = box->getXHalfLength();
        double halfY = box->getYHalfLength();
        double halfZ = box->getZHalfLength();
//        Trk::CuboidVolumeBounds* volBounds = new Trk::CuboidVolumeBounds(halfX1,halfY1,halfZ);
        return new Trk::CuboidVolumeBounds(halfX,halfY,halfZ);
//        return new Trk::CuboidVolumeBounds(200.,200.,200.);
  }

  if ( sh->type() == "Tube" ) {
    const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);

    std::cout <<"Found tube:" << tube->getRMax() <<"," << tube->getRMin() << "," << tube->getZHalfLength() <<std::endl;
    return new Trk::CylinderVolumeBounds(tube->getRMin(),tube->getRMax(),tube->getZHalfLength());
  }
  if ( sh->type() == "Pgon" ) {
    const GeoPgon* box=dynamic_cast<const GeoPgon*> (sh);

          std::cout <<"MW ECT Pgon "<<box->getSPhi() <<","<<box->getDPhi()<<","<<box->getNSides()<<","<< box->getNPlanes()<<std::endl;
          unsigned int i = box->getNPlanes();
          for (unsigned int j=0; j<i; j++){
            std::cout <<box->getZPlane(j)<<","<<box->getRMinPlane(j)<<","<<box->getRMaxPlane(j)<<","<<std::endl;
          }
//    return new Trk::CylinderVolumeBounds(box->getRMinPlane(0),box->getRMaxPlane(0),fabs(box->getZPlane(1)-box->getZPlane(0)));
    double hlz = fabs(box->getZPlane(1)-box->getZPlane(0))/2.;
    double hly =(box->getRMaxPlane(0)*cos(box->getDPhi()/2.)+box->getRMinPlane(0)*cos(box->getDPhi()/2.))/2.;
    double hlxmin =box->getRMinPlane(0)*sin(box->getDPhi()/2.);
    double hlxmax =box->getRMaxPlane(0)*sin(box->getDPhi()/2.);
    return new Trk::TrapezoidVolumeBounds(hlxmin,hlxmax,hly,hlz);
//    return new Trk::TrapezoidVolumeBounds(300.,300.,500.,500.);

  }

  return new Trk::CylinderVolumeBounds();
//  return new Trk::CylinderVolumeBounds(1000.,1200.,2000.);

}



///////////////////////////////////////////////////////

Trk::TrapezoidVolumeBounds* Muon::MuonInertMaterialBuilder::decodeColdSegment(const GeoShape* sh) const
{
  //std::cout << "  decoding shape " << sh->type() << std::endl;
  Trk::TrapezoidVolumeBounds* trap=0;

  if (sh->type() == "Trd") {
  	const GeoTrd* trapezoid = dynamic_cast<const GeoTrd*> (sh);
        std::cout << trapezoid->getXHalfLength1()<<" "
		  << trapezoid->getXHalfLength2()<<" "
		  << trapezoid->getYHalfLength1()<<" "
		  << trapezoid->getYHalfLength2()<<" "
		  << trapezoid->getZHalfLength()<<" "
		  <<std::endl;

        double halfX1 = trapezoid->getXHalfLength1();
        double halfX2 = trapezoid->getXHalfLength2();
        double halfY1 = trapezoid->getYHalfLength1();
        //double halfY2 = trapezoid->getYHalfLength2();
        double halfZ  = trapezoid->getZHalfLength();
//        Trk::CuboidVolumeBounds* volBounds = new Trk::CuboidVolumeBounds(halfX1,halfY1,halfZ);
        Trk::TrapezoidVolumeBounds* volBounds = new
 	        Trk::TrapezoidVolumeBounds(fmin(halfX1,halfX2),fmax(halfX1,halfX2),halfZ,halfY1);
                std::cout<<" Dump TrapezoidVolume Bounds "<<std::endl;
		volBounds->dump(std::cout);
		std::cout<<std::endl;

//         if (halfX1==halfX2 && halfY1==halfY2) volBounds
// 	 = new Trk::TrapezoidVolumeBounds(halfX1,halfX1,halfY1,halfZ);
//         if (halfX1==halfX2 && halfY1!=halfY2 ) {
//              HepTransform3D transf = transf*HepRotateY3D(90*deg)*HepRotateZ3D(90*deg);
//              volBounds = new
// 	        Trk::TrapezoidVolumeBounds(fmin(halfY1,halfY2),fmax(halfY1,halfY2),halfZ,halfX1);
//         }
//         if (halfX1!=halfX2 && halfY1==halfY2 ) {
//              volBounds = new
// 	        Trk::TrapezoidVolumeBounds(fmin(halfX1,halfX2),fmax(halfX1,halfX2),halfY1,halfZ);
//         }

        return volBounds;
//        double alpha = atan((trapezoid->getXHalfLength1()-trapezoid->getXHalfLength2())/2./trapezoid->getZHalfLength());
//        return new Trk::TrapezoidVolumeBounds (trapezoid->getXHalfLength2(),trapezoid->getXHalfLength1() ,trapezoid->getYHalfLength1() ,trapezoid->getZHalfLength() );
//        return new Trk::TrapezoidVolumeBounds (trapezoid->getXHalfLength2(),trapezoid->getYHalfLength1() ,trapezoid->getZHalfLength(), alpha,alpha   );
  }
//  GeoShape* box = 0;

//   while ( sh->type() == "Subtraction" ) {
//     const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
//     std::cout << sub->type()<<","<<sub->getOpA()->type()<<","<<sub->getOpB()->type() << std::endl;
//      if (sub->getOpB()->type()=="Shift"){
//       const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*>(sub->getOpB());
//       std::cout <<"Shift "<<shift->getOp()->type()<<","<<
//        (shift->getX()).getTranslation()<<(shift->getX()).getRotation() << std::endl;
//        if (shift->getOp()->type()=="Box") {
//           const GeoBox* box = dynamic_cast<const GeoBox*>(shift->getOp());
//           std::cout <<"Box "<<box->getXHalfLength()<<","<<box->getYHalfLength()<<","<<box->getZHalfLength()<<","<< std::endl;
//        }
//      }
//     sh = sub->getOpA();
//   }
//  std::cout << "  found the last shape " << sh->type() << std::endl;
//   if ( sh->type() == "Tube" ) {
//     const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);
//
//     std::cout <<"Found tube:" << tube->getRMax() <<"," << tube->getRMin() << "," << tube->getZHalfLength() <<std::endl;
// //  MW hard-coded angle, to be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// //    double theta = 0.3927;
//     double theta = 22.5*deg;
//     return new Trk::BevelledCylinderVolumeBounds(tube->getRMin(),tube->getRMax(),
//              (tube->getZHalfLength()),theta,theta);
//
//   }

  return trap;
}


//////////////////////////
  void Muon::MuonInertMaterialBuilder::checkObject( const std::vector<const Trk::DetachedTrackingVolume*>* trkVolumes) const{

//       const GlobalPosition* center = trkVolume->center();




       std::vector<const Trk::DetachedTrackingVolume*>::const_iterator msTypeIter = trkVolumes->begin();
//// printout loop
       int ij=0;
       for (; msTypeIter != trkVolumes->end(); ++msTypeIter) {
         ij++;
         std::cout<<"MW InertElement "<<ij<<" "<<(*msTypeIter)->name()<<std::endl;

         if (((*msTypeIter)->name()).substr(0,6)=="BTCold") {
           const Trk::TrapezoidVolumeBounds* tzVolumeBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));
//           const Trk::CuboidVolumeBounds* tzVolumeBounds = dynamic_cast<const Trk::CuboidVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));

           tzVolumeBounds->dump(std::cout);
           std::cout<<std::endl;
           const HepTransform3D& transf =(*msTypeIter)->trackingVolume()->transform();
	   std::cout<<transf.getTranslation()<<transf.getRotation()<<std::endl;

	   const std::vector<const Trk::Surface*>* tzSurfaces =  tzVolumeBounds->decomposeToSurfaces(*(new HepTransform3D()) );
           std::vector<const Trk::Surface*>::const_iterator Iter = tzSurfaces->begin();
           for (; Iter != tzSurfaces->end(); ++Iter) {
             (*Iter)->dump(std::cout);
	     std::cout<<std::endl;
/*	     HepVector3D Normal = ((*Iter)->normal());
             HepTransform3D transf2 = transf;
             HepVector3D  Normal2 = transf2*Normal;
             std::cout<<"Center "<<transf2*((*Iter)->center())<<std::endl;
             std::cout<<"Normal "<<Normal<<" "<<Normal2<<std::endl;*/
	   }
	 }


         if (((*msTypeIter)->name()).substr(0,10)=="BTBevelled") {
           const Trk::BevelledCylinderVolumeBounds* bcVolumeBounds = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));
//           const Trk::CylinderVolumeBounds* bcVolumeBounds = dynamic_cast<const Trk::CylinderVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));
           bcVolumeBounds->dump(std::cout);
           std::cout<<std::endl;
           const HepTransform3D& transf =(*msTypeIter)->trackingVolume()->transform();
	   std::cout<<transf.getTranslation()<<transf.getRotation()<<std::endl;

	   const std::vector<const Trk::Surface*>* bcSurfaces =  bcVolumeBounds->decomposeToSurfaces(*(new HepTransform3D()) );
           std::vector<const Trk::Surface*>::const_iterator Iter = bcSurfaces->begin();
           for (; Iter != bcSurfaces->end(); ++Iter) {
             (*Iter)->dump(std::cout);
	     std::cout<<std::endl;
	     HepVector3D Normal = ((*Iter)->normal());
             HepTransform3D transf2 = transf;
             HepVector3D  Normal2 = transf2*Normal;
             std::cout<<"Center "<<transf2*((*Iter)->center())<<std::endl;
             std::cout<<"Normal "<<Normal<<" "<<Normal2<<std::endl;
	   }
	 }
       }
//// end of printout loop

       std::fstream myfile;
       myfile.open ("example.txt",ios::out);

//        double Xmin= -10000. ;
//        double Xmax= 0. ;
//        double Ymin= -10000. ;
//        double Ymax= 0. ;
//        double Zmin= 11000. ;
//        double Zmax= 13000. ;

       double Xmin= -10000. ;
       double Xmax= -7000. ;
       double Ymin=  2500. ;
       double Ymax=  4500. ;
       double Zmin=  11000. ;
       double Zmax=  13000. ;

      int N = 10;
       for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
         for (int k=0;k<N;k++){
          double xx=Xmin+(double)(i)/(double)(N)*(Xmax-Xmin);
          double yy=Ymin+(double)(j)/(double)(N)*(Ymax-Ymin);
          double zz=Zmin+(double)(k)/(double)(N)*(Zmax-Zmin);
//
          const Trk::GlobalPosition* center2 = new Trk::GlobalPosition(xx,yy,zz);
//          std::cout<<"MW POINT3D "<<xx<<" "<<yy<<" "<<zz<<std::endl;

//          int ss=0;
          int is=0;

          std::vector<const Trk::DetachedTrackingVolume*>::const_iterator msTypeIter = trkVolumes->begin();
          for (; msTypeIter != trkVolumes->end(); ++msTypeIter) {

//           const Trk::BevelledCylinderVolumeBounds* bcVolumeBounds = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));

           const HepTransform3D* transf = new HepTransform3D((*msTypeIter)->trackingVolume()->transform().inverse());
           const HepPoint3D* centerHnew = new HepPoint3D((*transf)*(*center2));
// printout
//               if (i==10 && j==10 && k==10) {
//                   std::cout<<"MW POINT3D "<<xx<<" "<<yy<<" "<<zz<<std::endl;
//                   std::cout<<transf->getTranslation()<<transf->getRotation()<<std::endl;
//                   std::cout<<centerHnew->x()<<" "<<centerHnew->y()<<" "<<centerHnew->z()<<std::endl;
//                }
// end of printout
//  	   int aux = (int)((*bcVolumeBounds).inside(*(centerHnew),0.));


//mw           bool isInside = bcVolumeBounds->inside(*(centerHnew),0.);
//mw           bool isInside = (*msTypeIter)->trackingVolume()->inside(*(centerHnew),0.);
            bool isInside = false;
	    if (((*msTypeIter)->name()).substr(0,10)=="BTBevelled") {
                const Trk::BevelledCylinderVolumeBounds* bcVolumeBounds
		= dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));
                isInside = bcVolumeBounds->inside(*(centerHnew),0.);
            }
            if (((*msTypeIter)->name()).substr(0,6)=="BTCold") {
                const Trk::TrapezoidVolumeBounds* tzVolumeBounds
		= dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));
                isInside = tzVolumeBounds->inside(*(centerHnew),0.);
            }
//               std::cout<<(int)((*bcVolumeBounds).inside(*(centerHnew),1.))<<std::endl;
//           std::stringstream my_stringstream;
//           my_stringstream <<  aux;
//           std::string aux2 = my_stringstream.str();
//           if (aux2 != "0") {
            if (isInside) {
	       myfile <<xx<<" "<<yy<<" "<<zz<<" "<<is<<std::endl;
//               ss++; 
           }
           is++;
/*  printout
             if (ss>0 && aux2 != "0") {
                std::cout<<" MW inside "<<is<<" "<<ss<<" "<<centerHnew->x()<<" "<<centerHnew->y()<<" "<<centerHnew->z()<<" ie "<<xx<<" "<<yy<<" "<<zz<<std::endl;
                bcVolumeBounds->dump(std::cout);
                std::cout<<std::endl;
             }*/
           if (transf) delete transf;
           if (centerHnew) delete centerHnew;
          }
//          if (ss>0) {
//                std::cout<<" MW inside "<<ss<<" "<<centerHnew->x()<<" "<<centerHnew->y()<<" "<<centerHnew->z()<<std::endl;
//                myfile <<xx<<" "<<yy<<" "<<zz<<" "<<js<<std::endl;
//          }

         }
        }
       }
       myfile.close();

  }
