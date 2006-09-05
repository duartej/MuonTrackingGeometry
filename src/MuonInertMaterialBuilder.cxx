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
// temporary
#include "TrkParameters/Perigee.h"

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

    // retrieve muon objects from GeoModel
    if (msTypes->size()) {
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++)
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild)));
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
       if (vname.substr(0,10)=="BTBevelled" ||
	       (vname.substr(0,6)=="BTCold" && vname.substr(0,9)!="BTColdRib")
	       || vname.substr(0,3)=="XXX")     //  not to build endcap ECT
	{
	  HepTransform3D transf = top->getXToChildVol(ichild);
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
	}
      }
    }
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonObjects=new std::vector<const Trk::DetachedTrackingVolume*>(mInert);
  std::cout << "muonInertMaterialBuilder returns:" << (*muonObjects).size() << std::endl;
///////////////////////////////////////////////////////////////////////////////////////
  // checkObject(muonObjects);
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
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));

      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++) 
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild))); 
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
        if ((vname.size()>7 && vname.substr(vname.size()-7,7) !="Station") ||  vname.size()<8 )     // Inert element
	{
	  if (vname.substr(0,10)=="BTBevelled"){

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

	          // is it in objs vector?
		bool inObjs = false;
		for (unsigned int i=0; i< objs.size(); i++) {
		  if (vname == objs[i]->volumeName()) inObjs=true;
		}

                if ( !inObjs) {
                  Trk::TrapezoidVolumeBounds* envBounds = decodeColdSegment(clv->getShape());

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

	          // is it in objs vector?
		bool inObjs = false;
		for (unsigned int i=0; i< objs.size(); i++) {
		  if (vname == objs[i]->volumeName()) inObjs=true;
		}

                if ( !inObjs) {
                  Trk::VolumeBounds* envBounds = decodeECTSegment(clv->getShape());

		    if (envBounds) {
                        HepTransform3D* tTr=0;
                        if (vname.substr(0,21)=="ECTEndplateStdSegment" ){
			 const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (envBounds);
                         tTr = new HepTransform3D( HepTranslateY3D(
			 trd->halflengthY()*(trd->maxHalflengthX()+trd->minHalflengthX())
			 /(trd->maxHalflengthX()-trd->minHalflengthX())
			 )  *HepTranslateY3D(2000.));
			}
			else
			 tTr = new HepTransform3D(HepTranslateY3D(2000.));
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
  const GeoBox* box = 0;
  while ( sh->type() == "Subtraction" ) {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
     if (sub->getOpB()->type()=="Shift"){
      const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*>(sub->getOpB());
       if (shift->getOp()->type()=="Box") {
          box = dynamic_cast<const GeoBox*>(shift->getOp());
       }
     }
    sh = sub->getOpA();
  }
  if ( sh->type() == "Tube" ) {
    const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);

//  MW hard-coded angle, to be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    double theta = 0.3927;
    double theta = 22.5*deg;
    return new Trk::BevelledCylinderVolumeBounds(tube->getRMin(),tube->getRMax(),
             (tube->getZHalfLength()),theta,theta);
  }
  return new Trk::BevelledCylinderVolumeBounds()=0;
}



///////////////////////////////////////////////////////

Trk::VolumeBounds* Muon::MuonInertMaterialBuilder::decodeECTSegment(const GeoShape* sh) const
{
  return 0; 
}



///////////////////////////////////////////////////////

Trk::TrapezoidVolumeBounds* Muon::MuonInertMaterialBuilder::decodeColdSegment(const GeoShape* sh) const
{

  if (sh->type() == "Trd") {
  	const GeoTrd* trapezoid = dynamic_cast<const GeoTrd*> (sh);

        double halfX1 = trapezoid->getXHalfLength1();
        double halfX2 = trapezoid->getXHalfLength2();
        double halfY1 = trapezoid->getYHalfLength1();
        //double halfY2 = trapezoid->getYHalfLength2();
        double halfZ  = trapezoid->getZHalfLength();
        Trk::TrapezoidVolumeBounds* volBounds = new
 	        Trk::TrapezoidVolumeBounds(fmin(halfX1,halfX2),fmax(halfX1,halfX2),halfZ,halfY1);

        return volBounds;
  }

  return new Trk::TrapezoidVolumeBounds()=0;
}


//////////////////////////
  void Muon::MuonInertMaterialBuilder::checkObject( const std::vector<const Trk::DetachedTrackingVolume*>* trkVolumes) const{

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

      int N = 100;
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
