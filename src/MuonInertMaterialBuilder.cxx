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
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkDetDescrUtils/BinUtility1DR.h"
#include "TrkDetDescrUtils/BinUtility1DZ.h"
#include "TrkDetDescrUtils/BinnedArray.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkVolumes/SubtractedVolumeBounds.h"
#include "TrkVolumes/CombinedVolumeBounds.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkGeometry/DiscLayer.h"
#include "TrkGeometry/PlaneLayer.h"

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

#include "TrkSurfaces/EllipseBounds.h"

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
	HepTransform3D* transf = new HepTransform3D( top->getXToChildVol(ichild) );
	std::cout <<"MW volume: "<<vname<<" position:" << transf->getTranslation() << std::endl;
	std::vector<const Trk::TrackingVolume*>::const_iterator msTypeIter = msTypes->begin();
	
	for (; msTypeIter != msTypes->end(); ++msTypeIter) {
	  std::string msTypeName = (*msTypeIter)->volumeName();
	  if (msTypeName == vname) {
	    const Trk::TrackingVolume* msTV = *msTypeIter;
	    const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(msTypeName,msTV);
            const Trk::DetachedTrackingVolume* newStat = typeStat->clone(vname,*transf);
	      
	    mInert.push_back(newStat);
	    std::cout<<"MW building inert element: "<<vname<<" position "<<transf->getTranslation()<<std::endl;
	  }
	} // msType
      }
    }
  }
  const std::vector<const Trk::DetachedTrackingVolume*>* muonObjects=new std::vector<const Trk::DetachedTrackingVolume*>(mInert);
//  std::cout << "muonInertMaterialBuilder returns:" << (*muonObjects).size() << std::endl;

//   checkObject(muonObjects);   /debug function to printout object features

  return muonObjects;

}

const std::vector<const Trk::TrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumeTypes() const
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building muon object types" << endreq;
    std::vector<const Trk::TrackingVolume*> objs;

    std::vector<std::string> objName;
    std::vector<unsigned int> objCount;

    if (m_muonMgr){

      // link to top tree
//      std::cout << "number of top tree:" << m_muonMgr->getNumTreeTops() << std::endl;
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));


      for (unsigned int ichild =0; ichild< top->getNChildVols(); ichild++)
      {
        const GeoVPhysVol* cv = &(*(top->getChildVol(ichild)));
        const GeoLogVol* clv = cv->getLogVol();
        std::string vname = clv->getName();
        
        if ((vname.size()>7 && vname.substr(vname.size()-7,7) !="Station") ||  vname.size()<8 )     // Inert element
	{
//	  std::cout << " INERT muon object found:" << vname <<" "<<ichild<< std::endl;

	  bool found = false;
	  for (unsigned int ip=0; ip<objName.size();ip++) {
	    if (vname==objName[ip]) {
	      ++(objCount[ip]);
	      found = true; 
	    } 
	  }  
	  // if (!found && (vname.substr(0,3)!="ECT" || vname=="ECTWallStdSegment")) {
	  if (!found ) {
	    objName.push_back(vname);
	    objCount.push_back(1);
	    printInfo(cv);

	    const Trk::Volume* envelope = translateGeoShape(clv->getShape(),new HepTransform3D());
	    if (envelope) {  
	      Trk::MaterialProperties mat = m_materialConverter->convert( clv->getMaterial() );
	      const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *envelope, mat, m_muonMagneticField,
									   0,0,vname);
	      objs.push_back(newType);
	    }             

	  } // end new object

	}

      }
   }

    // print statistics
    for (unsigned int i=0;i<objName.size();i++) std::cout << "statistics:" << objName[i] << "," << objCount[i] << std::endl; 

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
  Trk::BevelledCylinderVolumeBounds* bevCylBounds=0;

  while ( sh->type() == "Subtraction" ) {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
    sh = sub->getOpA();
  }
  if ( sh->type() == "Tube" ) {
    const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);
//  MW hard-coded angle, to be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double theta = 22.5*deg;
    return new Trk::BevelledCylinderVolumeBounds(tube->getRMin(),tube->getRMax(),
              tube->getZHalfLength(),theta,theta);
  }
  return bevCylBounds;
}


//////////////////////////
/*

  void Muon::MuonInertMaterialBuilder::checkObject( const std::vector<const Trk::DetachedTrackingVolume*>* trkVolumes) const{





       std::vector<const Trk::DetachedTrackingVolume*>::const_iterator msTypeIter = trkVolumes->begin();
//// printout loop
       int ij=0;
       for (; msTypeIter != trkVolumes->end(); ++msTypeIter) {
         ij++;
         std::cout<<"MW InertElement "<<ij<<" "<<(*msTypeIter)->name()<<std::endl;

         if (((*msTypeIter)->name()).substr(0,6)=="BTCold") {
           const Trk::TrapezoidVolumeBounds* tzVolumeBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));

           tzVolumeBounds->dump(std::cout);
           std::cout<<std::endl;
           const HepTransform3D& transf =(*msTypeIter)->trackingVolume()->transform();
	   std::cout<<transf.getTranslation()<<transf.getRotation()<<std::endl;

	   const std::vector<const Trk::Surface*>* tzSurfaces =  tzVolumeBounds->decomposeToSurfaces(*(new HepTransform3D()) );
           std::vector<const Trk::Surface*>::const_iterator Iter = tzSurfaces->begin();
           for (; Iter != tzSurfaces->end(); ++Iter) {
             (*Iter)->dump(std::cout);
	     std::cout<<std::endl;

	   }
	 }


         if (((*msTypeIter)->name()).substr(0,10)=="BTBevelled") {
           const Trk::BevelledCylinderVolumeBounds* bcVolumeBounds = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&((*msTypeIter)->trackingVolume()->volumeBounds()));
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
       myfile.open ("example.txt",std::ios::out);


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
          const HepPoint3D* center2 = new HepPoint3D(xx,yy,zz);

          int is=0;

          std::vector<const Trk::DetachedTrackingVolume*>::const_iterator msTypeIter = trkVolumes->begin();
          for (; msTypeIter != trkVolumes->end(); ++msTypeIter) {


           const HepTransform3D* transf = new HepTransform3D((*msTypeIter)->trackingVolume()->transform().inverse());
           const HepPoint3D* centerHnew = new HepPoint3D((*transf)*(*center2));

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
            if (isInside) {
	       myfile <<xx<<" "<<yy<<" "<<zz<<" "<<is<<std::endl;

           }
           is++;

           if (transf) delete transf;
           if (centerHnew) delete centerHnew;
          }


         }
        }
       }
       myfile.close();

  }
*/

const void Muon::MuonInertMaterialBuilder::printInfo(const GeoVPhysVol* pv) const
{
  const GeoLogVol* lv = pv->getLogVol();
  std::cout << "New Muon Inert Object:"<<lv->getName()<<", made of"<<lv->getMaterial()->getName()<<","<<lv->getShape()->type()<<std::endl;
  decodeShape(lv->getShape());
  printChildren(pv);
}

const void Muon::MuonInertMaterialBuilder::decodeShape(const GeoShape* sh) const
{
  std::cout << "  " ;
  std::cout << "decoding shape:"<< sh->type() << std::endl;

  if ( sh->type()=="Pgon") {
    const GeoPgon* pgon = dynamic_cast<const GeoPgon*>(sh);
    std::cout << "polygon: "<<pgon->getNPlanes()<<" planes "<<pgon->getSPhi()<<" "<<pgon->getDPhi()<<" "<<pgon->getNSides()<<std::endl;
  }

  if ( sh->type()=="Trd") {
    const GeoTrd* trd = dynamic_cast<const GeoTrd*> (sh);
    //
    std::cout << "dimensions:"<< trd->getXHalfLength1() <<","
	      << trd->getXHalfLength2() <<","  
	      << trd->getYHalfLength1() <<","  
	      << trd->getYHalfLength2() <<","  
	      << trd->getZHalfLength() <<std::endl; 
    //
  } 
  if ( sh->type()=="Box") {
    const GeoBox* box = dynamic_cast<const GeoBox*> (sh);
    //
    std::cout << "dimensions:"<< box->getXHalfLength() <<","
	      << box->getYHalfLength() <<","  
	      << box->getZHalfLength() <<std::endl; 
    //
  } 

  if ( sh->type() == "Tube" ) {
    const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);
    std::cout<<"dimensions:"<< tube->getRMin() << ","
                            << tube->getRMax() << ","
                            << tube->getZHalfLength() << std::endl;
  }

  if ( sh->type() == "Tubs" ) {
    const GeoTubs* tubs=dynamic_cast<const GeoTubs*> (sh);
    std::cout<<"dimensions:"<< tubs->getRMin() << ","
                            << tubs->getRMax() << ","
	     << tubs->getZHalfLength() <<"," << tubs->getSPhi() <<"," << tubs->getDPhi()  << std::endl;
  }

  if ( sh->type() == "Cons" ) {
    const GeoCons* cons=dynamic_cast<const GeoCons*> (sh);
    std::cout<<"dimensions:"<< cons->getRMin1() << "," << cons->getRMin2() << ","
                            << cons->getRMax1() << "," << cons->getRMax2() << ","
	     << cons->getDZ() <<"," << cons->getSPhi() <<"," << cons->getDPhi()  << std::endl;
  }

  if ( sh->type()=="Subtraction") {
    while ( sh->type() == "Subtraction" ) {
      const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
      sh = sub->getOpA();
      const GeoShape* shs = sub->getOpB();
      std::cout << "decoding subtracted shape:" << std::endl;
      decodeShape(shs);
    }
    std::cout << "decoding shape A:" << std::endl;
    decodeShape(sh);     
  }

  if ( sh->type()=="Union") {
    const GeoShapeUnion* sub = dynamic_cast<const GeoShapeUnion*> (sh);
    const GeoShape* shA = sub->getOpA();
    const GeoShape* shB = sub->getOpB();
    std::cout << "decoding shape A:" << std::endl;
    decodeShape(shA);
    std::cout << "decoding shape B:" << std::endl;
    decodeShape(shB);    
  }
  if ( sh->type()=="Shift") {
    const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (sh);
    const GeoShape* shA = shift->getOp();
    const HepTransform3D transf = shift->getX();
    std::cout << "shifted by:transl:" <<transf.getTranslation() <<", rot:" 
              << transf[0][0]<<"," << transf[1][1] <<"," << transf[2][2] << std::endl;  
    decodeShape(shA);
  }
}

Trk::Volume* Muon::MuonInertMaterialBuilder::translateGeoShape(const GeoShape* sh, HepTransform3D* transf) const
{
  Trk::Volume* vol=0;
  double tol = 0.1;
 
  std::cout << "  " ;
  std::cout << "translating shape:"<< sh->type() << std::endl;
  if ( sh->type()=="Trap") {
    std::cout << "very general trapezoid, trying to convert into ordinary trapezoid" << std::endl;
    const GeoTrap* trap = dynamic_cast<const GeoTrap*> (sh);
    std::cout << "dimensions anyway:z,Theta,Phi,..." << trap->getZHalfLength() <<","
	      << trap->getTheta() <<","  << trap->getPhi() <<","
	      << trap->getDydzn() <<","  << trap->getDxdyndzn() <<","
	      << trap->getDxdypdzn() <<","  << trap->getAngleydzn() <<","
	      << trap->getDydzp() <<","  << trap->getDxdyndzp() <<","
	      << trap->getDxdypdzp() <<","  << trap->getAngleydzp() <<std::endl;

    Trk::TrapezoidVolumeBounds* volBounds=new Trk::TrapezoidVolumeBounds(trap->getDxdyndzp(),trap->getDxdyndzn(),
									 trap->getDydzn(),trap->getZHalfLength() );
    // vol = new Trk::Volume(new HepTransform3D(*transf * HepRotateZ3D(90*deg) ), volBounds );
    vol = new Trk::Volume(new HepTransform3D(*transf), volBounds );

    return vol;
  }

  if ( sh->type()=="Pgon") {
    std::cout << "Processing polygon:" << std::endl;
    const GeoPgon* pgon = dynamic_cast<const GeoPgon*>(sh);
    std::cout << "polygon: "<<pgon->getNPlanes()<<" planes "<<pgon->getSPhi()<<" "<<pgon->getDPhi()<<" "<<pgon->getNSides()<<std::endl;
    for (unsigned ii=0;ii<pgon->getNPlanes(); ii++) 
      std::cout << "polygon: "<<pgon->getZPlane(ii)<<" "<<pgon->getRMinPlane(ii)<<" "<<pgon->getRMaxPlane(ii)<<std::endl;

    double hlz = 0.5*fabs(pgon->getZPlane(1)-pgon->getZPlane(0));
    double phiH = pgon->getDPhi()/(2.*pgon->getNSides()); 
    double hly = 0.5*cos(phiH)*(pgon->getRMaxPlane(0)-pgon->getRMinPlane(0));
    double dly = 0.5*cos(phiH)*(pgon->getRMaxPlane(0)+pgon->getRMinPlane(0));
    double hlxmin =pgon->getRMinPlane(0)*sin(phiH);
    double hlxmax =pgon->getRMaxPlane(0)*sin(phiH);

    if (pgon->getDPhi()==2*M_PI) {

      /*
      std::cout<<"Polygon as 2-side trapezoid "<<hlxmin<<" "<<hlxmax<<" "<<hly<<" "<<hlz<<std::endl;

      Trk::Volume* volume = 0;
      std::vector<Trk::Volume*> volSeg;
      for (unsigned int i=0; i<pgon->getNSides(); i++) { 
        Trk::TrapezoidVolumeBounds* volBounds = new Trk::TrapezoidVolumeBounds(hlxmin,hlxmax,hly,hlz);
	Trk::Volume* volS = new Trk::Volume(new HepTransform3D(*transf*
		HepRotateZ3D(2*i*phiH)*HepTranslateX3D(dly)),volBounds);
        volSeg.push_back(volS);
	std::cout << "segment center:" <<i << ":" << volS->center() << std::endl;
      }
      for (unsigned int i=0; i<pgon->getNSides(); i++) { 
        if (volume) {
          Trk::CombinedVolumeBounds* combBounds = new Trk::CombinedVolumeBounds(volume,volSeg[i],false);
          volume = new Trk::Volume(new HepTransform3D(*transf),combBounds);
        } else {
          volume = volSeg[i];
        }       
      } 
      */ 
      
      Trk::CylinderVolumeBounds* volBounds = new Trk::CylinderVolumeBounds(pgon->getRMaxPlane(0),hlz);
      Trk::CuboidVolumeBounds* subBounds = new Trk::CuboidVolumeBounds(hlxmax+tol,hlxmax+tol,hlz+tol);
      Trk::Volume* volume = new Trk::Volume(new HepTransform3D(*transf),volBounds);
      for (unsigned int i=0; i<pgon->getNSides(); i++) { 
	Trk::Volume* volS = new Trk::Volume(new HepTransform3D(*transf*
                HepRotateZ3D(2*i*phiH)*HepTranslateX3D(hlxmax+cos(phiH)*pgon->getRMaxPlane(0))),subBounds);
        Trk::SubtractedVolumeBounds* combBounds = new Trk::SubtractedVolumeBounds(volume,volS);
        volume = new Trk::Volume(new HepTransform3D(*transf),combBounds);
      }      
      return volume; 
    }
      

   
 
    if (pgon->getNSides() == 1 ) {

        std::cout<<"Polygon as 1-side trapezoid "<<hlxmin<<" "<<hlxmax<<" "<<hly<<" "<<hlz<<std::endl;

	Trk::TrapezoidVolumeBounds* volBounds = new Trk::TrapezoidVolumeBounds(hlxmin,hlxmax,hlz,hly);
	//Trk::TrapezoidVolumeBounds* volBounds = new Trk::TrapezoidVolumeBounds(hlxmin,hlxmax,hly,hlz);
        std::cout<<"Polygon as 1-side trapezoid "<<hlxmin<<" "<<hlxmax<<" "<<hly<<" "<<hlz<<std::endl;
	vol = new Trk::Volume(new HepTransform3D(*transf*HepRotateZ3D(-90*deg)*HepTranslateY3D(dly)*HepRotateX3D(-90*deg)),
			      volBounds);
        return vol;
    }

       
    if (pgon->getNSides() == 2) {

	Trk::TrapezoidVolumeBounds* volBounds1 = new Trk::TrapezoidVolumeBounds(hlxmin,hlxmax,hlz,hly);
	Trk::TrapezoidVolumeBounds* volBounds2 = new Trk::TrapezoidVolumeBounds(hlxmin,hlxmax,hlz,hly);
        std::cout<<"Polygon as 2-side trapezoid "<<hlxmin<<" "<<hlxmax<<" "<<hly<<" "<<hlz<<std::endl;
	Trk::Volume* vol1 = new Trk::Volume(new HepTransform3D(*transf*HepRotateZ3D(+phiH-90*deg)*HepTranslateY3D(dly)*HepRotateX3D(-90*deg)),volBounds1);
	Trk::Volume* vol2 = new Trk::Volume(new HepTransform3D(*transf*HepRotateZ3D(-phiH-90*deg)*HepTranslateY3D(dly)*HepRotateX3D(-90*deg)),volBounds2);
       
	// Trk::CombinedVolumeBounds* combBounds = new Trk::CombinedVolumeBounds(vol1,vol2,false);
	// vol = new Trk::Volume(new HepTransform3D(), combBounds ); 
	Trk::CylinderVolumeBounds* cylBounds = new Trk::CylinderVolumeBounds(0,dly+hly,hlz);
	vol = new Trk::Volume(new HepTransform3D(), cylBounds ); 

        return vol;
    }
     
    return vol;

  }

  if ( sh->type()=="Trd") {
    const GeoTrd* trd = dynamic_cast<const GeoTrd*> (sh);
    //
    double x1= trd->getXHalfLength1();
    double x2= trd->getXHalfLength2();  
    double y1= trd->getYHalfLength1();  
    double y2= trd->getYHalfLength2();  
    double z = trd->getZHalfLength();  
    //
    if (y1==y2) {
      Trk::TrapezoidVolumeBounds* volBounds=new Trk::TrapezoidVolumeBounds(x1,x2,y1,z);
      //Trk::TrapezoidVolumeBounds* volBounds=new Trk::TrapezoidVolumeBounds(x1,x2,z,y1);
      vol = new Trk::Volume(new HepTransform3D(*transf),volBounds);
      return vol;
    } else if (x1==x2) {
      Trk::TrapezoidVolumeBounds* volBounds=new Trk::TrapezoidVolumeBounds(y1,y2,x1,z);
      vol = new Trk::Volume(new HepTransform3D(*transf*HepRotateZ3D(90*deg)), volBounds );
      return vol;
    } else {
      std::cout << "PROBLEM: translating trapezoid: not recognized:" 
		<< x1 << "," << x2 << "," << y1 << "," << y2 <<"," << z << std::endl; 
    }  

  } 
  if ( sh->type()=="Box") {
    const GeoBox* box = dynamic_cast<const GeoBox*> (sh);
    //
    double x = box->getXHalfLength();
    double y = box->getYHalfLength();  
    double z = box->getZHalfLength(); 
    Trk::CuboidVolumeBounds* volBounds=new Trk::CuboidVolumeBounds(x,y,z);
    vol = new Trk::Volume(new HepTransform3D(*transf), volBounds );
    return vol;
    //
  } 
  if ( sh->type()=="Para") {
    const GeoPara* para = dynamic_cast<const GeoPara*> (sh);
    //
    double x = para->getXHalfLength();
    double y = para->getYHalfLength();  
    double z = para->getZHalfLength(); 
    double alpha = para->getAlpha();
    double theta = para->getTheta();
    double phi   = para->getPhi();
    std::cout << " para:dim:" << x <<"," << y << "," << z << "," << alpha << "," << theta << "," << phi << std::endl;
    Trk::CuboidVolumeBounds* volBounds=new Trk::CuboidVolumeBounds(x,y,z);
    vol = new Trk::Volume(new HepTransform3D(*transf), volBounds );
    return vol;
    //
  } 
  if ( sh->type() == "Tube" ) {
    const GeoTube* tube=dynamic_cast<const GeoTube*> (sh);
    double rMin= tube->getRMin();
    double rMax= tube->getRMax();
    double z =   tube->getZHalfLength();
    Trk::CylinderVolumeBounds* volBounds=new Trk::CylinderVolumeBounds(rMin,rMax,z);
    vol = new Trk::Volume(new HepTransform3D(*transf), volBounds );
    return vol;
   }

  if ( sh->type() == "Tubs" ) {
    const GeoTubs* tubs=dynamic_cast<const GeoTubs*> (sh);
    double rMin= tubs->getRMin();
    double rMax= tubs->getRMax();
    double z =   tubs->getZHalfLength();
    double aPhi =   tubs->getSPhi();
    double dPhi =   tubs->getDPhi();
    Trk::CylinderVolumeBounds* volBounds=new Trk::CylinderVolumeBounds(rMin,rMax,0.5*dPhi,z);
    vol = new Trk::Volume(new HepTransform3D(*transf * HepRotateZ3D(aPhi+dPhi*2)), volBounds );
    /*
    const std::vector<const Trk::Surface*>* surf = vol->volumeBounds().decomposeToSurfaces(vol->transform());
    if (surf->size()==5) {
      std::cout <<"sector plane4:"<<(*surf)[3]->center() <<","<<(*surf)[3]->normal() << std::endl; 
      std::cout <<"sector plane5:"<<(*surf)[4]->center() <<","<<(*surf)[4]->normal() << std::endl; 
    }
    */
    return vol;
  }

  if ( sh->type() == "Cons" ) {
    const GeoCons* cons=dynamic_cast<const GeoCons*> (sh);
    double rMin1= cons->getRMin1();
    double rMin2= cons->getRMin2();
    double rMax1= cons->getRMax1();
    double rMax2= cons->getRMax2();
    double z    = cons->getDZ();
    double aPhi = cons->getSPhi();
    double dPhi = cons->getDPhi();
    // translate into tube with average radius
    if ( dPhi == 2*M_PI ) {
      Trk::CylinderVolumeBounds* volBounds=new Trk::CylinderVolumeBounds(0.5*(rMin1+rMin2),0.5*(rMax1+rMax2),z);
      vol = new Trk::Volume(new HepTransform3D(*transf), volBounds );
      return vol;
    } else {
      Trk::CylinderVolumeBounds* volBounds=new Trk::CylinderVolumeBounds(0.5*(rMin1+rMin2),0.5*(rMax1+rMax2),0.5*dPhi,z);
      vol = new Trk::Volume(new HepTransform3D(*transf * HepRotateZ3D(aPhi+dPhi*2)), volBounds );
      return vol;
    }    
  }

  if ( sh->type()=="Subtraction") {
    const GeoShapeSubtraction* sub = dynamic_cast<const GeoShapeSubtraction*> (sh);
    const GeoShape* shA = sub->getOpA();
    const GeoShape* shB = sub->getOpB();
    Trk::Volume* volA = translateGeoShape(shA, transf);
    Trk::Volume* volB = translateGeoShape(shB, transf);
    if (!volA || !volB ) return vol;
    std::cout << "Subtracting volumes:" << std::endl;   
    Trk::SubtractedVolumeBounds* volBounds = new Trk::SubtractedVolumeBounds(volA, volB);
    //                                                volA->transform().inverse() * volB->transform());
    std::cout << "volume bounds processed" << std::endl;
    vol = new Trk::Volume(new HepTransform3D(), volBounds );
    std::cout << "volume processed" << std::endl;
    return vol;
  }

  if ( sh->type()=="Union") {
    const GeoShapeUnion* uni = dynamic_cast<const GeoShapeUnion*> (sh);
    const GeoShape* shA = uni->getOpA();
    const GeoShape* shB = uni->getOpB();
    Trk::Volume* volA = translateGeoShape(shA, transf);
    Trk::Volume* volB = translateGeoShape(shB, transf);
    if (!volA || !volB ) return vol;
    std::cout << "Combining volumes:" << std::endl;
    Trk::CombinedVolumeBounds* volBounds = new Trk::CombinedVolumeBounds(volA,volB,false);
    //                                           volA->transform().inverse() * volB->transform());
    std::cout << "volume bounds processed" << std::endl;
    // vol = new Trk::Volume(new HepTransform3D(volA->transform()), volBounds );
    vol = new Trk::Volume(new HepTransform3D(), volBounds );
    std::cout << "volume processed" << std::endl;
    return vol;
  }

  if ( sh->type()=="Intersection") {
    const GeoShapeIntersection* intersect = dynamic_cast<const GeoShapeIntersection*> (sh);
    const GeoShape* shA = intersect->getOpA();
    const GeoShape* shB = intersect->getOpB();
    Trk::Volume* volA = translateGeoShape(shA, transf);
    Trk::Volume* volB = translateGeoShape(shB, transf);
    if (!volA || !volB ) return vol;
    std::cout << "Intersecting volumes:" << std::endl;
    Trk::CombinedVolumeBounds* volBounds = new Trk::CombinedVolumeBounds(volA,volB,true);
    //                                           volA->transform().inverse() * volB->transform());
    std::cout << "volume bounds processed" << std::endl;
    //vol = new Trk::Volume(new HepTransform3D(volA->transform()), volBounds );
    vol = new Trk::Volume(new HepTransform3D(), volBounds );
    std::cout << "volume processed" << std::endl;
    return vol;
  }

  if ( sh->type()=="Shift") {
    const GeoShapeShift* shift = dynamic_cast<const GeoShapeShift*> (sh);
    const GeoShape* shA = shift->getOp();
    const HepTransform3D tr = shift->getX();
    std::cout << "Moving volume:" << std::endl;
    Trk::Volume* vol = translateGeoShape(shA, new HepTransform3D(*transf * tr));
    return vol;
  }
  std::cout << "shape not recognized, return 0" << std::endl;
  return vol;
}

const void Muon::MuonInertMaterialBuilder::printChildren(const GeoVPhysVol* pv) const
{
  // subcomponents
  unsigned int nc = pv->getNChildVols();
  for (unsigned int ic=0; ic<nc; ic++) {
    HepTransform3D transf = pv->getXToChildVol(ic);
 
    //
    std::cout << " dumping transform to subcomponent" << std::endl;
    std::cout << transf[0][0]<<"," <<transf[0][1]<<"," <<transf[0][2]<<","<<transf[0][3] << std::endl;
    std::cout << transf[1][0]<<"," <<transf[1][1]<<"," <<transf[1][2]<<","<<transf[1][3] << std::endl;
    std::cout << transf[2][0]<<"," <<transf[2][1]<<"," <<transf[2][2]<<","<<transf[2][3] << std::endl;
    //
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    std::cout << "  ";
    std::cout << "subcomponent:"<<ic<<":"<<clv->getName()<<", made of"<<clv->getMaterial()->getName()<<","<<clv->getShape()->type()<< ","<< transf.getTranslation()<<std::endl;

    decodeShape(clv->getShape()); 	 

    printChildren(cv);
  }  
   
}
