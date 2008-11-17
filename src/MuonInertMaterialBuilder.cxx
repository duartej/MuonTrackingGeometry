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
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkVolumes/SubtractedVolumeBounds.h"
#include "TrkVolumes/CombinedVolumeBounds.h"
#include "TrkVolumes/SimplePolygonBrepVolumeBounds.h"
#include "TrkVolumes/PrismVolumeBounds.h"
#include "TrkVolumes/VolumeExcluder.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkSurfaces/TrapezoidBounds.h"
#include "TrkSurfaces/CylinderSurface.h"
#include "TrkGeometry/DiscLayer.h"
#include "TrkGeometry/PlaneLayer.h"
#include "TrkGeometry/CylinderLayer.h"
#include "TrkGeometry/SubtractedPlaneLayer.h"
#include "TrkGeometry/SubtractedCylinderLayer.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/TrackingGeometry.h"
#include "TrkMagFieldInterfaces/IMagneticFieldTool.h"
#include "TrkMagFieldUtils/MagneticFieldMode.h"
#include "TrkMagFieldUtils/MagneticFieldMap.h"
#include "TrkMagFieldUtils/MagneticFieldMapConstant.h"
#include "TrkMagFieldUtils/MagneticFieldMapGrid3D.h"
#include "TrkMagFieldUtils/MagneticFieldMapSolenoid.h"
#include "TrkGeometry/HomogenousLayerMaterial.h"

// StoreGate
#include "StoreGate/StoreGateSvc.h"

// BField
#include "BFieldAth/MagFieldAthena.h"

//CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

//#include "PathResolver/PathResolver.h"

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
#include "GeoModelKernel/GeoVolumeCursor.h"

#include "TrkSurfaces/EllipseBounds.h"

// constructor
Muon::MuonInertMaterialBuilder::MuonInertMaterialBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgTool(t,n,p),
  Trk::TrackingVolumeManipulator(),
  m_muonMgrLocation("MuonMgr"),
  m_simplify(false),
  m_simplifyToLayers(false),
  m_debugMode(false),
  m_buildBT(true),
  m_buildECT(true),
  m_buildFeets(true),
  m_buildRails(1),
  m_buildShields(true),
  m_magFieldTool("Trk::MagneticFieldTool/AtlasMagneticFieldTool"),
  m_rndmGenSvc("RndmGenSvc","randomGen"),
  m_flatDist(0)
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
  declareProperty("SimplifyGeometry",                 m_simplify);
  declareProperty("SimplifyGeometryToLayers",         m_simplifyToLayers);
  declareProperty("DebugMode",                        m_debugMode);
  declareProperty("BuildBarrelToroids",               m_buildBT);
  declareProperty("BuildEndcapToroids",               m_buildECT);
  declareProperty("BuildFeets",                       m_buildFeets);
  declareProperty("BuildRails",                       m_buildRails);
  declareProperty("BuildShields",                     m_buildShields);
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
    if (s.isFailure() ) log<< MSG::INFO << "failing to initialize?" << endreq;

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


    // Retrieve the magnetic field tool   ----------------------------------------------------    
    if (m_magFieldTool.retrieve().isFailure())
    {
      log << MSG::FATAL << "Failed to retrieve tool " << m_magFieldTool << endreq;
      return StatusCode::FAILURE;
    } else
      log << MSG::INFO << "Retrieved tool " << m_magFieldTool << endreq;



    Trk::MagneticFieldProperties muonMagneticFieldProperties(&(*m_magFieldTool), Trk::RealisticField);
    m_muonMagneticField = muonMagneticFieldProperties;

// mw
    m_materialConverter= new Trk::GeoMaterialConverter();
    m_geoShapeConverter= new Trk::GeoShapeConverter();

    // random number generator
    if ( m_rndmGenSvc.retrieve().isFailure() ){
      log << MSG::FATAL << "Could not retrieve " << m_rndmGenSvc << endreq;
      return StatusCode::FAILURE;
    }
    m_flatDist = new Rndm::Numbers( &(*m_rndmGenSvc), Rndm::Flat(0.,1.) );

    if (m_simplifyToLayers) {
      log << MSG::INFO  << name() <<" option Simplify(Muon)GeometryToLayers no longer maintained " << endreq;
    }

    log << MSG::INFO  << name() <<" initialize() successful" << endreq;
    
  return StatusCode::SUCCESS;
}

const std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumes()
 const
{
  MsgStream log(msgSvc(), name());
  std::vector<const Trk::DetachedTrackingVolume*> mInert;

  // retrieve muon station prototypes from GeoModel
  const std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >* msTypes = buildDetachedTrackingVolumeTypes();
  log << MSG::INFO << name() <<" obtained " << msTypes->size() << " prototypes" << endreq;
  
  std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >::const_iterator msTypeIter = msTypes->begin();
  
  for (; msTypeIter != msTypes->end(); ++msTypeIter) {
    std::string msTypeName = (*msTypeIter).first->name();
    const Trk::DetachedTrackingVolume* msTV = (*msTypeIter).first;
    for (unsigned int it=0 ;it<(*msTypeIter).second.size(); it++)  { 
      HepTransform3D combTr((*msTypeIter).second[it]); 
      const Trk::DetachedTrackingVolume* newStat = msTV->clone(msTypeName,combTr);
      mInert.push_back(newStat);
    }
  }

  // clean up prototypes
  for (unsigned int it = 0; it < msTypes->size(); it++) delete (*msTypes)[it].first;
  delete msTypes;   
 
  const std::vector<const Trk::DetachedTrackingVolume*>* muonObjects=new std::vector<const Trk::DetachedTrackingVolume*>(mInert);

  log << MSG::INFO << name() << " returns  " << (*muonObjects).size() << " objects (detached volumes)" << endreq;

  return muonObjects;

}

const std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumeTypes() const
{
    MsgStream log( msgSvc(), name() );

    log << MSG::INFO  << name() <<" building muon object types" << endreq;
    std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > > objs;

    std::vector<std::string> objName;
    
    if (m_muonMgr){
      
      // link to top tree
      const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
      GeoVolumeCursor vol (top);
      while (!vol.atEnd()) {
	const GeoVPhysVol* cv = &(*(vol.getVolume()));
	const GeoLogVol* clv = cv->getLogVol();
	std::string vname = clv->getName();
	
	if ((vname.size()>7 && vname.substr(vname.size()-7,7) !="Station") ||  vname.size()<8 ) {        // Inert element
	  
	  //std::cout << " INERT muon object found:" << vname <<" "<<ichild<< std::endl;
	  
	  bool accepted = true ;
	  if (  vname.substr(0,3)=="BAR" || vname.substr(0,2)=="BT" || vname.substr(0,6) == "EdgeBT" 
		|| vname.substr(0,6) == "HeadBT" ) accepted = m_buildBT ? true : false;
	  else if ( vname.substr(0,3)=="ECT" ) accepted = m_buildECT ? true : false; 
	  else if ( vname.substr(0,4)=="Feet" || ( vname.size()>7 && 
			 (vname.substr(3,4)=="Feet" || vname.substr(4,4)=="Feet" ) ) ) accepted = m_buildFeets ? true : false; 
	  else if ( vname.substr(0,4)=="Rail" ) accepted = m_buildRails>0 ? true : false; 
	  //else accepted = m_buildShields ? true : false;

          if ( vname=="EdgeBTVoussoir" && accepted && m_simplify ) accepted = false;
	  
	  if (!accepted) { vol.next(); continue; }  

          // update to accomodate AGDD structures

	  //printInfo(cv);
	    
	  std::vector<const GeoShape*> input_shapes;
	  std::vector<std::pair<const GeoLogVol*,std::vector<HepTransform3D> > > vols;

          bool simpleTree = false;
          if ( !cv->getNChildVols() ) {
	    std::vector<HepTransform3D > volTr;
	    volTr.push_back(vol.getTransform()); 
	    vols.push_back(std::pair<const GeoLogVol*,std::vector<HepTransform3D> > (clv,volTr) );
            simpleTree = true;
          } else {
	    getObjsForTranslation(cv,HepTransform3D(),vols );
	  }
	  input_shapes.resize(vols.size());             
	  for (unsigned int i=0;i<vols.size();i++) input_shapes[i]=vols[i].first->getShape();

	  for (unsigned int ish=0; ish < vols.size(); ish++) { 
	    
	    std::string protoName = vname;
	    if (!simpleTree) protoName = vname+(vols[ish].first->getName());
	    
	    bool found = false;
	    for (unsigned int ip=0; ip<objs.size();ip++) {
	      if (protoName==objs[ip].first->name()) {
		found = true;
                if (simpleTree) objs[ip].second.push_back(vol.getTransform());
                //else objs[ip].second.insert(objs[ip].second.end(),vols[ish].second.begin(),vols[ish].second.end());
	      } 
	    }  
	    if (found) continue;
            // m_geoShapeConverter->decodeShape(input_shapes[ish]);
            HepTransform3D ident;
	    const Trk::Volume* trObject = m_geoShapeConverter->translateGeoShape(input_shapes[ish],&ident);
	    if (trObject) {  
	      Trk::MaterialProperties mat = m_materialConverter->convert( vols[ish].first->getMaterial() );
	      const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *trObject, mat, m_muonMagneticField,0,0,protoName);
	      const Trk::TrackingVolume* simType = simplifyShape(newType);
	      const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(protoName,simType);
              typeStat->saveConstituents(m_constituents.back());
	      objs.push_back(std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> >(typeStat,vols[ish].second));
              delete trObject;
	    }  else {
	      log << MSG::WARNING << name()<< " volume not translated: " << vname << std::endl;
	    }            
	  } // end new object
	}
	vol.next();      	
      }
    }

    const std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >* mObjects = new std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >(objs);

    int count = 0;
    //double massEstimate = 0.;
    for (unsigned int i=0;i<mObjects->size();i++) { 
       log << MSG::DEBUG << i<< ":"<< (*mObjects)[i].first->name() << ","<< (*mObjects)[i].second.size()<< endreq;
       const Trk::MaterialProperties* mat = (*mObjects)[i].first->trackingVolume()->confinedDenseVolumes() ?
	 (*(*mObjects)[i].first->trackingVolume()->confinedDenseVolumes())[0] :(*mObjects)[i].first->trackingVolume();
       /*
	 for (unsigned int ic = 0; ic<(*mObjects)[i].first->constituents()->size(); ic++) {         
	   double protMass = calculateVolume((*(*mObjects)[i].first->constituents())[ic].first)
	   *(*(*mObjects)[i].first->constituents())[ic].second.first*(*mat).zOverAtimesRho();
           massEstimate += protMass*(*mObjects)[i].second.size();
         }
       */
       for (unsigned int j=0;j<(*mObjects)[i].second.size();j++) 
	 log << MSG::DEBUG << j<< "th  position at "<<((*mObjects)[i].second)[j].getTranslation()<<","<<((*mObjects)[i].second)[j].getRotation()<< endreq;      
       count += (*mObjects)[i].second.size();
    }

    log << MSG::INFO << name() << " returns " << mObjects->size() << " prototypes, to be cloned into "<< count <<" objects" << endreq;   

    return mObjects;
}

// finalize
StatusCode Muon::MuonInertMaterialBuilder::finalize()
{
    MsgStream log(msgSvc(), name());
    delete m_materialConverter;
    delete m_geoShapeConverter;
    delete m_flatDist;
    for (unsigned int i=0;i<m_constituents.size();i++) delete m_constituents[i];
    log << MSG::INFO  << name() <<" finalize() successful" << endreq;
    return StatusCode::SUCCESS;
}
//

const void Muon::MuonInertMaterialBuilder::printInfo(const GeoVPhysVol* pv) const
{
  const GeoLogVol* lv = pv->getLogVol();
  std::cout << "New Muon Inert Object:"<<lv->getName()<<", made of"<<lv->getMaterial()->getName()<<","<<lv->getShape()->type()<<std::endl;
  m_geoShapeConverter->decodeShape(lv->getShape());
  printChildren(pv);
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

    m_geoShapeConverter->decodeShape(clv->getShape()); 	 

    printChildren(cv);
  }  
   
}

const Trk::TrackingVolume* Muon::MuonInertMaterialBuilder::simplifyShape(const Trk::TrackingVolume* trVol) const
{
  // envelope
  const Trk::Volume* envelope = 0;
  // resolve composed volumes (returns constituents with material fraction accounting for subtractions)
  std::vector<std::pair<const Trk::Volume*,double> > constituents = splitComposedVolume(trVol,true);

  //std::cout << "simplifying shape for:"<<trVol->volumeName()<<","<< constituents[0].second<< std::endl;
  //for (unsigned int i=0;i<constituents.size(); i++) std::cout << "constituent:"<< calculateVolume(constituents[i].first)<<","<<
  //						       constituents[i].second<< std::endl;   
  
  int simpleMode = 0;
  
  if (constituents.size() == 1) {  // easy case
    //const Trk::SimplePolygonBrepVolumeBounds* spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*> (&(constituents[0].first->volumeBounds()));
    //envelope = spb ? new Trk::Volume(*(spb->envelope()),trVol->transform()) 
    //  : new Trk::Volume(*(constituents.front().first), trVol->transform());
    envelope =  new Trk::Volume(*(constituents.front().first), trVol->transform());
    if ( constituents.front().second > 0.999 ) simpleMode=1;          // no need to simplify nor envelope
  } else { // construct envelope using constituent edges 
    envelope = createEnvelope(trVol->transform(),constituents);
  }
  if (!envelope) envelope = new Trk::Volume(*trVol);
  
  // simplification
  
  const Trk::TrackingVolume* newVol = 0;
  std::vector<const Trk::TrackingVolume*>* confinedVols = new std::vector<const Trk::TrackingVolume*>;

  std::string envName=trVol->volumeName();
  
  if ( simpleMode == 1 ) {
    //std::cout << "simple mode 1, no envelope for object:"<< trVol->volumeName() << std::endl;
    newVol = trVol;
    delete confinedVols;    
  } else if (m_simplify || envName.substr(0,3)=="BAR" ) { 
    if (constituents.size()==1) {   // simplified volume
      double fraction = constituents.front().second;
      Trk::MaterialProperties mat(1.,trVol->x0()/fraction,fraction*trVol->zOverAtimesRho());
      newVol = new Trk::TrackingVolume( *envelope, mat, m_muonMagneticField, 0, 0, envName);  
      //std::cout << "dense envelope for object:"<< trVol->volumeName() << std::endl;
      delete trVol;  
      delete confinedVols;
    } else {  // enclose simplified constituents
      for (unsigned int ic=0;ic<constituents.size();ic++) { 
	double fraction = constituents[ic].second;
	Trk::MaterialProperties mat(1.,trVol->x0()/fraction,fraction*trVol->zOverAtimesRho());
	Trk::TrackingVolume* trc = new Trk::TrackingVolume(*(constituents[ic].first),mat,m_muonMagneticField, 0, 0, trVol->volumeName());
        confinedVols->push_back(trc);
      }  
      envName=trVol->volumeName()+"_envelope";
      newVol = new Trk::TrackingVolume( *envelope, m_muonMaterial, m_muonMagneticField, confinedVols, envName);    
      for (unsigned int iv = 0; iv < confinedVols->size(); iv++)
	Trk::TrackingVolumeManipulator::confineVolume(*((*confinedVols)[iv]),newVol);
      delete trVol;  
    }
  } else {    // enclose the exact transcript
    confinedVols->push_back(trVol);
    envName=trVol->volumeName()+"_envelope";
    newVol = new Trk::TrackingVolume( *envelope, m_muonMaterial, m_muonMagneticField, confinedVols, envName);    
    Trk::TrackingVolumeManipulator::confineVolume(*trVol,newVol);
    //std::cout << "enclosing object:"<< trVol->volumeName() << std::endl;
  }
  
  
  // save calculable volumes for blending
  std::vector<std::pair<const Trk::Volume*,std::pair<double,double> > > confinedConst;
  for (unsigned int ic=0;ic<constituents.size();ic++) {
    confinedConst.push_back(std::pair<const Trk::Volume*,std::pair<double,double> >
			    ( new Trk::Volume(*(constituents[ic].first),newVol->transform().inverse()),
			      std::pair<double,double>(constituents[ic].second,-1.) ) );
  }
  m_constituents.push_back(new std::vector<std::pair<const Trk::Volume*,std::pair<double,double> > >(confinedConst));
  //for (unsigned int ic=0;ic<constituents.size();ic++) delete constituents[ic].first; 

  delete envelope;
  return newVol;
}
 
double  Muon::MuonInertMaterialBuilder::calculateVolume( const Trk::Volume* envelope) const
{
  double envVol = 0.;
  
  if (!envelope) return 0.;
  
  const Trk::CylinderVolumeBounds*  cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::CuboidVolumeBounds*    box = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::BevelledCylinderVolumeBounds*  bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::PrismVolumeBounds* prism = dynamic_cast<const Trk::PrismVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::SimplePolygonBrepVolumeBounds* spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::CombinedVolumeBounds*  comb = dynamic_cast<const Trk::CombinedVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::SubtractedVolumeBounds*  sub = dynamic_cast<const Trk::SubtractedVolumeBounds*> (&(envelope->volumeBounds()));  

  if ( cyl ) envVol = 2*cyl->halfPhiSector()*(cyl->outerRadius()*cyl->outerRadius()-cyl->innerRadius()*cyl->innerRadius())*cyl->halflengthZ();
  if ( box ) envVol = (8*box->halflengthX()*box->halflengthY()*box->halflengthZ());
  if ( trd ) envVol = (4*(trd->minHalflengthX()+trd->maxHalflengthX())*trd->halflengthY()*trd->halflengthZ());
  if ( bcyl ) {
    int type = bcyl->type();
    if ( type<1 ) envVol = 2*bcyl->halfPhiSector()*(bcyl->outerRadius()*bcyl->outerRadius()-bcyl->innerRadius()*bcyl->innerRadius())*bcyl->halflengthZ(); 
    if ( type==1 ) envVol = 2*bcyl->halflengthZ()*( bcyl->halfPhiSector()*bcyl->outerRadius()*bcyl->outerRadius()
						    -bcyl->innerRadius()*bcyl->innerRadius()*tan(bcyl->halfPhiSector()) ); 
    if ( type==2 ) envVol = 2*bcyl->halflengthZ()*( -bcyl->halfPhiSector()*bcyl->innerRadius()*bcyl->innerRadius()
						    +bcyl->outerRadius()*bcyl->outerRadius()*tan(bcyl->halfPhiSector()) ); 
    if ( type==3 ) envVol = 2*bcyl->halflengthZ()*tan(bcyl->halfPhiSector())*( bcyl->outerRadius()*bcyl->outerRadius() 
									       -bcyl->innerRadius()*bcyl->innerRadius()); 
  }
  if ( prism ) {
    std::vector<std::pair<double,double> > v=prism->xyVertices();
    double a2 = v[1].first*v[1].first+v[1].second*v[1].second
               +v[0].first*v[0].first+v[0].second*v[0].second
	    -2*(v[0].first*v[1].first+v[0].second*v[1].second);
    double c2 = v[2].first*v[2].first+v[2].second*v[2].second
               +v[0].first*v[0].first+v[0].second*v[0].second
	    -2*(v[0].first*v[2].first+v[0].second*v[2].second);
    double ca = v[1].first*v[2].first+v[1].second*v[2].second
               +v[0].first*v[0].first+v[0].second*v[0].second
	       -v[0].first*v[1].first-v[0].second*v[1].second
	       -v[0].first*v[2].first-v[0].second*v[2].second;
    double vv = sqrt(c2-ca*ca/a2);
    envVol = vv*sqrt(a2)*prism->halflengthZ();
  }
  if ( spb ) {
    envVol = calculateVolume(spb->combinedVolume());    // exceptional use of combined volume (no intersections)
  }   
  if ( comb ) {
    envVol = calculateVolume(comb->first()) + calculateVolume(comb->second());
  }
  if ( sub ) {
    return -1;
  }

  return envVol;
}

void Muon::MuonInertMaterialBuilder::getObjsForTranslation(const GeoVPhysVol* pv,HepTransform3D transform, std::vector<std::pair<const GeoLogVol*, std::vector<HepTransform3D> > >& vols ) const
{
  // subcomponents 
  unsigned int nc = pv->getNChildVols();
  //std::cout << "getObjsForTranslation from:"<< pv->getLogVol()->getName()<<","<<pv->getLogVol()->getMaterial()->getName()<<", looping over "<< nc << " children" << std::endl;
  for (unsigned int ic=0; ic<nc; ic++) {
    HepTransform3D transf = pv->getXToChildVol(ic);
    const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
    const GeoLogVol* clv = cv->getLogVol();
    if (!cv->getNChildVols()) {
      bool found = false;
      for (unsigned int is = 0; is < vols.size(); is++) {
	if (clv->getName() == vols[is].first->getName()) {
	  found = true; 
	  vols[is].second.push_back(transform*transf);
	  break;
	}
      }
      if (!found) {
	std::vector<HepTransform3D > volTr;
	volTr.push_back(transform*transf); 
	vols.push_back(std::pair<const GeoLogVol*,std::vector<HepTransform3D> > (clv,volTr) );
	//std::cout << "new volume added:"<< clv->getName() <<","<<clv->getMaterial()->getName()<<std::endl;
	//printInfo(cv);
      }
    } else {
      getObjsForTranslation(cv, transform*transf, vols);
    }
  }
}

const Trk::Volume* Muon::MuonInertMaterialBuilder::createEnvelope(const HepTransform3D transf, std::vector<std::pair<const Trk::Volume*,double> > constituents ) const 
{
  Trk::Volume* envelope = 0;
    
  std::vector<std::pair<const Trk::Volume*,double> >::iterator sIter = constituents.begin();
  std::vector<Trk::GlobalPosition> edges;
  bool cylinder = false;

  double cVol = 0.;
  while (sIter!= constituents.end()) {
    cVol +=(*sIter).second * calculateVolume((*sIter).first);
    const Trk::SimplePolygonBrepVolumeBounds* spbBounds = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*> (&((*sIter).first->volumeBounds()));
    const Trk::CylinderVolumeBounds*          cylBounds = dynamic_cast<const Trk::CylinderVolumeBounds*>  (&((*sIter).first->volumeBounds()));
    const Trk::CuboidVolumeBounds*            cubBounds = dynamic_cast<const Trk::CuboidVolumeBounds*>    (&((*sIter).first->volumeBounds()));
    const Trk::TrapezoidVolumeBounds*         trdBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&((*sIter).first->volumeBounds()));
 
    if (cylBounds) {  
      double rOut = cylBounds->outerRadius();
      double hZ   = cylBounds->halflengthZ();
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( rOut, rOut, hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-rOut, rOut, hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( rOut,-rOut, hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-rOut,-rOut, hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( rOut, rOut,-hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-rOut, rOut,-hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( rOut,-rOut,-hZ));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-rOut,-rOut,-hZ));
      cylinder = true;
    } else if (cubBounds) {  
      double x = cubBounds->halflengthX();
      double y = cubBounds->halflengthY();
      double z = cubBounds->halflengthZ();
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x, y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x, y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x,-y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x,-y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x, y,-z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x, y,-z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x,-y,-z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x,-y,-z));
    } else if (trdBounds) {  
      double x1 = trdBounds->minHalflengthX();
      double x2 = trdBounds->maxHalflengthX();
      double y  = trdBounds->halflengthY();
      double z  = trdBounds->halflengthZ();
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x2, y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x2, y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x1,-y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x1,-y, z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x2, y,-z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x2, y,-z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition( x1,-y,-z));
      edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(-x1,-y,-z));
    } else if (spbBounds) {
      std::vector<std::pair<double,double> > xyVtx = spbBounds->xyVertices();
      double z  = spbBounds->halflengthZ();
      for (unsigned int iv=0;iv<xyVtx.size();iv++) {
	edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(xyVtx[iv].first,xyVtx[iv].second,-z));
	edges.push_back(transf.inverse()*(*sIter).first->transform()*Trk::GlobalPosition(xyVtx[iv].first,xyVtx[iv].second, z));
      }
    } else {
      std::cout << "bounds not recognized"<< std::endl;
      return 0;
    }
    sIter++;
  }
  
  double xMin = 13000.;
  double xMax =-13000.;
  double yMin = 13000.;
  double yMax =-13000.;
  double zMin = 25000.;
  double zMax =-25000.;
  //double rSize = -1.;
  for (unsigned int ie=0; ie < edges.size(); ie++) {
    xMin = fmin( xMin, edges[ie].x() );
    xMax = fmax( xMax, edges[ie].x() );
    yMin = fmin( yMin, edges[ie].y() );
    yMax = fmax( yMax, edges[ie].y() );
    zMin = fmin( zMin, edges[ie].z() );
    zMax = fmax( zMax, edges[ie].z() );
  }

  double xSize = 0.5*(xMax-xMin);
  double ySize = 0.5*(yMax-yMin);
  double zSize = 0.5*(zMax-zMin);

  //std::cout << "envelope parameters:"<< xSize<<","<<ySize<<","<<zSize<< std::endl;
  //std::cout << "envelope position:"<<0.5*(xMin+xMax)<<","<<0.5*(yMin+yMax)<<","<<0.5*(zMin+zMax) << std::endl;
  
  if ( cylinder && fabs(xSize-ySize)/fmax(xSize,ySize)<0.1) { // make it a cylinder
    envelope = new Trk::Volume(new HepTransform3D(transf*HepTranslate3D(Trk::GlobalPosition(0.5*(xMin+xMax),0.5*(yMin+yMax),0.5*(zMin+zMax)))),
			     new Trk::CylinderVolumeBounds(fmax(xSize,ySize),zSize));
  } else {
    envelope = new Trk::Volume(new HepTransform3D(transf*HepTranslate3D(Trk::GlobalPosition(0.5*(xMin+xMax),0.5*(yMin+yMax),0.5*(zMin+zMax)))),
			     new Trk::CuboidVolumeBounds(xSize,ySize,zSize));
  }  
  /*
  // check if all edges really confined:
  for (unsigned int ie=0; ie < edges.size(); ie++) {
    if (!envelope->inside(edges[ie],0.001)) {
      std::cout << "something's wrong with the envelope, edge not confined"<< std::endl;    
      envelope = 0;
    }  
  }
  */
  //std::cout << "volume/envelope fraction:"<< cVol/calculateVolume(envelope)<< std::endl;
  
  return envelope;
}

Trk::GlobalPosition Muon::MuonInertMaterialBuilder::getScanPoint(const Trk::Volume* vol) const
{
  Trk::GlobalPosition gp(0.,0.,0.);
  
  const Trk::CuboidVolumeBounds*     cub = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::CylinderVolumeBounds*   cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::TrapezoidVolumeBounds*  trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::SimplePolygonBrepVolumeBounds*  spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*> (&(vol->volumeBounds()));
  
  if ( !cub && !cyl && !trd && !spb ) return (vol->transform()*gp);
  
  std::vector<double> rndm(3);  
  // generate random numbers
  for (unsigned int ir=0;ir<3;ir++) rndm[ir]= m_flatDist->shoot(); 
  
  if (cub) {
    double x = cub->halflengthX();
    double y = cub->halflengthY();
    double z = cub->halflengthZ();
    gp = Trk::GlobalPosition(-x+2*x*rndm[0],-y+2*y*rndm[1],-z+2*z*rndm[2]);
  } else if (trd) {
    double x1 = trd->minHalflengthX();
    double x2 = trd->maxHalflengthX();
    double y = trd->halflengthY();
    double z = trd->halflengthZ();
    gp = Trk::GlobalPosition(-x2+(x1+x2)*rndm[0],-y+2*y*rndm[1],-z+2*z*rndm[2]);
    if (!vol->inside( vol->transform()*gp,0.001)) {
      gp = Trk::GlobalPosition(x1+(x1+x2)*rndm[0],y-2*y*rndm[1],-z+2*z*rndm[2]);
      if (!vol->inside( vol->transform()*gp,0.001))
	std::cout << "trapezoid hit not correct:"<<-y+2*y*rndm[1]<< std::endl; 
    }
  } else if (cyl) {
    double x1 = cyl->innerRadius();
    double x2 = cyl->outerRadius();
    double y = cyl->halfPhiSector();
    double z = cyl->halflengthZ();   
    double r = sqrt(x1*x1+rndm[0]*(x2*x2-x1*x1));
    double phi = -y+2*y*rndm[1];
    gp = Trk::GlobalPosition(r*cos(phi),r*sin(phi),-z+2*z*rndm[2]);
  } else if (spb) {
    double z = spb->halflengthZ();
    std::vector<const Trk::Volume*> subVols;
    std::vector<double> subVolumes;
    const Trk::Volume* comVol = spb->combinedVolume();
    const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*> (&(comVol->volumeBounds()));
    while (comb) {
      const Trk::CombinedVolumeBounds* comb1 = dynamic_cast<const Trk::CombinedVolumeBounds*> (&(comb->first()->volumeBounds()));
      const Trk::CombinedVolumeBounds* comb2 = dynamic_cast<const Trk::CombinedVolumeBounds*> (&(comb->second()->volumeBounds()));
      if (!comb1 && !comb2) {
        subVols.push_back(comb->second());
        comVol =comb->first(); comb = 0; 
      } else if (comb1) {
        comb = comb1;
        subVols.push_back(comb->second());
      } else {
        comb = comb2;
        subVols.push_back(comb->first());
      }
    }
    subVols.push_back(comVol);
    subVolumes.resize(subVols.size());
    double spbVol=0.; 
    for (unsigned int i=0;i<subVols.size();i++) { subVolumes[i] = calculateVolume(subVols[i]); spbVol += subVolumes[i]; } 
    // first rndm defines prism
    unsigned int iPrism = 0;
    double vCount = 0.;
    while ( iPrism < subVols.size() && vCount < rndm[0]*spbVol ) {
      vCount +=subVolumes[iPrism];
      iPrism++;
    }       
    vCount = subVolumes[iPrism-1] - (vCount-rndm[0]*spbVol);
    double zfr = vCount/subVolumes[iPrism-1];
    //std::cout << "prism picked:"<< iPrism-1<<","<< subVolumes.size()<<": z fraction:"<<zfr<< std::endl;
    // xy
    const Trk::PrismVolumeBounds* prism = dynamic_cast<const Trk::PrismVolumeBounds*> (&(subVols[iPrism-1]->volumeBounds()));
    if (prism) {
      std::vector<std::pair<double,double> > xy = prism->xyVertices();      
      //std::cout << "prism vertices:"<< xy.size()<< std::endl;
      gp = Trk::GlobalPosition(xy[2].first +sqrt(rndm[1])*( xy[0].first -xy[2].first 
							    + rndm[2]*( xy[1].first -xy[0].first  ) ),
			       xy[2].second+sqrt(rndm[1])*( xy[0].second-xy[2].second 
							    + rndm[2]*( xy[1].second-xy[0].second ) ),
			       -z + zfr*2*z);
    }
    
  } else {
    std::cout << "volume bounds not recognized in scan"<< std::endl;
  }
  if (!vol->inside(vol->transform()*gp,0.001)) std::cout << "test hit:wrong scan hit:"<<gp<< std::endl;
  
  return (vol->transform()*gp);
}


std::vector<std::pair<const Trk::Volume*,double> > 
Muon::MuonInertMaterialBuilder::splitComposedVolume(const Trk::Volume* trVol, bool estimateVol) const
{
  std::vector<const Trk::Volume*> garbage;
  std::vector<std::pair<const Trk::Volume*,const Trk::Volume*> > constituents;
  constituents.push_back(std::pair<const Trk::Volume*,const Trk::Volume*>(trVol,0));
  std::vector<std::pair<const Trk::Volume*,const Trk::Volume*> >::iterator sIter= constituents.begin(); 
  const Trk::Volume* subVol = 0;
  while (sIter!= constituents.end()) {
    const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*>
      (&((*sIter).first->volumeBounds()));
    const Trk::SubtractedVolumeBounds* sub = dynamic_cast<const Trk::SubtractedVolumeBounds*>
      (&((*sIter).first->volumeBounds()));
    if (comb) {
      subVol = (*sIter).second;
      sIter = constituents.erase(sIter);
      if (comb->intersection()) {
	Trk::Volume* newSubVol = new Trk::Volume(0,new Trk::SubtractedVolumeBounds(comb->first()->clone(),comb->second()->clone()));
        if (subVol) { 
	  Trk::Volume* newCSubVol = new Trk::Volume(0,new Trk::CombinedVolumeBounds(subVol->clone(),newSubVol,false));
	  constituents.insert(sIter,std::pair<const Trk::Volume*,const Trk::Volume*> (comb->first(),newCSubVol));        
	  garbage.push_back(newCSubVol);
	} else {
	  constituents.insert(sIter,std::pair<const Trk::Volume*,const Trk::Volume*> (comb->first(),newSubVol));
	  garbage.push_back(newSubVol);
	}
      } else {
	constituents.insert(sIter,std::pair<const Trk::Volume*,const Trk::Volume*> (comb->first(),subVol));
	constituents.insert(sIter,std::pair<const Trk::Volume*,const Trk::Volume*> (comb->second(),subVol));
      }
      sIter=constituents.begin();
    } else if (sub) {
      subVol = (*sIter).second;
      sIter = constituents.erase(sIter);
      if (subVol) {
	Trk::Volume* newSubVol = new Trk::Volume(0,new Trk::CombinedVolumeBounds(subVol->clone(),sub->inner()->clone(),false));
	constituents.insert(sIter,std::pair<const Trk::Volume*,const Trk::Volume*> (sub->outer(),newSubVol));
        garbage.push_back(newSubVol);
      } else {
	constituents.insert(sIter,std::pair<const Trk::Volume*,const Trk::Volume*> (sub->outer(),sub->inner()));
      }
      sIter = constituents.begin();
    } else {
      sIter++; 
    }    
  } 

  std::vector<std::pair<const Trk::Volume*,double> > wConst;
  for (unsigned int i=0;i<constituents.size();i++) 
    wConst.push_back(std::pair<const Trk::Volume*,double>(constituents[i].first,1.));     

  if (estimateVol && (constituents.size()>1 || constituents[0].second) ) {  
    for (unsigned int iv=0;iv<constituents.size();iv++) {
      const Trk::Volume* replaceVol = 0;
      if (constituents[iv].second) {   
	const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(constituents[iv].first->volumeBounds()));
	if (cyl && cyl->innerRadius()<0.001) {
	  // recalculate inner radius for better efficiency
	  std::vector<std::pair<const Trk::Volume*,double> > subtr=splitComposedVolume(constituents[iv].second,false);
	  for (unsigned int is=0;is<subtr.size();is++) {
	    const Trk::CylinderVolumeBounds* cyls = dynamic_cast<const Trk::CylinderVolumeBounds*> 
	      (&(subtr[is].first->volumeBounds()));
	    if (cyls) {
	      Trk::GlobalPosition dc = subtr[is].first->transform().inverse()*constituents[iv].first->center(); 
	      if (cyls->outerRadius()<cyl->outerRadius() && dc.perp()<0.001 && fabs(dc.z())<10.
		  && fabs(cyl->halflengthZ()-cyls->halflengthZ())<10.){
		replaceVol = new Trk::Volume(new HepTransform3D(constituents[iv].first->transform()),
					     new Trk::CylinderVolumeBounds(cyls->outerRadius(),
									   cyl->outerRadius(),
									   cyl->halfPhiSector(),
									   cyl->halflengthZ()   ) );
		break;
	      }
	    }
	  }
	}
      }
      double fr = 0.;
      const Trk::Volume* scanVol = replaceVol ? replaceVol : constituents[iv].first;
      double nHits = calculateVolume(scanVol)/1.e6;
      if (nHits<100.) nHits=100.;
      if (nHits>1000.) nHits=1000.;
      for (unsigned int ih=0;ih<nHits;ih++) {
	Trk::GlobalPosition gp=getScanPoint(scanVol);
	double w=1.;
	if (constituents[iv].second && constituents[iv].second->inside(gp,0.)) w=0.;
	else {
	  for (unsigned int icv=0;icv<constituents.size();icv++)
	    if (icv!=iv && constituents[icv].first->inside(gp,0.)
		&& (!constituents[icv].second || !constituents[icv].second->inside(gp,0.)))
	      w = 1./(1./w+1.); 
	}
	fr += w;
      }
      if (replaceVol) wConst[iv].first = replaceVol; 
      wConst[iv].second = fr/nHits; 
    }
  }

  for (size_t i=0; i<garbage.size();i++) delete garbage[i];

  return wConst;
}

