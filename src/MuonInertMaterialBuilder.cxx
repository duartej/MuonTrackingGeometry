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
  m_layerThicknessLimit(2.),
  m_debugMode(false),
  m_buildBT(true),
  m_buildECT(true),
  m_buildFeets(true),
  m_buildRails(1),
  m_buildShields(true),
  m_magFieldTool("Trk::MagneticFieldTool/AtlasMagneticFieldTool")
{
  declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
  declareProperty("MuonDetManagerLocation",           m_muonMgrLocation);
  declareProperty("SimplifyGeometry",                 m_simplify);
  declareProperty("SimplifyGeometryToLayers",         m_simplifyToLayers);
  declareProperty("LayerThicknessLimit",              m_layerThicknessLimit);
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

    log << MSG::INFO  << name() <<" initialize() successful" << endreq;

    getVolumeFractions();
    
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
      HepTransform3D* combTr = new HepTransform3D((*msTypeIter).second[it]); 
      const Trk::DetachedTrackingVolume* newStat = msTV->clone(msTypeName,*combTr);
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
	  
	  //	  std::cout << " INERT muon object found:" << vname <<" "<<ichild<< std::endl;
	  
	  bool accepted = true ;
	  if (  vname.substr(0,2)=="BT" || vname.substr(0,6) == "EdgeBT" || vname.substr(0,6) == "HeadBT" ) accepted = m_buildBT ? true : false;
	  else if ( vname.substr(0,3)=="ECT" ) accepted = m_buildECT ? true : false; 
	  else if ( vname.size()>7 && (vname.substr(3,4)=="Feet" || vname.substr(4,4)=="Feet" ) ) accepted = m_buildFeets ? true : false; 
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
                if (found && !simpleTree) {
		  //std::cout << "object repeated?:"<< protoName<<","<<vols[ish].second.size() << std::endl;
		  //std::cout << "already booked at positions:"<< std::endl;
                  //for (unsigned int j=0;j<objs[ip].second.size();j++) 
                  //  std::cout << j <<":"<<objs[ip].second[j].getTranslation()<< std::endl;
		  //std::cout << "additionally requested at positions:"<< std::endl;
                  //for (unsigned int j=0;j<vols[ish].second.size();j++) 
                  //  std::cout << j <<":"<<vols[ish].second[j].getTranslation()<< std::endl;
		}
                if (simpleTree) objs[ip].second.push_back(vol.getTransform());
                //else objs[ip].second.insert(objs[ip].second.end(),vols[ish].second.begin(),vols[ish].second.end());
	      } 
	    }  
	    if (found) continue;
	    //std::cout << "decoding shape:"<< protoName << std::endl;
            //decodeShape(input_shapes[ish]);
	    const Trk::Volume* trObject = m_geoShapeConverter->translateGeoShape(input_shapes[ish],new HepTransform3D());
	    if (trObject) {  
	      Trk::MaterialProperties mat = m_materialConverter->convert( vols[ish].first->getMaterial() );
	      const Trk::TrackingVolume* newType= new Trk::TrackingVolume( *trObject, mat, m_muonMagneticField,0,0,protoName);
	      const Trk::TrackingVolume* simType = simplifyShape(newType);
	      const Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(protoName,simType);
	      objs.push_back(std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> >(typeStat,vols[ish].second));
	    }  else {
	      log << MSG::WARNING << name()<< " volume not translated: " << vname << std::endl;
	    }            
	  } // end new object
	}
	vol.next();      	
      }
    }

    const std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >* mObjects = new std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >(objs);

    log << MSG::INFO << name() << " returns " << mObjects->size() << " prototypes " << endreq;   
    
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
  const Trk::Volume* envelope = 0;
  int simpleMode = 0;
  double maxD    = 45.;
  // retrieve fraction of material in the envelope
  double fraction = 0.;
  for (unsigned int ivol=0; ivol<m_volFractions.size(); ivol++) {
    if (m_volFractions[ivol].first == trVol->volumeName() ) { fraction = m_volFractions[ivol].second.second ; break; }
  }         
  // 
  std::vector<const Trk::Volume*>* constituents = new std::vector<const Trk::Volume*>;
  std::vector<const Trk::Volume*>* subtractions = new std::vector<const Trk::Volume*>;
  constituents->push_back(trVol);
  std::vector<const Trk::Volume*>::iterator sIter= constituents->begin(); 
  while (sIter!= constituents->end()) {
    const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*> (&((*sIter)->volumeBounds()));
    const Trk::SubtractedVolumeBounds* sub = dynamic_cast<const Trk::SubtractedVolumeBounds*> (&((*sIter)->volumeBounds()));
    if (comb) {
      sIter = constituents->erase(sIter);
      constituents->insert(sIter,comb->first());
      constituents->insert(sIter,comb->second());
      sIter=constituents->begin();
    } else if (sub) {
      sIter = constituents->erase(sIter);
      constituents->insert(sIter,sub->outer());
      subtractions->push_back(sub->inner());
      sIter = constituents->begin();
    } else {
      sIter++; 
    }    
  } 
  
  if (constituents->size() == 1) {  // easy case

    envelope = new Trk::Volume(*(constituents->front()));

    if ( m_simplify || (trVol->volumeName().substr(0,4)=="Rail" && m_buildRails<2) ) { 
      if ( subtractions->size()==0 || trVol->volumeName().substr(0,4)=="ECTC" ) {
        simpleMode = 2;                         // keep original volume
      } else {
	double thin = thinDim(envelope);
	if ( thin < maxD ) simpleMode = 1;      // layers
	else simpleMode = 3;                    // envelope with diluted material 
      }
    }

  } else {
    
    double xSize = 0.;
    double ySize = 0.;
    double zSize = 0.; 
    double rSize = -1.;

    sIter = constituents->begin();
    while (sIter!= constituents->end()) {
      const Trk::CylinderVolumeBounds*  cylBounds = dynamic_cast<const Trk::CylinderVolumeBounds*>  (&((*sIter)->volumeBounds()));
      const Trk::CuboidVolumeBounds*    cubBounds = dynamic_cast<const Trk::CuboidVolumeBounds*>    (&((*sIter)->volumeBounds()));
      const Trk::TrapezoidVolumeBounds* trdBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&((*sIter)->volumeBounds()));
      if (cylBounds) {  
        double rOut = cylBounds->outerRadius();
        double hZ   = cylBounds->halflengthZ();
        rSize = fmax(rSize,rOut);
	Trk::GlobalPosition gpp = (*sIter)->transform() * Trk::GlobalPosition( 0., 0., hZ );
	Trk::GlobalPosition gpm = (*sIter)->transform() * Trk::GlobalPosition( 0., 0.,-hZ );
	Trk::GlobalDirection gpd = (*sIter)->transform().getRotation() * Trk::GlobalDirection( 0., 0.,1. );
        xSize = fmax( xSize, fabs(gpp.x())+rOut*(1.-fabs(gpd[0])) );
        ySize = fmax( ySize, fabs(gpp.y())+rOut*(1.-fabs(gpd[1])) );
        zSize = fmax( zSize, fabs(gpp.z())+rOut*(1.-fabs(gpd[2])) );
        xSize = fmax( xSize, fabs(gpm.x())+rOut*(1.-fabs(gpd[0])) );
        ySize = fmax( ySize, fabs(gpm.y())+rOut*(1.-fabs(gpd[1])) );
        zSize = fmax( zSize, fabs(gpm.z())+rOut*(1.-fabs(gpd[2])) );
      }
      if (cubBounds) {  
        double x = cubBounds->halflengthX();
        double y = cubBounds->halflengthY();
        double z = cubBounds->halflengthZ();
	std::vector<Trk::GlobalPosition> edges;
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x, y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x, y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x,-y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x,-y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x, y,-z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x, y,-z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x,-y,-z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x,-y,-z));
	for (unsigned int ie=0; ie < edges.size(); ie++) {
	  xSize = fmax( xSize, fabs(edges[ie].x()) );
	  ySize = fmax( ySize, fabs(edges[ie].y()) );
	  zSize = fmax( zSize, fabs(edges[ie].z()) );
        }
      }
      if (trdBounds) {  
        double x1 = trdBounds->minHalflengthX();
        double x2 = trdBounds->maxHalflengthX();
        double y  = trdBounds->halflengthY();
        double z  = trdBounds->halflengthZ();
	std::vector<Trk::GlobalPosition> edges;
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x2, y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x2, y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x1,-y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x1,-y, z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x2, y,-z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x2, y,-z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition( x1,-y,-z));
	edges.push_back((*sIter)->transform()*Trk::GlobalPosition(-x1,-y,-z));
	for (unsigned int ie=0; ie < edges.size(); ie++) {
	  xSize = fmax( xSize, fabs(edges[ie].x()) );
	  ySize = fmax( ySize, fabs(edges[ie].y()) );
	  zSize = fmax( zSize, fabs(edges[ie].z()) );
        }
      }
      sIter++;
    }

    if (rSize>0. && xSize <=rSize && ySize <=rSize ) {
      if ( m_simplify && fraction > 0.1 ) {
        simpleMode = 3;            // envelope with diluted material
      }
      envelope = new Trk::Volume(new HepTransform3D(trVol->transform()),new Trk::CylinderVolumeBounds(rSize,zSize));
    } else { 
      if ( (m_simplify || (trVol->volumeName().substr(0,4)=="Rail" && m_buildRails<2))  && fraction > 0.1 ) {
        simpleMode = 3;           // envelope with diluted material 
      }
      if (m_simplify && trVol->volumeName()=="BTVoussoir" ) { // merge with EdgeBTVoussoir
	std::pair<double,double> fedge;
	for (unsigned int ivol=0; ivol<m_volFractions.size(); ivol++) {
	  if (m_volFractions[ivol].first == "EdgeBTVoussoir" ) { fedge = m_volFractions[ivol].second ; break; }
	}
	envelope = new Trk::Volume(new HepTransform3D(trVol->transform()),new Trk::TrapezoidVolumeBounds(2026.,2270.,ySize,zSize));
	fraction = ( 8*fraction*xSize*ySize*zSize + fedge.first*fedge.second )/8/2148/ySize/zSize;
      } else {
	envelope = new Trk::Volume(new HepTransform3D(trVol->transform()),new Trk::CuboidVolumeBounds(xSize,ySize,zSize));
      }
    }
  }

  if (!envelope) envelope = new Trk::Volume(*trVol);

  const Trk::TrackingVolume* newVol;
  std::vector<const Trk::TrackingVolume*>* confinedVols = new std::vector<const Trk::TrackingVolume*>;

  //if ( m_simplifyToLayers && !simpleMode ) simpleMode = 1;
  if ( m_simplifyToLayers ) simpleMode = 1;

  if ( m_simplify  && !simpleMode ) {
    if (trVol->volumeName().substr(0,6)=="ECTWal" || trVol->volumeName().substr(0,6)=="ECTKey" ||
        trVol->volumeName().substr(0,6)=="BTRibE" ) simpleMode=5;
    else simpleMode = 4;
  }

  if ( simpleMode==1 ) {
    confinedVols->push_back(trVol);
    std::pair<std::vector<const Trk::Layer*>*,std::vector<const Trk::TrackingVolume*>* >  confinedObjs = translateToLayers(confinedVols, simpleMode-1 );
    std::string envName=trVol->volumeName()+"_envelope";
    newVol = new Trk::TrackingVolume( *envelope, confinedObjs.first, confinedObjs.second, m_muonMaterial,m_muonMagneticField,envName);
    if (confinedObjs.second) {
      for (unsigned int iv = 0; iv < confinedObjs.second->size(); iv++)
	Trk::TrackingVolumeManipulator::confineVolume(*((*(confinedObjs.second))[iv]),newVol);
    }
  } else if ( simpleMode==2 ) {
    newVol = trVol;
    delete confinedVols;    
  } else if ( simpleMode==3 ) {
    std::string envName=trVol->volumeName();
    Trk::MaterialProperties mat(1.,trVol->x0()/fraction,fraction*trVol->zOverAtimesRho());
    newVol = new Trk::TrackingVolume( *envelope, mat, m_muonMagneticField, 0, 0, envName);  
    delete trVol;  
    delete confinedVols;
  } else if ( simpleMode==4 ) {
    std::vector<const Trk::Layer*>* confinedLays = translateBoundariesToLayers(envelope,trVol,fraction);
    std::string envName=trVol->volumeName()+"_envelope";
    newVol = new Trk::TrackingVolume( *envelope, m_muonMaterial, m_muonMagneticField, confinedLays, envName);    
    delete trVol;
    delete confinedVols;
  } else if ( simpleMode==5 ) {
    confinedVols->push_back(trVol);
    std::pair<std::vector<const Trk::Layer*>*,std::vector<const Trk::TrackingVolume*>* >  confinedObjs = translateToLayers(confinedVols, simpleMode-1 );
    std::string envName=trVol->volumeName()+"_envelope";
    newVol = new Trk::TrackingVolume( *envelope, confinedObjs.first, confinedObjs.second, m_muonMaterial,m_muonMagneticField,envName);
    if (confinedObjs.second) {
      for (unsigned int iv = 0; iv < confinedObjs.second->size(); iv++)
	Trk::TrackingVolumeManipulator::confineVolume(*((*(confinedObjs.second))[iv]),newVol);
    }
  } else {    
    confinedVols->push_back(trVol);
    std::string envName=trVol->volumeName()+"_envelope";
    newVol = new Trk::TrackingVolume( *envelope, m_muonMaterial, m_muonMagneticField, confinedVols, envName);    
  }

  delete envelope;
  return newVol;
}

std::vector<const Trk::Layer*>*  Muon::MuonInertMaterialBuilder::translateBoundariesToLayers(const Trk::Volume* envelope, const Trk::TrackingVolume* trVol, double fraction) const
{
  std::vector<const Trk::Layer*>* lays = 0;
  if (!envelope) return lays; 
  double thickness = 0.;
  
  const Trk::CylinderVolumeBounds*   cyl  = dynamic_cast<const Trk::CylinderVolumeBounds*>   (&(envelope->volumeBounds()));

  const std::vector<const Trk::Surface*>* bSurfs = envelope->volumeBounds().decomposeToSurfaces(envelope->transform());

  if (cyl) {
    for (unsigned int ib=0; ib< bSurfs->size(); ib++ ){
      const Trk::Surface* surf = (*bSurfs)[ib];

      if (ib<2) thickness = (envelope->transform().inverse()*surf->center()).mag()*fraction/2;
      else thickness = (cyl->outerRadius()-cyl->innerRadius())*fraction/2;
  
      const Trk::Layer* lay = boundarySurfaceToLayer(*surf,trVol, thickness);
      if (lay) {
	lay->setLayerType(0); 
	if (!lays) lays = new std::vector<const Trk::Layer*>;
	lays->push_back( lay );
      }
    }

  } else {
    for (unsigned int ib=0; ib< bSurfs->size(); ib++ ){
      const Trk::Surface* surf = (*bSurfs)[ib];
      thickness = (envelope->transform().inverse()*surf->center()).mag()*fraction/3; 
      const Trk::Layer* lay = boundarySurfaceToLayer(*surf,trVol, thickness); 
      if (lay) {
	lay->setLayerType(0); 
	if (!lays) lays = new std::vector<const Trk::Layer*>;
	lays->push_back( lay );
      }
    }
  }
  for (unsigned int ib=0; ib< bSurfs->size(); ib++ ) delete (*bSurfs)[ib];
  delete bSurfs;
  
  return lays; 
}

std::pair<std::vector<const Trk::Layer*>*,std::vector<const Trk::TrackingVolume*>* >  Muon::MuonInertMaterialBuilder::translateToLayers(const std::vector<const Trk::TrackingVolume*>* vols, int mode) const
{
  std::vector<const Trk::Layer*>* lays = 0;
  std::vector<const Trk::TrackingVolume*>* dVols = 0;
  if (!vols) return std::pair<std::vector<const Trk::Layer*>*,std::vector<const Trk::TrackingVolume*>* > (0,0); 
  double tol = 0.001;
    
  for (unsigned int i=0; i< vols->size(); i++) {
    
    if ( mode!=0 ) {
      double thickness = 2000.;
      const std::vector< Trk::SharedObject<const Trk::BoundarySurface<Trk::TrackingVolume> > > bounds = (*vols)[i]->boundarySurfaces();
      for (unsigned int ib=0; ib< bounds.size(); ib++ ){
	thickness = 2000.;
	const Trk::Surface& surf = (bounds[ib].getPtr())->surfaceRepresentation();
	for (unsigned int ib2=ib+1; ib2< bounds.size(); ib2++ ){
	  const Trk::Surface& surf2 = (bounds[ib2].getPtr())->surfaceRepresentation();
          if ( fabs(surf.normal().dot(surf2.normal()))>0.99 ) {
            double dist = (surf.center()-surf2.center()).dot(surf.normal());
            if ( !bounds[ib].getPtr()->attachedVolume(surf.center(),surf.normal(),Trk::alongMomentum) && dist > tol && dist < thickness ) thickness = dist; 
            if ( bounds[ib].getPtr()->attachedVolume(surf.center(),surf.normal(),Trk::alongMomentum) && dist < tol && dist > -thickness ) thickness = -dist; 
	  }
	}
        if ( thickness > tol && thickness < 100 ) {
	  const Trk::Layer* lay = boundarySurfaceToLayer(surf,(*vols)[i], thickness); 
	  if (lay) {
            lay->setLayerType(0);
	    if (!lays) lays = new std::vector<const Trk::Layer*>;
	    lays->push_back( lay );
	  }
	} else if ((*vols)[i]->volumeName().substr(0,3)=="ECT" ) {
          const Trk::CylinderBounds* cyl = dynamic_cast<const Trk::CylinderBounds*> (&(surf.bounds()));
          if (cyl) {
	    const Trk::Layer* lay = boundarySurfaceToLayer(surf,(*vols)[i], 70.); 
	    if (lay) {
	      lay->setLayerType(0);
	      if (!lays) lays = new std::vector<const Trk::Layer*>;
	      lays->push_back( lay );
	    }
	  }
	}
      }
      //delete (*vols)[i];
    } else {
      std::vector<const Trk::Layer*>* temp_layers= new std::vector<const Trk::Layer*>;
      double thick = volumeToLayers(*temp_layers,(*vols)[i],0,(*vols)[i]);
      if ( m_layerThicknessLimit<0. || thick<m_layerThicknessLimit ) {
	if (!lays && temp_layers->size()) lays = new std::vector<const Trk::Layer*>;
	for (unsigned int il=0;il<temp_layers->size();il++) { (*temp_layers)[il]->setLayerType(0); lays->push_back((*temp_layers)[il]);}
	//delete (*vols)[i];
      } else {
        if (!dVols) dVols = new std::vector<const Trk::TrackingVolume*>;
	dVols->push_back((*vols)[i]);
	for (unsigned int il=0;il<temp_layers->size();il++) delete (*temp_layers)[il];
      }
      delete temp_layers;
    }
  }
  
  return std::pair<std::vector<const Trk::Layer*>*,std::vector<const Trk::TrackingVolume*>* > (lays,dVols);
}

double Muon::MuonInertMaterialBuilder::volumeToLayers(std::vector<const Trk::Layer*>& lays, const Trk::Volume* vol, Trk::Volume* subtrVol, const Trk::MaterialProperties* mat) const
{
  double thickness = 0.;  
  if (!vol) return thickness;
  double thInX0 = 0.;
  
  const Trk::CombinedVolumeBounds*   comb = dynamic_cast<const Trk::CombinedVolumeBounds*>   (&(vol->volumeBounds()));
  const Trk::SubtractedVolumeBounds* sub  = dynamic_cast<const Trk::SubtractedVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::TrapezoidVolumeBounds*  trap = dynamic_cast<const Trk::TrapezoidVolumeBounds*>  (&(vol->volumeBounds()));
  const Trk::CuboidVolumeBounds*     box  = dynamic_cast<const Trk::CuboidVolumeBounds*>     (&(vol->volumeBounds()));
  const Trk::CylinderVolumeBounds*   cyl  = dynamic_cast<const Trk::CylinderVolumeBounds*>   (&(vol->volumeBounds()));
  
  if (!comb && !sub && !trap && !box && !cyl ) std::cout << "ERROR: unknown volume boundaries!!!" << std::endl;
  
  if (comb) {
    thInX0 = fmax(volumeToLayers(lays,comb->first(),subtrVol,mat),volumeToLayers(lays,comb->second(),subtrVol,mat));
  } 
  if (sub) {
    if (subtrVol) {
      // here is necessary to combine subtracted volumes
      Trk::Volume* dsub = new Trk::Volume(0,new Trk::CombinedVolumeBounds(sub->inner(),subtrVol,false));
      thInX0 = volumeToLayers(lays,sub->outer(),dsub,mat);
    } else {
      thInX0 = volumeToLayers(lays,sub->outer(),sub->inner(),mat);
    }
  } 
  if (box) {
    //if (!checkVolume(vol)) std::cout << "problems in volume boundaries" << std::endl;
    double hx = box->halflengthX();
    double hy = box->halflengthY();
    double hz = box->halflengthZ();
    const std::vector<const Trk::Surface*>* box_surf = box->decomposeToSurfaces(vol->transform());    
      
    Trk::PlaneSurface* plane=0;
    if ( hz<=hx && hz<=hy ) { // x-y plane
      const Trk::PlaneSurface* pl = dynamic_cast<const Trk::PlaneSurface*> ((*box_surf)[0]);
      if (pl) plane = new Trk::PlaneSurface(*pl,HepTransform3D(pl->transform().inverse()*vol->transform()) );  
      thickness = 2*hz;
    } else if ( hx<=hy && hx<=hz ) {
      const Trk::PlaneSurface* pl = dynamic_cast<const Trk::PlaneSurface*> ((*box_surf)[2]);
      if (pl) plane = new Trk::PlaneSurface(*pl,HepTransform3D(pl->transform().inverse()*vol->transform()
								   *HepRotateY3D(90.*deg)*HepRotateZ3D(90.*deg)) );
      thickness = 2*hx;
    } else if ( hy<=hx && hy<=hz ) {
      const Trk::PlaneSurface* pl = dynamic_cast<const Trk::PlaneSurface*> ((*box_surf)[4]);
      if (pl) plane = new Trk::PlaneSurface(*pl,HepTransform3D(pl->transform().inverse()*vol->transform()
								   *HepRotateY3D(-90.*deg)*HepRotateX3D(-90.*deg)) );
      thickness = 2*hy;
    }
    //material
    Trk::MaterialProperties material(thickness,mat->x0(),mat->zOverAtimesRho(),mat->averageZ(),mat->dEdX());  
    Trk::HomogenousLayerMaterial some_mat(material, Trk::oppositePre);
    thInX0 = material.thicknessInX0();  
    
    if (subtrVol) {
      Trk::Volume* subVol = createSubtractedVolume(plane->transform(), subtrVol);
      Trk::SubtractedPlaneSurface* subtrSurf = new Trk::SubtractedPlaneSurface(*plane,
									       new Trk::VolumeExcluder(subVol), false );
      lays.push_back(new Trk::SubtractedPlaneLayer(subtrSurf, some_mat,thickness));
      delete subtrSurf;
    } else {
      lays.push_back(new Trk::PlaneLayer(plane,some_mat, thickness));
    }
    
    for (unsigned int ib=0; ib< box_surf->size(); ib++ ) delete (*box_surf)[ib];
    delete box_surf;
  }
  
  if(trap){
    //if (!checkVolume(vol)) std::cout << "problems in volume boundaries" << std::endl;
    double hxmin = trap->minHalflengthX();
    double hxmax = trap->maxHalflengthX();
    double hy = trap->halflengthY();
    double hz = trap->halflengthZ();
    const std::vector<const Trk::Surface*>* trap_surf = trap->decomposeToSurfaces(vol->transform());    
      
    Trk::PlaneSurface* plane=0;
    if ( hz<=0.5*(hxmin+hxmax) && hz<=hy ) { // x-y plane
      const Trk::PlaneSurface* pl = dynamic_cast<const Trk::PlaneSurface*> ((*trap_surf)[0]);
      if (pl) plane = new Trk::PlaneSurface(*pl,HepTranslate3D(vol->center()- pl->center()) );  
      thickness = 2*hz;
    } else if ( 0.5*(hxmin+hxmax)<=hy && 0.5*(hxmin+hxmax)<=hz ) {
      const Trk::PlaneSurface* pl = dynamic_cast<const Trk::PlaneSurface*> ((*trap_surf)[2]);
      if (pl) plane = new Trk::PlaneSurface(*pl,HepTransform3D(pl->transform().inverse()*vol->transform()
								   *HepRotateY3D(90.*deg)*HepRotateZ3D(90.*deg)) );
      thickness = hxmin+hxmax;
    } else if ( hy<=0.5*(hxmin+hxmax) && hy<=hz ) {
      const Trk::PlaneSurface* pl = dynamic_cast<const Trk::PlaneSurface*> ((*trap_surf)[4]);
      if (pl) plane = new Trk::PlaneSurface(*pl,HepTransform3D(pl->transform().inverse()*vol->transform()
								   *HepRotateY3D(-90.*deg)*HepRotateX3D(-90.*deg)) );
      thickness = 2*hy;
    }
    //material
    Trk::MaterialProperties material(thickness,mat->x0(),mat->zOverAtimesRho(),mat->averageZ(),mat->dEdX());  
    Trk::HomogenousLayerMaterial some_mat(material, Trk::oppositePre);
    thInX0 = material.thicknessInX0();  
    if (subtrVol) {
      Trk::Volume* subVol = createSubtractedVolume(plane->transform(), subtrVol);
      Trk::SubtractedPlaneSurface* subtrSurf = new Trk::SubtractedPlaneSurface(*plane,
									       new Trk::VolumeExcluder(subVol), false );
      lays.push_back(new Trk::SubtractedPlaneLayer(subtrSurf, some_mat,thickness));
      delete subtrSurf;
    } else {
      lays.push_back(new Trk::PlaneLayer(plane, some_mat, thickness));
    }
    for (unsigned int ib=0; ib< trap_surf->size(); ib++ ) delete (*trap_surf)[ib];
    delete trap_surf;
  }
  
  if(cyl){
    const std::vector<const Trk::Surface*>* cyl_surf = cyl->decomposeToSurfaces(vol->transform());    
    //double radius = cyl->mediumRadius();
    double drad   = cyl->deltaRadius();
    double hz     = cyl->halflengthZ();

    if ( hz < drad ) {       // disc/ellipse surface
      const Trk::PlaneSurface*   plane = dynamic_cast<const Trk::PlaneSurface*> ((*cyl_surf)[0]);
      const Trk::DiscSurface*     disc = dynamic_cast<const Trk::DiscSurface*> ((*cyl_surf)[0]);
      //material
      Trk::MaterialProperties material(hz,mat->x0(),mat->zOverAtimesRho(),mat->averageZ(),mat->dEdX());  
      Trk::HomogenousLayerMaterial some_mat(material, Trk::oppositePre);
      thInX0 = material.thicknessInX0();  
      
      if ( plane ) {
        Trk::PlaneSurface* pl = new Trk::PlaneSurface(*plane,HepTranslate3D(vol->center()- plane->center()) );  
	if (subtrVol) {
	  Trk::Volume* subVol = createSubtractedVolume(pl->transform(), subtrVol);
	  Trk::SubtractedPlaneSurface* subtrSurf = new Trk::SubtractedPlaneSurface(*pl, new Trk::VolumeExcluder(subVol), false );
	  lays.push_back(new Trk::SubtractedPlaneLayer(subtrSurf,some_mat,hz));
          delete subtrSurf;
	} else {
	  lays.push_back(new Trk::PlaneLayer(pl,some_mat,hz));
	} 
      }	else if (disc) {
	const Trk::DiscBounds* db = dynamic_cast<const Trk::DiscBounds*> (&(disc->bounds())); 
	Trk::PlaneSurface* di = new Trk::PlaneSurface(new HepTransform3D(HepTranslate3D(vol->center())
									 *HepRotate3D(disc->transform().getRotation())),  
						      new Trk::EllipseBounds(db->rMin(),db->rMin(),db->rMax(), db->rMax(),
									     db->halfPhiSector()));
 	if (subtrVol) {
	  Trk::Volume* subVol = createSubtractedVolume(di->transform(), subtrVol);
	  Trk::SubtractedPlaneSurface* subtrSurf = new Trk::SubtractedPlaneSurface(*di, new Trk::VolumeExcluder(subVol), false );
	  lays.push_back(new Trk::SubtractedPlaneLayer(subtrSurf,some_mat,hz));
          delete subtrSurf;
	} else {
	  lays.push_back(new Trk::PlaneLayer(di,some_mat,hz));
	}
      }
    } else if ( hz > 200. && drad > 100. ) {   // set of disc/ellipse surface      
      const Trk::PlaneSurface*   plane = dynamic_cast<const Trk::PlaneSurface*> ((*cyl_surf)[0]);
      const Trk::DiscSurface*     disc = dynamic_cast<const Trk::DiscSurface*> ((*cyl_surf)[0]);
      unsigned int numSlice = int(hz/200.)+1;
      
      for (unsigned int islice = 0; islice < numSlice; islice++) {
	// distance from center
	double d = hz * ( 2*islice/numSlice -1 ); 
	//material
	Trk::MaterialProperties material(2*hz/numSlice,mat->x0(),mat->zOverAtimesRho(),mat->averageZ(),mat->dEdX());  
	Trk::HomogenousLayerMaterial some_mat(material, Trk::oppositePre);
	thInX0 = material.thicknessInX0();  
	if ( plane ) {
	  Trk::PlaneSurface* pl = new Trk::PlaneSurface(*plane,HepTranslate3D((1.+d/hz)*(vol->center()- plane->center())) );  
	  if (subtrVol) {
	    Trk::Volume* subVol = createSubtractedVolume(pl->transform(), subtrVol);
	    Trk::SubtractedPlaneSurface* subtrSurf = new Trk::SubtractedPlaneSurface(*pl, new Trk::VolumeExcluder(subVol), false );
	    lays.push_back(new Trk::SubtractedPlaneLayer(subtrSurf,some_mat,2*hz/numSlice));
            delete subtrSurf;
	  } else {
	    lays.push_back(new Trk::PlaneLayer(pl,some_mat,2*hz/numSlice));
	  } 
	} else if (disc) {
	  const Trk::DiscBounds* db = dynamic_cast<const Trk::DiscBounds*> (&(disc->bounds())); 
	  Trk::PlaneSurface* di = new Trk::PlaneSurface(new HepTransform3D(HepRotate3D(disc->transform().getRotation())
									   *HepTranslate3D(vol->center())*HepTranslateZ3D(d)),  
							new Trk::EllipseBounds(db->rMin(),db->rMin(),db->rMax(), db->rMax(),
									       db->halfPhiSector()));
	  if (subtrVol) {
	    Trk::Volume* subVol = createSubtractedVolume(di->transform(), subtrVol);
	    Trk::SubtractedPlaneSurface* subtrSurf = new Trk::SubtractedPlaneSurface(*di, new Trk::VolumeExcluder(subVol), false );
	    lays.push_back(new Trk::SubtractedPlaneLayer(subtrSurf,some_mat,2*hz/numSlice));
            delete subtrSurf; 
	  } else {
	    lays.push_back(new Trk::PlaneLayer(di,some_mat,2*hz/numSlice));
	  } 
	}
      }
    } else {           // cylinder layer	
      //material
      Trk::MaterialProperties material(drad,mat->x0(),mat->zOverAtimesRho(),mat->averageZ(),mat->dEdX());  
      Trk::HomogenousLayerMaterial some_mat(material, Trk::oppositePre);
      thInX0 = material.thicknessInX0();  

      const Trk::CylinderSurface*   cylSurf = dynamic_cast<const Trk::CylinderSurface*> ((*cyl_surf)[2]);
      Trk::CylinderSurface* pSurf = new Trk::CylinderSurface(*cylSurf);  
      
      if (subtrVol) {
	Trk::Volume* subVol = createSubtractedVolume(pSurf->transform(), subtrVol);
	Trk::SubtractedCylinderSurface* subtrSurf = new Trk::SubtractedCylinderSurface(*pSurf, new Trk::VolumeExcluder(subVol),
										       false );
	lays.push_back(new Trk::SubtractedCylinderLayer(subtrSurf,some_mat,drad));
        delete subtrSurf;
      } else {
	lays.push_back(new Trk::CylinderLayer(pSurf,some_mat,drad));
      }
    }

    for (unsigned int ib=0; ib< cyl_surf->size(); ib++ ) delete (*cyl_surf)[ib];
    delete cyl_surf;
  }
  return thInX0;
}

const Trk::Layer* Muon::MuonInertMaterialBuilder::boundarySurfaceToLayer( const Trk::Surface& surf, const Trk::MaterialProperties* mat, double thickness) const
{
  // the method creates copy of the input surface

  const Trk::Layer* layer = 0;

  Trk::MaterialProperties material(thickness,mat->x0(),mat->zOverAtimesRho(),mat->averageZ(),mat->dEdX());  
  Trk::HomogenousLayerMaterial layMat(material, Trk::oppositePre);

  const Trk::SubtractedPlaneSurface* subPlane = dynamic_cast<const Trk::SubtractedPlaneSurface*> (&surf);
  const Trk::SubtractedCylinderSurface* subCyl = dynamic_cast<const Trk::SubtractedCylinderSurface*> (&surf);
  const Trk::PlaneSurface* plane = dynamic_cast<const Trk::PlaneSurface*> (&surf);
  const Trk::CylinderSurface* cyl = dynamic_cast<const Trk::CylinderSurface*> (&surf);

  if (subCyl)  {
    layer = new Trk::SubtractedCylinderLayer(new Trk::SubtractedCylinderSurface(*subCyl),layMat,thickness); 
    return layer;
  } else if (subPlane) {
    layer = new Trk::SubtractedPlaneLayer(new Trk::SubtractedPlaneSurface(*subPlane),layMat,thickness); 
    return layer;
  } else if (cyl) {
    layer = new Trk::CylinderLayer(new Trk::CylinderSurface(*cyl),layMat,thickness); 
    return layer;
  } else if (plane) {
    layer = new Trk::PlaneLayer(new Trk::PlaneSurface(*plane),layMat,thickness); 
    return layer;
  }
  return layer;
}
 
Trk::Volume* Muon::MuonInertMaterialBuilder::createSubtractedVolume(const HepTransform3D& transf, Trk::Volume* subtrVol) const
{
  Trk::Volume* subVol = 0;
  if (!subtrVol) return subVol;

  subVol = new Trk::Volume( *subtrVol, transf );
  
  return subVol;
}

const bool  Muon::MuonInertMaterialBuilder::checkVolume(const Trk::Volume* chVol) const
{
  const std::vector<const Trk::Surface*>* surf = chVol->volumeBounds().decomposeToSurfaces(chVol->transform());
  std::vector<const Trk::GlobalPosition*> bd_apexes;
  bool cylinder = false;
  
  for( unsigned int i = 0; i< surf->size(); i++ ){
    const Trk::PlaneSurface* plane = dynamic_cast<const Trk::PlaneSurface*> ((*surf)[i]);
    const Trk::CylinderSurface* cyl = dynamic_cast<const Trk::CylinderSurface*> ((*surf)[i]);
    if (cyl) { cylinder = true; break; }
    if (plane) {
      const Trk::TrapezoidBounds* trd  = dynamic_cast<const Trk::TrapezoidBounds*> (&(plane->bounds()));
      const Trk::RectangleBounds* rect = dynamic_cast<const Trk::RectangleBounds*> (&(plane->bounds()));
      if(trd){
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition( trd->maxHalflengthX(), trd->halflengthY())));  
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition(-trd->maxHalflengthX(), trd->halflengthY())));  
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition( trd->minHalflengthX(),-trd->halflengthY())));  
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition(-trd->minHalflengthX(),-trd->halflengthY()))); 
      }
      if(rect){
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition( rect->halflengthX(), rect->halflengthY())));  
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition(-rect->halflengthX(), rect->halflengthY())));  
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition( rect->halflengthX(),-rect->halflengthY())));  
	bd_apexes.push_back(plane->localToGlobal(Trk::LocalPosition(-rect->halflengthX(),-rect->halflengthY()))); 
      }
    }
  }

  unsigned int n_apexes = 4; //nr of apexes for one surface
  double tol = 0.001;
  unsigned int n_surf = surf->size(); 
  int apex_s[16];
  for( unsigned int i = 0; i < n_surf; i++) apex_s[i] = 0;
 
  for( unsigned int i = 0; i < n_apexes; i++){
    for( unsigned int j =  2*n_apexes; j < bd_apexes.size(); j++){
      if( (fabs( bd_apexes[i]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i]->z() - bd_apexes[j]->z()) < tol)) apex_s[0]++;
      if( (fabs( bd_apexes[i+n_apexes]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->z() - bd_apexes[j]->z()) < tol)) apex_s[1]++;
    }
  }

  for( unsigned int i = 2*n_apexes; i < 3*n_apexes; i++){
    for( unsigned int j =  4*n_apexes; j < bd_apexes.size(); j++){
      if( (fabs( bd_apexes[i]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i]->z() - bd_apexes[j]->z()) < tol)) apex_s[2]++;
      if( (fabs( bd_apexes[i+n_apexes]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->z() - bd_apexes[j]->z()) < tol)) apex_s[3]++;
    }
  }

  for( unsigned int i = 2*n_apexes; i < 3*n_apexes; i++){
    for( unsigned int j =  0; j < 2*n_apexes; j++){
      if( (fabs( bd_apexes[i]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i]->z() - bd_apexes[j]->z()) < tol)) apex_s[2]++;
      if( (fabs( bd_apexes[i+n_apexes]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->z() - bd_apexes[j]->z()) < tol)) apex_s[3]++;
    }
  }

  for( unsigned int i = 4*n_apexes; i < 5*n_apexes; i++){
    for( unsigned int j =  0; j < 4*n_apexes;  j++){
      if( (fabs( bd_apexes[i]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i]->z() - bd_apexes[j]->z()) < tol)) apex_s[4]++;
      if( (fabs( bd_apexes[i+n_apexes]->x() - bd_apexes[j]->x()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->y() - bd_apexes[j]->y()) < tol) &&
	  (fabs( bd_apexes[i+n_apexes]->z() - bd_apexes[j]->z()) < tol)) apex_s[5]++;
    }
  }

  
  for( unsigned int i = 0; i < n_surf; i++){
    if( apex_s[i] != 8 ){
      for( unsigned int j = 0; j < n_surf; j++) {
	const Trk::PlaneSurface* plane = dynamic_cast<const Trk::PlaneSurface*> ((*surf)[j]);
	if (plane) {
	  const Trk::TrapezoidBounds* trd  = dynamic_cast<const Trk::TrapezoidBounds*> (&(plane->bounds()));
	  const Trk::RectangleBounds* rect = dynamic_cast<const Trk::RectangleBounds*> (&(plane->bounds()));
	  std::cout<<"\tSurf["<<j<<"] common appexes::(trd,rect)" << trd <<"," << rect <<","<<apex_s[j]<<std::endl;
        }
      } 
      return false;
    }
  }	  
  return true;
}

void  Muon::MuonInertMaterialBuilder::splitShape(const GeoShape* sh, std::vector<const GeoShape*>& shapes) const
{
  if ( sh->type()=="Union" ) {
    const GeoShapeUnion* sub = dynamic_cast<const GeoShapeUnion*> (sh);
    const GeoShape* shA = sub->getOpA();
    const GeoShape* shB = sub->getOpB();
    splitShape(shA,shapes);
    splitShape(shB,shapes);
  } else {
    shapes.push_back(sh);
  }
  return;
}

void  Muon::MuonInertMaterialBuilder::getVolumeFractions( ) const
{
  m_volFractions.clear();
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTBevelledLongTubeIn", std::pair<double,double> (1.20078e+09,0.9367) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTBevelledLongTubeOut", std::pair<double,double> (1.20078e+09,0.9566) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTVoussoirAttWing", std::pair<double,double> (2.37221e+09,0.0073) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTBevelledShortTube", std::pair<double,double> (1.78598e+08,0.8698) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTBevelledCornerTube", std::pair<double,double> (6.34357e+07,0.6432) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTRibEnvelope", std::pair<double,double> (5.09211e+08,0.1497) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTWingRib", std::pair<double,double> (5.60805e+08,0.0479) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTVoussoirAttachment", std::pair<double,double> (2.036e+08,0.031) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTColdLongSegment", std::pair<double,double> (5.95526e+09,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTColdShortSegment", std::pair<double,double> (8.37822e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTColdCornerSegment", std::pair<double,double> (1.69662e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTColdRib", std::pair<double,double> (2.4552e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTVoussoir", std::pair<double,double> (1.47921e+09,0.307) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("EdgeBTVoussoir", std::pair<double,double> (2.0092e+09,0.0509) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("HeadBTVoussoir", std::pair<double,double> (3.54772e+08,0.5719) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTStrut", std::pair<double,double> (2.0905e+09,0.1766) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTWingStrut", std::pair<double,double> (1.56065e+09,0.0489) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("BTCryoring", std::pair<double,double> (2.00278e+08,0.9803) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("Rail", std::pair<double,double> (8.01099e+09,0.4905) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTEndplateStdSegment", std::pair<double,double> (9.42675e+08,0.8750) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTEndplateRailSegment6", std::pair<double,double> (9.42675e+08,0.8748) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTEndplateRailSegment8", std::pair<double,double> (9.42675e+08,0.8736) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTWallStdSegment", std::pair<double,double> (5.59807e+10,0) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTWallRailSegment", std::pair<double,double> (5.59807e+10,0) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTCentralTube", std::pair<double,double> (2.90421e+09,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTStaytube", std::pair<double,double> (2.44805e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTConductorBox", std::pair<double,double> (5.3088e+09,0.9619) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTKeystoneBox", std::pair<double,double> (1.91453e+10,0.0569) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTTower", std::pair<double,double> (1.55166e+11,0.0422) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTBottomTower", std::pair<double,double> (2.06925e+08,0.7120) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTBackTowerWall", std::pair<double,double> (9.25155e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ECTServiceTurretTower", std::pair<double,double> (4.72478e+10,0.0108) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetSidePlate", std::pair<double,double> (2.72438e+09,0.206) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetLargePlate", std::pair<double,double> (1.31733e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetCoverPlate", std::pair<double,double> (2.24968e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetInnerPlate", std::pair<double,double> (5.34198e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetBottomPlate", std::pair<double,double> (1.1682e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetSidePlate", std::pair<double,double> (2.83923e+09,0.1966) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetLargePlate", std::pair<double,double> (1.38125e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetInclinedPlate", std::pair<double,double> (7.11822e+07,0.8284) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetCoverPlate", std::pair<double,double> (3.02684e+08,0.9619) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetInnerPlate", std::pair<double,double> (7.41336e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetInnerTopVoussPlate", std::pair<double,double> (5.695e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetBottomPlate", std::pair<double,double> (1.57176e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetGirder01", std::pair<double,double> (9.98861e+08,0.1605) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetGirder02", std::pair<double,double> (1.26939e+09,0.1545) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetGirder03", std::pair<double,double> (1.40019e+09,0.1517) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetRailSupport", std::pair<double,double> (5.84777e+08,0.5076) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetMinusRailSupport", std::pair<double,double> (8.1122e+08,0.4129) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ExtrFeetPlusRailSupport", std::pair<double,double> (7.9211e+08,0.4238) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("StdFeetVoussoir", std::pair<double,double> (1.26638e+09,0.3078) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ConnFeetVoussoir", std::pair<double,double> (1.9295e+08,0.4927) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingFrontDisk", std::pair<double,double> (1.05675e+09,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingBackDisk", std::pair<double,double> (4.89468e+09,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingPlugsExtension", std::pair<double,double> (9.07988e+06,0) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingPlugs", std::pair<double,double> (1.21083e+09,0.8779) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingTubeBackDisk", std::pair<double,double> (7.18286e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingMainTube", std::pair<double,double> (6.46201e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingBrassCone", std::pair<double,double> (1.17168e+08,0.4441) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingPolyCone", std::pair<double,double> (2.86637e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingLeadCone", std::pair<double,double> (1.2051e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingHubBrassCone", std::pair<double,double> (7.90909e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingPolyCladding", std::pair<double,double> (2.59232e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("DiskShieldingLeadCladding", std::pair<double,double> (1.17697e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ToroidShieldingOuterPlugs", std::pair<double,double> (4.5584e+09,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ToroidShieldingInnerPlugs", std::pair<double,double> (3.47561e+09,0.8579) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ToroidShieldingPolyRingsTube", std::pair<double,double> (9.50078e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ToroidShieldingPolyRingsCone", std::pair<double,double> (1.28304e+07,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ForwardShieldingPlug", std::pair<double,double> (1.50796e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ForwardShieldingMainCylinder", std::pair<double,double> (5.25941e+10,0.9762) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ForwardShieldingOctogon", std::pair<double,double> (4.86944e+10,0.4814) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ForwardShieldingTX1STMTube", std::pair<double,double> (8.63938e+08,1) ) );
m_volFractions.push_back(std::pair<std::string,std::pair<double,double> > ("ForwardShieldingTX1STMCone", std::pair<double,double> (6.28319e+08,1) ) );
  return;
}

void Muon::MuonInertMaterialBuilder::removeTV( const Trk::Volume* vol ) const
{ 
  const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::SubtractedVolumeBounds* sub  = dynamic_cast<const Trk::SubtractedVolumeBounds*> (&(vol->volumeBounds()));

  if (comb) {
    removeTV(comb->first());
    removeTV(comb->second());
  } else if (sub) {
    removeTV(sub->outer());
    removeTV(sub->inner());
  } else {
    delete vol;
  }
  return;
}

double Muon::MuonInertMaterialBuilder::thinDim( const Trk::Volume* vol) const
{ 
  // returns thinnest dimension
  double dim = 50;
 
  const Trk::CylinderVolumeBounds* cylBounds = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(vol->volumeBounds()));
  const Trk::CuboidVolumeBounds*   cubBounds = dynamic_cast<const Trk::CuboidVolumeBounds*>   (&(vol->volumeBounds()));
  const Trk::TrapezoidVolumeBounds* trdBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(vol->volumeBounds()));
	  
  if (cylBounds) {
    double z = cylBounds->halflengthZ(); 
    double ri = cylBounds->innerRadius(); 
    double ro = cylBounds->outerRadius(); 
    if ( z <= (ro-ri) ) return z;
    else return (ro-ri);
  }

  if (cubBounds ) {
    double x =  cubBounds->halflengthX();
    double y =  cubBounds->halflengthY(); 
    double z = 	cubBounds->halflengthZ();

    if ( z<=x && z<=y ) return z;
    else if ( y<=x && y<=z ) return y;
    else return x;
  }

  if (trdBounds ) {
    double y =  trdBounds->halflengthY(); 
    double z = 	trdBounds->halflengthZ();
    double minX = trdBounds->minHalflengthX();
    double maxX = trdBounds->maxHalflengthX();

    if ( z<=y ) return z;
    else return y;
  }
   
  return dim;
}

double  Muon::MuonInertMaterialBuilder::calculateVolume( const Trk::Volume* envelope) const
{
  double envVol = 0.;
  
  if (!envelope) return 0.;
  
  const Trk::CylinderVolumeBounds*  cyl = dynamic_cast<const Trk::CylinderVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::CuboidVolumeBounds*    box = dynamic_cast<const Trk::CuboidVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*> (&(envelope->volumeBounds()));
  const Trk::BevelledCylinderVolumeBounds*  bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*> (&(envelope->volumeBounds()));
  
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
	std::cout << "new volume added:"<< clv->getName() <<","<<clv->getMaterial()->getName()<<std::endl;
	printInfo(cv);
      }
    } else {
      getObjsForTranslation(cv, transform*transf, vols);
    }
  }
}
