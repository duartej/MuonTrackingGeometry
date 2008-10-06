///////////////////////////////////////////////////////////////////
// MuonInertMaterialBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONINERTMATERIALBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONINERTMATERIALBUILDER_H

//Trk
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkGeometry/MaterialProperties.h"
#include "TrkGeometry/DetachedTrackingVolume.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkGeometry/TrackingVolumeManipulator.h"
#include "TrkDetDescrGeoModelCnv/GeoShapeConverter.h"
// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#include "GeoModelKernel/GeoVPhysVol.h"
//mw
#include "TrkDetDescrGeoModelCnv/GeoMaterialConverter.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/BevelledCylinderVolumeBounds.h"
#include "TrkVolumes/CylinderVolumeBounds.h"

namespace Trk {
 class TrackingGeometry;
 class TrackingVolume;
 class Volume;
 class Layer;
 class ITrackingVolumeBuilder;
 class ITrackingVolumeArrayCreator;
 class ILayerBuilder;
 class ILayerArrayCreator;
 class IMagneticFieldTool;

//mw
 class MaterialProperties;

}

namespace MuonGM {
  class MuonDetectorManager;
  class MuonStation;
}
 
namespace Muon {
        
  /** @class MuonInertMaterialBuilder
  
      The Muon::MuonInertMaterialBuilder retrieves muon stations from Muon Geometry Tree
      
      by Sarka.Todorova@cern.ch, Marcin.Wolter@cern.ch
    */
    
  class MuonInertMaterialBuilder : public AlgTool,
                                   public Trk::TrackingVolumeManipulator,
                             virtual public Trk::IDetachedTrackingVolumeBuilder {
  public:
      /** Constructor */
      MuonInertMaterialBuilder(const std::string&,const std::string&,const IInterface*);
      /** Destructor */
      ~MuonInertMaterialBuilder();
      /** AlgTool initailize method.*/
      StatusCode initialize();
      /** AlgTool finalize method */
      StatusCode finalize();

      /** Method returning cloned and positioned material objects */
      const std::vector<const Trk::DetachedTrackingVolume*>* buildDetachedTrackingVolumes() const; 

    private:

      /** Method creating material object prototypes */
      const std::vector<std::pair<const Trk::DetachedTrackingVolume*,std::vector<HepTransform3D> > >* buildDetachedTrackingVolumeTypes() const; 
      /** Method extracting material objects from GeoModel tree */
      void getObjsForTranslation(const GeoVPhysVol* pv,HepTransform3D , std::vector<std::pair<const GeoLogVol*,std::vector<HepTransform3D> > >& vols ) const;
      /** Dump from GeoModel tree  */
      const void printInfo(const GeoVPhysVol* pv) const;
      const void printChildren(const GeoVPhysVol* pv) const;
      /** Simplification of GeoModel object + envelope */
      const Trk::TrackingVolume* simplifyShape(const Trk::TrackingVolume* tr) const;
      
      const Trk::Layer* boundarySurfaceToLayer(const Trk::Surface&, const Trk::MaterialProperties*, double) const;
      Trk::Volume* createSubtractedVolume(const HepTransform3D& tr, Trk::Volume* subtrVol) const;
      void  splitShape(const GeoShape* sh, std::vector<const GeoShape*>& shapes) const;

      std::pair<std::vector<const Trk::Layer*>*, std::vector<const Trk::TrackingVolume*>* >
	translateToLayers(const std::vector<const Trk::TrackingVolume*>* vols, int mode) const;
      std::vector<const Trk::Layer*>*  translateBoundariesToLayers(const Trk::Volume* vol, const Trk::TrackingVolume* trVol, double) const;
      double volumeToLayers(std::vector<const Trk::Layer*>& lays, const Trk::Volume* vol, 
			  Trk::Volume* subtrVol, const Trk::MaterialProperties* mat) const;
      const bool checkVolume(const Trk::Volume*) const;
      void getVolumeFractions() const;
      void removeTV(const Trk::Volume*) const;
      double thinDim( const Trk::Volume* vol) const;
      double calculateVolume( const Trk::Volume* envelope) const;

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      bool                                m_simplify;              // switch geometry simplification on/off 
      bool                                m_simplifyToLayers;      // switch geometry simplification to layers on/off 
      double                              m_layerThicknessLimit;   // maximal thickness of layer allowed (in X0)   
      bool                                m_debugMode;             // build layers & dense volumes in parallel - double counting material !!! 
      bool                                m_buildBT;                    // build barrel toroids 
      bool                                m_buildECT;                   // build endcap toroids 
      bool                                m_buildFeets;                 // build feets 
      int                                 m_buildRails;                 // build rails 
      bool                                m_buildShields;               // build shieldings 
      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the material
      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system
 
      ToolHandle<Trk::IMagneticFieldTool> m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      mutable Trk::MagneticFieldProperties m_muonMagneticField;          //!< the magnetic Field
//mw
      Trk::GeoMaterialConverter*           m_materialConverter;          //!< material converter
      Trk::GeoShapeConverter*              m_geoShapeConverter;          //!< shape converter

      mutable std::vector<std::pair<std::string,std::pair<double,double> > >   m_volFractions;
      
      mutable std::vector<Trk::Volume*>    m_garbage;      
    };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONINERTMATERIALBUILDER_H


