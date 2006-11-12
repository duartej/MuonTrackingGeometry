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
// Gaudi
#include "GaudiKernel/AlgTool.h"
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

      const std::vector<const Trk::DetachedTrackingVolume*>* buildDetachedTrackingVolumes() const; 
      const std::vector<const Trk::TrackingVolume*>* buildDetachedTrackingVolumeTypes() const; 

    private:

//        void checkObject( const Trk::TrackingVolume* objs) const;
      void checkObject( const std::vector<const Trk::DetachedTrackingVolume*>*trkVolumes) const;
      const void printInfo(const GeoVPhysVol* pv) const;
      const void printChildren(const GeoVPhysVol* pv) const;
      const void decodeShape(const GeoShape* sh) const;
      Trk::Volume* translateGeoShape(const GeoShape* sh, HepTransform3D* tr) const;

      Trk::TrapezoidVolumeBounds* decodeColdSegment(const GeoShape* sh) const;
      Trk::VolumeBounds* decodeECTSegment(const GeoShape* sh) const;
      Trk::BevelledCylinderVolumeBounds* decodeBevelledCylinder(const GeoShape* sh) const;

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      Trk::IMagneticFieldTool*            m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      std::string                         m_magFieldToolName;            //!< Name of the Tracking Magnetic Field Svc
      std::string                         m_magFieldToolInstanceName;    //!< Instance Name of Tracking Magnetic Field Svc
      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the material
      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      Trk::MagneticFieldProperties        m_muonMagneticField;          //!< the magnetic Field

//mw
      Trk::GeoMaterialConverter*           m_materialConverter;          //!< material converter
    };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONINERTMATERIALBUILDER_H


