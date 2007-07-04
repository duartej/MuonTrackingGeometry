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

      const std::vector<const Trk::DetachedTrackingVolume*>* buildDetachedTrackingVolumes() const; 
      const std::vector<const Trk::DetachedTrackingVolume*>* buildDetachedTrackingVolumeTypes() const; 

    private:

//        void checkObject( const Trk::TrackingVolume* objs) const;
      void checkObject( const std::vector<const Trk::DetachedTrackingVolume*>*trkVolumes) const;
      const void printInfo(const GeoVPhysVol* pv) const;
      const void printChildren(const GeoVPhysVol* pv) const;
      const void decodeShape(const GeoShape* sh) const;
      Trk::Volume* translateGeoShape(const GeoShape* sh, HepTransform3D* tr) const;
      const Trk::TrackingVolume* simplifyShape(const Trk::TrackingVolume* tr) const;
      const Trk::Volume* findEnvelope(const Trk::TrackingVolume* tr) const;
      const Trk::TrackingVolume* findECTEnvelope(const std::vector<const Trk::TrackingVolume*>* vols) const;

      Trk::TrapezoidVolumeBounds* decodeColdSegment(const GeoShape* sh) const;
      Trk::VolumeBounds* decodeECTSegment(const GeoShape* sh) const;
      Trk::BevelledCylinderVolumeBounds* decodeBevelledCylinder(const GeoShape* sh) const;

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      bool                                m_simplify;              // switch geometry simplification on/off 
      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the material
      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system
 
      ToolHandle<Trk::IMagneticFieldTool> m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      mutable Trk::MagneticFieldProperties m_muonMagneticField;          //!< the magnetic Field

//mw
      Trk::GeoMaterialConverter*           m_materialConverter;          //!< material converter
    };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONINERTMATERIALBUILDER_H


