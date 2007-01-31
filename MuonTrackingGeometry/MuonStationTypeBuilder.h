///////////////////////////////////////////////////////////////////
// MuonStationTypeBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H

//Trk
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkGeometry/MaterialProperties.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkDetDescrGeoModelCnv/GeoMaterialConverter.h"
// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GeoModelKernel/GeoVPhysVol.h"
#include "GeoModelKernel/GeoMaterial.h"
// CLHEP
#include "CLHEP/Geometry/Transform3D.h"

namespace Trk {
 class TrackingGeometry;
 class TrackingVolume;
 class DetachedTrackingVolume;
 class MaterialProperties;
 class LayerMaterialProperties;
 class Volume;
 class Layer;
 class ITrackingVolumeBuilder;
 class ITrackingVolumeArrayCreator;
 class ILayerBuilder;
 class ILayerArrayCreator;
 class IMagneticFieldTool;
 class CuboidVolumeBounds;
 class TrapezoidVolumeBounds;
 class DoubleTrapezoidVolumeBounds;
 class PlaneLayer;   
}

namespace MuonGM {
  class MuonDetectorManager;
  class MuonStation;
}
 
namespace Muon {

  typedef std::pair<Trk::SharedObject<const Trk::Layer>,const HepTransform3D*> LayTr;
     
  /** @class MuonStationTypeBuilder
  
      The Muon::MuonStationTypeBuilder retrieves components of muon stations from Muon Geometry Tree,
      builds 'prototype' object (TrackingVolume with NameType)
      
      by Sarka.Todorova@cern.ch
    */
    
  class MuonStationTypeBuilder : public AlgTool{
  public:
      /** Constructor */
      MuonStationTypeBuilder(const std::string&,const std::string&,const IInterface*);
      /** Destructor */
      virtual ~MuonStationTypeBuilder();
      /** AlgTool initailize method.*/
      StatusCode initialize();
      /** AlgTool finalize method */
      StatusCode finalize();
      /** Interface methode */
      static const InterfaceID& interfaceID();
      /** steering routine */
      const Trk::TrackingVolumeArray* processBoxStationComponents(const GeoVPhysVol* cv, Trk::CuboidVolumeBounds* envBounds) const; 
      const Trk::TrackingVolumeArray* processTrdStationComponents(const GeoVPhysVol* cv, Trk::TrapezoidVolumeBounds* envBounds) const; 
      const Trk::TrackingVolume* processCscStation(const GeoVPhysVol* cv, std::string name) const; 
      std::vector<const Trk::TrackingVolume*> processTgcStation(const GeoVPhysVol* cv) const; 
      /** components */
      const Trk::TrackingVolume* processMdtBox(Trk::Volume*&,const GeoVPhysVol*&, HepTransform3D*&) const;
      const Trk::TrackingVolume* processMdtTrd(Trk::Volume*&,const GeoVPhysVol*&, HepTransform3D*&) const;
      const Trk::TrackingVolume* processRpc(Trk::Volume*&,std::vector<const GeoVPhysVol*>,std::vector<HepTransform3D*>) const;
      const Trk::TrackingVolume* processSpacer(Trk::Volume*&,std::vector<const GeoVPhysVol*>,std::vector<HepTransform3D*>) const;
      const Trk::LayerArray* processCSCTrdComponent(const GeoVPhysVol*&, Trk::TrapezoidVolumeBounds*&, HepTransform3D*&) const;
      const Trk::LayerArray* processCSCDiamondComponent(const GeoVPhysVol*&, Trk::DoubleTrapezoidVolumeBounds*&, HepTransform3D*&) const;
      const Trk::LayerArray* processTGCComponent(const GeoVPhysVol*&, Trk::TrapezoidVolumeBounds*&, HepTransform3D*&) const;
      const Trk::MaterialProperties getLayerMaterial(std::string mat, double thickness) const;
      const Trk::Layer* createLayerRepresentation(const Trk::TrackingVolume* trVol) const; 
   
      const void printChildren(const GeoVPhysVol*) const ;
  private:
      const double get_x_size(const GeoVPhysVol*) const ;
      const double decodeX(const GeoShape*) const ;
      Trk::MaterialProperties* getAveragedLayerMaterial(const GeoVPhysVol*,double,double) const;
      Trk::MaterialProperties collectMaterial(const GeoVPhysVol*,Trk::MaterialProperties ,double) const;
      /** Private method to fill default material */
      //void fillDefaultServiceMaterial();

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      Trk::ITrackingVolumeArrayCreator*   m_trackingVolumeArrayCreator;             //!< Helper Tool to create TrackingVolume Arrays
      std::string                         m_trackingVolumeArrayCreatorName;         //!< Name of the helper tool
      std::string                         m_trackingVolumeArrayCreatorInstanceName; //!< Instance of the helper tool
      Trk::IMagneticFieldTool*            m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      std::string                         m_magFieldToolName;            //!< Name of the Tracking Magnetic Field Svc
      std::string                         m_magFieldToolInstanceName;    //!< Instance Name of Tracking Magnetic Field Svc
      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the material
      Trk::MagneticFieldProperties        m_muonMagneticField;          //!< the magnetic Field
      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      mutable Trk::MaterialProperties*    m_mdtTubeMat;                  //
      mutable Trk::MaterialProperties*    m_mdtFoamMat;                  //
      mutable Trk::MaterialProperties*    m_rpc46;                  
      mutable Trk::MaterialProperties*    m_rpcDed50;                
      mutable Trk::MaterialProperties*    m_rpcLayer;                  
      mutable Trk::MaterialProperties*    m_rpcExtPanel;                  
      mutable Trk::MaterialProperties*    m_rpcMidPanel;                  
      mutable Trk::MaterialProperties*    m_matCSC01;                  //
      mutable Trk::MaterialProperties*    m_matCSCspacer1;                  //
      mutable Trk::MaterialProperties*    m_matCSC02;                  //
      mutable Trk::MaterialProperties*    m_matCSCspacer2;                  //
      mutable Trk::MaterialProperties*    m_matTGC01;                  //
      mutable Trk::MaterialProperties*    m_matTGC06;                  //
      Trk::GeoMaterialConverter*          m_materialConverter;

    };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H


