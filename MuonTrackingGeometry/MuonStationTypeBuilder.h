//////////////////////////////////////////////////////////////////
// MuonStationTypeBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H

//Trk
#include "TrkGeometry/MagneticFieldProperties.h"
//#include "TrkGeometry/ExtendedMaterialProperties.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkDetDescrGeoModelCnv/GeoMaterialConverter.h"
// Gaudi
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#include "GeoModelKernel/GeoVPhysVol.h"
#include "GeoModelKernel/GeoMaterial.h"
// CLHEP
#include "CLHEP/Geometry/Transform3D.h"

namespace Trk {
 class TrackingGeometry;
 class TrackingVolume;
 class DetachedTrackingVolume;
 class ExtendedMaterialProperties;
 class LayerMaterialProperties;
 class Volume;
 class Layer;
 class ITrackingVolumeArrayCreator;
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

  typedef std::pair<Trk::SharedObject<const Trk::Layer>,const HepGeom::Transform3D*> LayTr;
     
  /** @class MuonStationTypeBuilder
  
      The Muon::MuonStationTypeBuilder retrieves components of muon stations from Muon Geometry Tree,
      builds 'prototype' object (TrackingVolume with NameType)
      
      by Sarka.Todorova@cern.ch
    */
    
  class MuonStationTypeBuilder : public AthAlgTool{
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
      const Trk::TrackingVolume* processMdtBox(Trk::Volume*&,const GeoVPhysVol*&, HepGeom::Transform3D*, double ) const;
      const Trk::TrackingVolume* processMdtTrd(Trk::Volume*&,const GeoVPhysVol*&, HepGeom::Transform3D*) const;
      const Trk::TrackingVolume* processRpc(Trk::Volume*&,std::vector<const GeoVPhysVol*>,std::vector<HepGeom::Transform3D>) const;
      const Trk::TrackingVolume* processSpacer(Trk::Volume*&,std::vector<const GeoVPhysVol*>,std::vector<HepGeom::Transform3D>) const;
      const Trk::TrackingVolume* processNSW(std::vector<const Trk::Layer*>) const;
      const Trk::LayerArray* processCSCTrdComponent(const GeoVPhysVol*&, Trk::TrapezoidVolumeBounds*&, HepGeom::Transform3D*&) const;
      const Trk::LayerArray* processCSCDiamondComponent(const GeoVPhysVol*&, Trk::DoubleTrapezoidVolumeBounds*&, HepGeom::Transform3D*&) const;
      const Trk::LayerArray* processTGCComponent(const GeoVPhysVol*&, Trk::TrapezoidVolumeBounds*&, HepGeom::Transform3D*&) const;
      std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>* > createLayerRepresentation(const Trk::TrackingVolume* trVol) const; 
      const Trk::Layer* createLayer(const Trk::TrackingVolume* trVol,Trk::ExtendedMaterialProperties*, HepGeom::Transform3D&) const; 
      Identifier identifyNSW( std::string, HepGeom::Transform3D ) const;
   
      void printChildren(const GeoVPhysVol*) const ;
      // used to be private ..
      double get_x_size(const GeoVPhysVol*) const ;
      double decodeX(const GeoShape*) const ;
      double getVolume(const GeoShape*) const;
      Trk::ExtendedMaterialProperties* getAveragedLayerMaterial(const GeoVPhysVol*,double,double) const;
      void collectMaterial(const GeoVPhysVol*,Trk::ExtendedMaterialProperties*& ,double) const;
      Trk::ExtendedMaterialProperties collectStationMaterial(const Trk::TrackingVolume* trVol,double) const; 

  private:
     /** Private method to fill default material */
      //void fillDefaultServiceMaterial();

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      bool                                m_multilayerRepresentation;   
      bool                                m_resolveSpacer;   

      ToolHandle<Trk::ITrackingVolumeArrayCreator>   m_trackingVolumeArrayCreator;  //!< Helper Tool to create TrackingVolume Arrays
      unsigned int                                   m_magFieldMode;                //!< set the magnetic field mode
      ToolHandle<Trk::IMagneticFieldTool>            m_magFieldTool;                //!< Tracking Interface to Magnetic Field

      mutable Trk::MaterialProperties*     m_muonMaterial;               //!< the material
      mutable Trk::MagneticFieldProperties m_muonMagneticField;          //!< the magnetic Field
 
      mutable Trk::ExtendedMaterialProperties*    m_mdtTubeMat;                  //
      mutable std::vector<Trk::ExtendedMaterialProperties*>    m_mdtFoamMat;                  //
      mutable Trk::ExtendedMaterialProperties*    m_rpc46;                  
      mutable std::vector<Trk::ExtendedMaterialProperties*>    m_rpcDed;                
      mutable Trk::ExtendedMaterialProperties*    m_rpcLayer;                  
      mutable Trk::ExtendedMaterialProperties*    m_rpcExtPanel;                  
      mutable Trk::ExtendedMaterialProperties*    m_rpcMidPanel;                  
      mutable Trk::ExtendedMaterialProperties*    m_matCSC01;                  //
      mutable Trk::ExtendedMaterialProperties*    m_matCSCspacer1;                  //
      mutable Trk::ExtendedMaterialProperties*    m_matCSC02;                  //
      mutable Trk::ExtendedMaterialProperties*    m_matCSCspacer2;                  //
      mutable Trk::ExtendedMaterialProperties*    m_matTGC01;                  //
      mutable Trk::ExtendedMaterialProperties*    m_matTGC06;                  //
      Trk::GeoMaterialConverter*                  m_materialConverter;

    };

} // end of namespace

#endif //MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H
