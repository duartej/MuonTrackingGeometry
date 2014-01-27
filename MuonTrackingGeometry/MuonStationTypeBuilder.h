//////////////////////////////////////////////////////////////////
// MuonStationTypeBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H

//Trk
#include "TrkGeometry/TrackingVolume.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkDetDescrGeoModelCnv/GeoMaterialConverter.h"
// Gaudi
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#include "GeoModelKernel/GeoVPhysVol.h"
#include "GeoModelKernel/GeoMaterial.h"

namespace Trk {
 class TrackingGeometry;
 class TrackingVolume;
 class DetachedTrackingVolume;
 class MaterialProperties;
 class Volume;
 class Layer;
 class ITrackingVolumeArrayCreator;
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

  struct LayerMaterial{
    Trk::MaterialProperties material;
    float                   thickness;

    LayerMaterial(){
      material = Trk::MaterialProperties(10.e10,10.e10,0.,0.,0.);
      thickness = 0.;
    }

    void addMaterial(Trk::MaterialProperties mat, float d) {
  
      thickness += d;
      float fnew = d/thickness;
      float fold = 1.-fnew;
      
      material = Trk::MaterialProperties(1./(fnew/mat.x0()+fold/material.x0()),
					 1./(fnew/mat.l0()+fold/material.l0()),
				         fnew*mat.averageA()+fold*material.averageA(),
				         fnew*mat.averageZ()+fold*material.averageZ(),
				         fnew*mat.averageRho()+fold*material.averageRho());
        
    }

    void scale( float factor ) {     // scaling of material thickness

      thickness *= factor;
      material   = Trk::MaterialProperties( factor*material.x0(), factor*material.l0(),
                                            material.averageA(), material.averageZ(),
                                            material.averageRho()/factor );        
    }

    LayerMaterial operator=( const LayerMaterial layMat) {
      material = layMat.material;
      thickness = layMat.thickness;
      return (*this);
    } 

    bool empty() { return !(thickness>0.) ;  }

  };

  typedef std::pair<Trk::SharedObject<const Trk::Layer>,const Amg::Transform3D*> LayTr;
     
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
      const Trk::TrackingVolume* processMdtBox(Trk::Volume*&,const GeoVPhysVol*&, Amg::Transform3D*, double ) const;
      const Trk::TrackingVolume* processMdtTrd(Trk::Volume*&,const GeoVPhysVol*&, Amg::Transform3D*) const;
      const Trk::TrackingVolume* processRpc(Trk::Volume*&,std::vector<const GeoVPhysVol*>,std::vector<Amg::Transform3D>) const;
      const Trk::TrackingVolume* processSpacer(Trk::Volume*&,std::vector<const GeoVPhysVol*>,std::vector<Amg::Transform3D>) const;
      const Trk::TrackingVolume* processNSW(std::vector<const Trk::Layer*>) const;
      const Trk::LayerArray* processCSCTrdComponent(const GeoVPhysVol*&, Trk::TrapezoidVolumeBounds*&, Amg::Transform3D*&) const;
      const Trk::LayerArray* processCSCDiamondComponent(const GeoVPhysVol*&, Trk::DoubleTrapezoidVolumeBounds*&, Amg::Transform3D*&) const;
      const Trk::LayerArray* processTGCComponent(const GeoVPhysVol*&, Trk::TrapezoidVolumeBounds*&, Amg::Transform3D*&) const;
      std::pair<const Trk::Layer*,const std::vector<const Trk::Layer*>* > createLayerRepresentation(const Trk::TrackingVolume* trVol) const; 
      const Trk::Layer* createLayer(const Trk::TrackingVolume* trVol,Trk::MaterialProperties*, Amg::Transform3D&) const; 
      Identifier identifyNSW( std::string, Amg::Transform3D ) const;
   
      void printChildren(const GeoVPhysVol*) const ;
      // used to be private ..
      double get_x_size(const GeoVPhysVol*) const ;
      double decodeX(const GeoShape*) const ;
      double getVolume(const GeoShape*) const;
      LayerMaterial getAveragedLayerMaterial(const GeoVPhysVol*,double,double) const;
      void collectMaterial(const GeoVPhysVol*,LayerMaterial&,double) const;
      LayerMaterial collectStationMaterial(const Trk::TrackingVolume* trVol,double) const; 

  private:
     /** Private method to fill default material */
      //void fillDefaultServiceMaterial();

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      bool                                m_multilayerRepresentation;   
      bool                                m_resolveSpacer;   

      ToolHandle<Trk::ITrackingVolumeArrayCreator>   m_trackingVolumeArrayCreator;  //!< Helper Tool to create TrackingVolume Arrays

      mutable Trk::MaterialProperties*     m_muonMaterial;               //!< the material
 
      mutable LayerMaterial    m_mdtTubeMat;                  //
      mutable std::vector<LayerMaterial>   m_mdtFoamMat;                  //
      mutable LayerMaterial    m_rpc46;                  
      mutable std::vector<LayerMaterial>   m_rpcDed;                
      mutable LayerMaterial    m_rpcLayer;                  
      mutable LayerMaterial    m_rpcExtPanel;                  
      mutable LayerMaterial    m_rpcMidPanel;                  
      mutable LayerMaterial    m_matCSC01;                  //
      mutable LayerMaterial    m_matCSCspacer1;                  //
      mutable LayerMaterial    m_matCSC02;                  //
      mutable LayerMaterial    m_matCSCspacer2;                  //
      mutable LayerMaterial    m_matTGC01;                  //
      mutable LayerMaterial    m_matTGC06;                  //
      Trk::GeoMaterialConverter*                  m_materialConverter;

    };

} // end of namespace

#endif //MUONTRACKINGGEOMETRY_MUONSTATIONTYPEBUILDER_H
