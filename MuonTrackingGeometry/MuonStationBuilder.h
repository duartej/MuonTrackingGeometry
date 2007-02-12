///////////////////////////////////////////////////////////////////
// MuonStationBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONSTATIONBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONSTATIONBUILDER_H

//Trk
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkGeometry/MaterialProperties.h"
#include "TrkGeometry/DetachedTrackingVolume.h"
#include "TrkGeometry/TrackingVolume.h"
// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GeoModelKernel/GeoVPhysVol.h"

class MdtIdHelper;
class RpcIdHelper;
class CscIdHelper;
class TgcIdHelper;

namespace Trk {
 class TrackingGeometry;
 class TrackingVolume;
 class Volume;
 class Layer;
 class ITrackingVolumeBuilder;
 class ITrackingVolumeArrayCreator;
 class ITrackingVolumeHelper;
 class ILayerBuilder;
 class ILayerArrayCreator;
 class IMagneticFieldTool;
}

namespace MuonGM {
  class MuonDetectorManager;
  class MuonStation;
}
 
namespace Muon {

  class MuonStationTypeBuilder;      
  
  /** @class MuonStationBuilder
  
      The Muon::MuonStationBuilder retrieves muon stations from Muon Geometry Tree
      prototypes built with help of Muon::MuonStationTypeBuilder
      
      by Sarka.Todorova@cern.ch
    */
    
  class MuonStationBuilder : public AlgTool,
                             virtual public Trk::IDetachedTrackingVolumeBuilder {
  public:
      /** Constructor */
      MuonStationBuilder(const std::string&,const std::string&,const IInterface*);
      /** Destructor */
      ~MuonStationBuilder();
      /** AlgTool initailize method.*/
      StatusCode initialize();
      /** AlgTool finalize method */
      StatusCode finalize();

      const std::vector<const Trk::DetachedTrackingVolume*>* buildDetachedTrackingVolumes() const; 
      const std::vector<const Trk::TrackingVolume*>* buildDetachedTrackingVolumeTypes() const; 

    private:

      void glueComponents(const Trk::DetachedTrackingVolume* ) const;    
      void encloseLayers( const Trk::DetachedTrackingVolume* ) const; 
      void identifyLayers(const Trk::DetachedTrackingVolume*, int, int ) const;
 
      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      const MdtIdHelper*            m_mdtIdHelper;           //!< 
      const RpcIdHelper*            m_rpcIdHelper;           //!< 
      const CscIdHelper*            m_cscIdHelper;           //!< 
      const TgcIdHelper*            m_tgcIdHelper;           //!< 
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager

      Trk::IMagneticFieldTool*            m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      std::string                         m_magFieldToolName;            //!< Name of the Tracking Magnetic Field Svc
      std::string                         m_magFieldToolInstanceName;    //!< Instance Name of Tracking Magnetic Field Svc

      Muon::MuonStationTypeBuilder*       m_muonStationTypeBuilder;             //!< Helper Tool to create TrackingVolume Arrays
      std::string                         m_muonStationTypeBuilderName;         //!< Name of the helper tool
      std::string                         m_muonStationTypeBuilderInstanceName; //!< Instance of the helper tool

      Trk::ITrackingVolumeHelper*         m_trackingVolumeHelper;             //!< Helper Tool to create TrackingVolumes
      std::string                         m_trackingVolumeHelperName;         //!< Name of the helper tool
      std::string                         m_trackingVolumeHelperInstanceName; //!< Instance of the helper tool

      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the material
      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      Trk::MagneticFieldProperties        m_muonMagneticField;          //!< the magnetic Field
      bool                                m_buildBarrel;
      bool                                m_buildEndcap;
      bool                                m_buildCsc;
      bool                                m_buildTgc;  
      bool                                m_identifyActive;
    };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONSTATIONBUILDER_H


