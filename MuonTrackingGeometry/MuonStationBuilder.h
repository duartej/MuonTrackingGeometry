///////////////////////////////////////////////////////////////////
// MuonStationBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONSTATIONBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONSTATIONBUILDER_H

//Trk
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkEventPrimitives/GlobalPosition.h"
#include "TrkEventPrimitives/GlobalMomentum.h"
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkGeometry/MaterialProperties.h"
#include "TrkGeometry/DetachedTrackingVolume.h"
#include "TrkGeometry/TrackingVolume.h"
// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GeoModelKernel/GeoVPhysVol.h"

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
}

namespace MuonGM {
  class MuonDetectorManager;
  class MuonStation;
}
 
namespace Muon {

  class MuonStationTypeBuilder;      
  
  /** @class MuonStationBuilder
  
      The Muon::MuonStationBuilder retrieves muon stations from Muon Geometry Tree
      
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

      const MuonGM::MuonDetectorManager*  m_muonMgr;               //!< the MuonDetectorManager
      std::string                         m_muonMgrLocation;       //!< the location of the Muon Manager
      Muon::MuonStationTypeBuilder*       m_muonStationTypeBuilder;             //!< Helper Tool to create TrackingVolume Arrays
      std::string                         m_muonStationTypeBuilderName;         //!< Name of the helper tool
      std::string                         m_muonStationTypeBuilderInstanceName; //!< Instance of the helper tool
      Trk::IMagneticFieldTool*            m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      std::string                         m_magFieldToolName;            //!< Name of the Tracking Magnetic Field Svc
      std::string                         m_magFieldToolInstanceName;    //!< Instance Name of Tracking Magnetic Field Svc
      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the material
      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      Trk::MagneticFieldProperties        m_muonMagneticField;          //!< the magnetic Field
    };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONSTATIONBUILDER_H


