///////////////////////////////////////////////////////////////////
// MuonTrackingGeometryBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H

//Trk
#include "TrkDetDescrInterfaces/IGeometryBuilder.h"
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkGeometry/MaterialProperties.h"
// Gaudi
#include "GaudiKernel/AlgTool.h"

namespace Trk {
 class TrackingGeometry;
 class VolumeBounds;
 class ITrackingVolumeBuilder;
 class ITrackingVolumeHelper;
 class ITrackingVolumeDisplayer;
 class ITrackingVolumeArrayCreator;
 class IMagneticFieldTool;
}
 
namespace Muon {
     
  /** @class MuonTrackingGeometryBuilder
  
      The Muon::MuonTrackingGeometryBuilder retrieves LayerBuilders and VolumeBuilders
      for the Muon Detector sub detectors and combines the given Volumes to a full Trk::TrackingGeometry.
      
      Inheriting directly from IGeometryBuilder it can use the protected member functions of the IGeometryBuilder
      to glue Volumes together and exchange BoundarySurfaces, as the IGeometryBuilder has friend privileges to Trk::Volume
      and Trk::BoundarySurface.
      
      @author Andreas.Salzburger@cern.ch
       adapted from InDet for Muons by Sarka.Todorova@cern.ch
    */
    
  class MuonTrackingGeometryBuilder : public AlgTool,
                               virtual public Trk::IGeometryBuilder {
  public:
      /** Constructor */
      MuonTrackingGeometryBuilder(const std::string&,const std::string&,const IInterface*);
      /** Destructor */
      virtual ~MuonTrackingGeometryBuilder();
      /** AlgTool initailize method.*/
      StatusCode initialize();
      /** AlgTool finalize method */
      StatusCode finalize();
      /** TrackingGeometry Interface methode */
      const Trk::TrackingGeometry* trackingGeometry(const Trk::TrackingVolume* tvol = 0) const; 

    private:
      /** Private method to fill default material */
      //void fillDefaultServiceMaterial();

           
      Trk::IMagneticFieldTool*            m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      std::string                         m_magFieldToolName;            //!< Name of the Tracking Magnetic Field Svc
      std::string                         m_magFieldToolInstanceName;    //!< Instance Name of Tracking Magnetic Field Svc
      
      /*
      Trk::ILayerArrayCreator*            m_layerArrayCreator;            //!< A Tool for coherent LayerArray creation
      std::string                         m_layerArrayCreatorName;        //!< Name of the LayerCreator implementation
      std::string                         m_layerArrayCreatorInstanceName;//!< Instance Name of the Layer Creator
      */
      
      Trk::ITrackingVolumeArrayCreator*   m_trackingVolumeArrayCreator;             //!< Helper Tool to create TrackingVolume Arrays
      std::string                         m_trackingVolumeArrayCreatorName;         //!< Name of the helper tool
      std::string                         m_trackingVolumeArrayCreatorInstanceName; //!< Instance of the helper tool

      Trk::ITrackingVolumeHelper*        m_trackingVolumeHelper;             //!< Helper Tool to create TrackingVolumes
      std::string                        m_trackingVolumeHelperName;         //!< Name of the helper tool
      std::string                        m_trackingVolumeHelperInstanceName; //!< Instance of the helper tool

      Trk::ITrackingVolumeDisplayer*     m_trackingVolumeDisplayer;             //!< Displayer Tool to create TrackingVolumes
      std::string                        m_trackingVolumeDisplayerName;         //!< Name of the helper tool
      std::string                        m_trackingVolumeDisplayerInstanceName; //!< Instance of the helper tool
   
      bool                                m_muonSimple;
      //bool                                m_muonStandalone;

      bool                                m_enclosingVolumes;             //!< build the enclosing volumes for glueing with outer detectors

      // Overall Dimensions
      mutable double                      m_innerBarrelRadius;             //!< minimal extend in radial dimension of the muon barrel
      double                              m_outerBarrelRadius;             //!< maximal extend in radial dimension of the muon barrel
      mutable double                      m_barrelZ;                  //!< maximal extend in z of the muon barrel
      double                              m_innerEndcapZ;             //!< maximal extend in z of the inner part of muon endcap 
      double                              m_outerEndcapZ;             //!< maximal extend in z of the outer part of muon endcap
      double                              m_beamPipeRadius;
      double                              m_outerEndcapRadius;

      mutable Trk::MaterialProperties                  m_muonMaterial;               //!< the material
      Trk::MagneticFieldProperties        m_muonMagneticField;          //!< the magnetic Field

      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      mutable Trk::TrackingVolume*    m_standaloneTrackingVolume;   // muon standalone tracking volume                 
 };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H


