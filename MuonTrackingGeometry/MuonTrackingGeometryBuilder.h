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
#include "MuonTrackingGeometry/MuonStationBuilder.h"
#include "MuonTrackingGeometry/MuonInertMaterialBuilder.h"

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

  typedef std::vector<double> Span;
     
  /** @class MuonTrackingGeometryBuilder
  
      The Muon::MuonTrackingGeometryBuilder retrieves MuonStationBuilder and MuonInertMaterialBuilder
      for the Muon Detector sub detectors and combines the given Volumes to a full Trk::TrackingGeometry.
      
      Inheriting directly from IGeometryBuilder it can use the protected member functions of the IGeometryBuilder
      to glue Volumes together and exchange BoundarySurfaces
      
      @author Sarka.Todorova@cern.ch
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
      /** Private method to find z/phi span of detached volumes */
      const Span* findVolumeSpan(const Trk::VolumeBounds* volBounds, HepTransform3D transf, double zTol, double phiTol) const;
      const std::vector<const Span*>* findVolumesSpan(const std::vector<const Trk::DetachedTrackingVolume*>*& objs, double zTol, double phiTol) const;
      /** Private method to define subvolumes and fill them with detached volumes */
      const Trk::TrackingVolume* processVolume( const Trk::Volume*, int, int, std::string) const; 
      /** Private method to check volume properties */
      void checkVolume(const Trk::TrackingVolume*) const;
      /** Private method to find detached volumes */
      std::vector<const Trk::DetachedTrackingVolume*>* getDetachedObjects(const Trk::Volume*) const;
           
      Trk::IMagneticFieldTool*            m_magFieldTool;                //!< Tracking Interface to Magnetic Field
      std::string                         m_magFieldToolName;            //!< Name of the Tracking Magnetic Field Svc
      std::string                         m_magFieldToolInstanceName;    //!< Instance Name of Tracking Magnetic Field Svc

      Trk::IDetachedTrackingVolumeBuilder*      m_stationBuilder;            //!< A Tool for station type creation
      std::string                         m_stationBuilderName;        //!< Name of the station type implementation
      std::string                         m_stationBuilderInstanceName;//!< Instance Name of the station type builder
      Trk::IDetachedTrackingVolumeBuilder*      m_inertBuilder;            //!< A Tool for inert object  creation
      std::string                         m_inertBuilderName;        //!< Name of the inert object  implementation
      std::string                         m_inertBuilderInstanceName;//!< Instance Name of the inert object type builder
      
      Trk::ITrackingVolumeArrayCreator*   m_trackingVolumeArrayCreator;             //!< Helper Tool to create TrackingVolume Arrays
      std::string                         m_trackingVolumeArrayCreatorName;         //!< Name of the helper tool
      std::string                         m_trackingVolumeArrayCreatorInstanceName; //!< Instance of the helper tool

      Trk::ITrackingVolumeHelper*         m_trackingVolumeHelper;             //!< Helper Tool to create TrackingVolumes
      std::string                         m_trackingVolumeHelperName;         //!< Name of the helper tool
      std::string                         m_trackingVolumeHelperInstanceName; //!< Instance of the helper tool

      Trk::ITrackingVolumeDisplayer*      m_trackingVolumeDisplayer;             //!< Displayer Tool to create TrackingVolumes
      std::string                         m_trackingVolumeDisplayerName;         //!< Name of the helper tool
      std::string                         m_trackingVolumeDisplayerInstanceName; //!< Instance of the helper tool
   
      bool                                m_muonSimple;
      bool                                m_muonActive;
      bool                                m_muonInert;

      // Overall Dimensions
      mutable double                      m_innerBarrelRadius;             //!< minimal extend in radial dimension of the muon barrel
      double                              m_outerBarrelRadius;             //!< maximal extend in radial dimension of the muon barrel
      mutable double                      m_barrelZ;                  //!< maximal extend in z of the muon barrel
      double                              m_innerEndcapZ;             //!< maximal extend in z of the inner part of muon endcap 
      double                              m_outerEndcapZ;             //!< maximal extend in z of the outer part of muon endcap
      double                              m_beamPipeRadius;
      double                              m_outerEndcapRadius;

      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the (empty) material
      Trk::MagneticFieldProperties        m_muonMagneticField;          //!< the magnetic Field

      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      mutable Trk::TrackingVolume*        m_standaloneTrackingVolume;   // muon standalone tracking volume                 
      int                                 m_barrelEtaPartition;
      int                                 m_innerEndcapEtaPartition;
      int                                 m_outerEndcapEtaPartition;
      int                                 m_phiPartition;
      mutable const std::vector<const Trk::DetachedTrackingVolume*>*    m_stations;    // muon chambers 
      mutable const std::vector<const Trk::DetachedTrackingVolume*>*    m_inertObjs;   // muon inert material 
      mutable const std::vector<const Span*>*                     m_stationSpan; 
      mutable const std::vector<const Span*>*                     m_inertSpan; 
 };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H


