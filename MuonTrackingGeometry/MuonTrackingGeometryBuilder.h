///////////////////////////////////////////////////////////////////
// MuonTrackingGeometryBuilder.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H
#define MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H

//Trk
#include "TrkDetDescrInterfaces/IGeometryBuilder.h"
#include "TrkDetDescrInterfaces/ITrackingVolumesSvc.h"
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkGeometry/MaterialProperties.h"
// Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"
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
      /** Private methods to define subvolumes and fill them with detached volumes */
      const Trk::TrackingVolume* processVolume( const Trk::Volume*, int, int, std::string) const; 
      const Trk::TrackingVolume* processVolume( const Trk::Volume*, int, std::string) const; 
      /** Private method to check volume properties */
      void checkVolume(const Trk::TrackingVolume*) const;
      /** Private method to find detached volumes */
      std::vector<const Trk::DetachedTrackingVolume*>* getDetachedObjects(const Trk::Volume*) const;
      /** Private method to retrieve z/phi/h partition */
      void getZParts() const;
      void getPhiParts() const;
      void getHParts() const;
      double calculateVolume(const Trk::Volume*) const;
           
      ToolHandle<Trk::IMagneticFieldTool>                  m_magFieldTool;                  //!< Tracking Interface to Magnetic Field

      ToolHandle<Trk::IDetachedTrackingVolumeBuilder>      m_stationBuilder;                //!< A Tool for station type creation

      ToolHandle<Trk::IDetachedTrackingVolumeBuilder>      m_inertBuilder;                  //!< A Tool for inert object  creation
      
      ToolHandle<Trk::ITrackingVolumeArrayCreator>         m_trackingVolumeArrayCreator;    //!< Helper Tool to create TrackingVolume Arrays

      ToolHandle<Trk::ITrackingVolumeHelper>               m_trackingVolumeHelper;          //!< Helper Tool to create TrackingVolumes

      ToolHandle<Trk::ITrackingVolumeDisplayer>            m_trackingVolumeDisplayer;       //!< Displayer Tool to create TrackingVolumes

      ServiceHandle<ITrackingVolumesSvc>                   m_tvSvc;       //!< service to provide input volume size

   
      bool                                m_muonSimple;
      bool                                m_loadMSentry;
      bool                                m_muonActive;
      bool                                m_muonInert;

      // Overall Dimensions
      mutable double                      m_innerBarrelRadius;             //!< minimal extend in radial dimension of the muon barrel
      mutable double                      m_outerBarrelRadius;             //!< maximal extend in radial dimension of the muon barrel
      mutable double                      m_barrelZ;                  //!< maximal extend in z of the muon barrel
      double                              m_innerEndcapZ;             //!< maximal extend in z of the inner part of muon endcap 
      double                              m_outerEndcapZ;             //!< maximal extend in z of the outer part of muon endcap
      double                              m_beamPipeRadius;
      double                              m_innerShieldRadius;
      double                              m_outerShieldRadius;
      double                              m_diskShieldZ;

      mutable Trk::MaterialProperties     m_muonMaterial;               //!< the (empty) material
      mutable Trk::MagneticFieldProperties  m_muonMagneticField;          //!< the magnetic Field

      mutable std::vector< double >       m_muonMaterialProperties;     //!< The material properties of the created muon system 
      mutable Trk::TrackingVolume*        m_standaloneTrackingVolume;   // muon standalone tracking volume                 
      int                                 m_barrelEtaPartition;
      int                                 m_innerEndcapEtaPartition;
      int                                 m_outerEndcapEtaPartition;
      int                                 m_phiPartition;
      bool                                m_adjustStatic;
      bool                                m_static3d;
      bool                                m_blendInertMaterial; 
      mutable double                      m_alignTolerance;

      mutable const std::vector<const Trk::DetachedTrackingVolume*>*    m_stations;    // muon chambers 
      mutable const std::vector<const Trk::DetachedTrackingVolume*>*    m_inertObjs;   // muon inert material 
      mutable const std::vector<const Span*>*                     m_stationSpan; 
      mutable const std::vector<const Span*>*                     m_inertSpan; 
      mutable std::vector<double>                                 m_zPartitions;
      mutable std::vector<int>                                    m_zPartitionsType;
      mutable std::vector<double>                                 m_adjustedPhi;
      mutable std::vector<int>                                    m_adjustedPhiType;
      mutable std::vector<std::vector<std::vector<std::vector<std::pair<int,double> > > > > m_hPartitions;
 };


} // end of namespace


#endif //MUONTRACKINGGEOMETRY_MUONTRACKINGGEOMETRYBUILDER_H


