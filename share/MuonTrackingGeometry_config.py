
# muon active/passive geometry setup
from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonStationBuilder
MuonStationBuilder= Muon__MuonStationBuilder(name = 'MuonStationBuilder',
                                             BuildBarrelStations = True,     
                                             BuildEndcapStations = True,
                                             BuildCSCStations = True,
                                             BuildTGCStations = True,
                                             IdentifyActiveLayers = True )
ToolSvc += MuonStationBuilder

from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonInertMaterialBuilder
MuonInertMaterialBuilder= Muon__MuonInertMaterialBuilder(name = 'MuonInertMaterialBuilder',
                                                         SimplifyGeometry = True,
                                                         SimplifyGeometryToLayers = False )
ToolSvc += MuonInertMaterialBuilder 

# muon tracking geometry builder
from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonTrackingGeometryBuilder
MuonTrackingGeometryBuilder= Muon__MuonTrackingGeometryBuilder(name = 'MuonTrackingGeometryBuilder')

ToolSvc += MuonTrackingGeometryBuilder
ToolSvc.MuonTrackingGeometryBuilder.BuildActiveMaterial = True
ToolSvc.MuonTrackingGeometryBuilder.BuildInertMaterial = True
ToolSvc.MuonTrackingGeometryBuilder.AdjustStatic = True
ToolSvc.MuonTrackingGeometryBuilder.StaticPartition3D = True

print MuonTrackingGeometryBuilder
