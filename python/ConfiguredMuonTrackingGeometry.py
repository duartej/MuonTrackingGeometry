
from AthenaCommon.AppMgr import ToolSvc

from TrkDetDescrSvc.TrkDetDescrJobProperties import TrkDetFlags

from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonStationTypeBuilder
MuonStationTypeBuilder= Muon__MuonStationTypeBuilder(name = 'MuonStationTypeBuilder')
MuonStationTypeBuilder.MagneticFieldMode = TrkDetFlags.MagneticFieldMode()
ToolSvc += MuonStationTypeBuilder

# muon active/passive geometry setup
from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonStationBuilder
MuonStationBuilder= Muon__MuonStationBuilder(name = 'MuonStationBuilder')
MuonStationBuilder.StationTypeBuilder = MuonStationTypeBuilder
MuonStationBuilder.MagneticFieldMode  = TrkDetFlags.MagneticFieldMode()
ToolSvc += MuonStationBuilder

from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonInertMaterialBuilder
MuonInertMaterialBuilder= Muon__MuonInertMaterialBuilder(name = 'MuonInertMaterialBuilder')
ToolSvc += MuonInertMaterialBuilder 

# muon tracking geometry builder
from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonTrackingGeometryBuilder
MuonTrackingGeometryBuilder= Muon__MuonTrackingGeometryBuilder(name = 'MuonTrackingGeometryBuilder')
MuonTrackingGeometryBuilder.EntryVolumeName = TrkDetFlags.MuonSystemEntryVolumeName()
MuonTrackingGeometryBuilder.ExitVolumeName  = TrkDetFlags.MuonSystemContainerName()
MuonTrackingGeometryBuilder.MagneticFieldMode = TrkDetFlags.MagneticFieldMode()
ToolSvc += MuonTrackingGeometryBuilder

#print MuonTrackingGeometryBuilder
