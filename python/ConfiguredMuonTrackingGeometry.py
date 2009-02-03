
from AthenaCommon.AppMgr import ToolSvc

# muon active/passive geometry setup
from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonStationBuilder
MuonStationBuilder= Muon__MuonStationBuilder(name = 'MuonStationBuilder')
ToolSvc += MuonStationBuilder

from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonInertMaterialBuilder
MuonInertMaterialBuilder= Muon__MuonInertMaterialBuilder(name = 'MuonInertMaterialBuilder')
ToolSvc += MuonInertMaterialBuilder 

from TrkDetDescrSvc.TrkDetDescrJobProperties import TrkDetFlags
# muon tracking geometry builder
from MuonTrackingGeometry.MuonTrackingGeometryConf import Muon__MuonTrackingGeometryBuilder
MuonTrackingGeometryBuilder= Muon__MuonTrackingGeometryBuilder(name = 'MuonTrackingGeometryBuilder')
MuonTrackingGeometryBuilder.EntryVolumeName = TrkDetFlags.MuonSystemEntryVolumeName.get_Value()
MuonTrackingGeometryBuilder.ExitVolumeName  = TrkDetFlags.MuonSystemContainerName.get_Value()
ToolSvc += MuonTrackingGeometryBuilder

print MuonTrackingGeometryBuilder
