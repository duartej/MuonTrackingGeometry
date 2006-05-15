#include "GaudiKernel/DeclareFactoryEntries.h"
#include "MuonTrackingGeometry/MuonTrackingGeometryBuilder.h"
#include "MuonTrackingGeometry/MuonStationBuilder.h"
#include "MuonTrackingGeometry/MuonStationTypeBuilder.h"


using namespace Muon;

DECLARE_TOOL_FACTORY( MuonTrackingGeometryBuilder );
DECLARE_TOOL_FACTORY( MuonStationBuilder );
DECLARE_TOOL_FACTORY( MuonStationTypeBuilder );

DECLARE_FACTORY_ENTRIES( MuonTrackingGeometry )
{
  DECLARE_TOOL( MuonTrackingGeometryBuilder );  
  DECLARE_TOOL( MuonStationBuilder );  
  DECLARE_TOOL( MuonStationTypeBuilder );  
}


