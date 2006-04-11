#include "GaudiKernel/DeclareFactoryEntries.h"
#include "MuonTrackingGeometry/MuonTrackingGeometryBuilder.h"


using namespace Muon;

DECLARE_TOOL_FACTORY( MuonTrackingGeometryBuilder )

DECLARE_FACTORY_ENTRIES( MuonTrackingGeometry )
{
    DECLARE_TOOL( MuonTrackingGeometryBuilder )  
};


