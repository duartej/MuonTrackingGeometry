from AthenaCommon.GlobalFlags import GlobalFlags
###################################################################################
# ConfiguredMuonGeometryBuilder
###################################################################################
# Python Setup Class for Muon::MuonGeometryBuilder
#                    and Muon::MuonCTB_GeometryBuilder
#
# Author: Andreas.Salzburger@cern.ch
#
# Date: 16/07/2005
#
# adapted by Sarka.Todorova@cern.ch
###################################################################################
#  TODO: write documentation
###################################################################################

class ConfiguredMuonGeometryBuilder :
   def __init__(self, instname         = None,  
                      magfieldtool     = None,
                      args = 0):
                      
        #load the Dll for this
        if 'MuonTrackingGeometry' not in theApp.Dlls:
           theApp.Dlls += [ 'MuonTrackingGeometry' ]
        
        #load the Dll for the TrkDetDescrTools
        if 'TrkDetDescrTools' not in theApp.Dlls:
           theApp.Dlls += [ 'TrkDetDescrTools' ]
           
        # set the service name and reset it if CTB 2004 setup
        self.__svcname__      = 'Muon::MuonTrackingGeometryBuilder'
        if GlobalFlags.DetGeo.is_ctb() :
          self.__svcname__       = 'Muon::MuonCTB_TrackingGeometryBuilder'
        
        # set (a default) instance name for the GeometryBuilder
        if instname is None:
          instname = 'MuonTrackingGeometryBuilder'
        self.__instname__     = instname
        
        # create the link to the ToolSvc name
        self.__toolname__ = 'ToolSvc.'+self.__instname__
        # self.__thisGeoBuilder__ = Service( self.__toolname__ )
        self.__thisGeoBuilder = Service( self.__toolname__ )
                
        # set the magnetic field tool with appropriate method
        self.setMagneticFieldTool(magfieldtool)
        # set the beampipe builder with the appropriate method
        #if not GlobalFlags.DetGeo.is_ctb() :
        #   self.setBeamPipeBuilder(bpbuilder)
        # set the pixel volume builder with the appropriate method
        #self.setPixelVolumeBuilder(pixvolumebuilder, magfieldtool)
        # set the SCT volume builder with the appropriate method
        #self.setSCT_VolumeBuilder(sctvolumebuilder, magfieldtool)
        # set the TRT volume builder with the appropriate method
        #self.setTRT_VolumeBuilder(trtvolumebuilder, magfieldtool)
        
           
   # Set the MsgStream level      
   def msgStreamLevel(self, level):
        self.__thisGeoBuilder.OutputLevel = level
        
   def setMagneticFieldTool(self, magfieldtool=None) :
        # set (a default) magnetic field tool to be used by the GeometryBuilder
        if magfieldtool is None:
           include ( 'TrkMagFieldTools/ConfiguredMagneticFieldTool.py' )
           magfieldtool = ConfiguredMagneticFieldTool()
        # declare the magnetic field tool
        self.__thisGeoBuilder.MagneticFieldTool           = magfieldtool.name()
        self.__thisGeoBuilder.MagneticFieldToolInstance   = magfieldtool.instance()       

          
   # Return method for service name
   def name(self):
       return self.__svcname__
   
   # Return method for instance name
   def instance(self):
       return self.__instname__        
               
   # output
   def printInfo(self):
       print '***** ConfiguredGeometryBuilder *******************************************************'
       print '*'
          
               
