set (TOOLS_SRCS    
  CubicSpline.cc   
  MuffinTin.cc     
  ParseCommand.cc  
  Respline.cc      
  ExportCasino.cc  
  ParserClass.cc   
  ParseCommand.cc  
  AtomicOrbital.cc 
  APWClass.cc      
  CellClass.cc     
  CenterClass.cc   
  CubicSpline.cc   
  Gaussian.cc      
  LocalizeOrbs.cc  
  MuffinTin.cc     
  OrbitalClass.cc  
  OrbitalSetClass.cc
  ParserClass.cc   
  ParseCommand.cc  
  ReadABINIT.cc    
  Read_ESHDF.cc    
  Read_FPMD.cc     
  ReadLAPW.cc      
  TileMapClass.cc  
  WriteCASINO.cc   
  Write_qmcPACK.cc 
  Write_ESHDF.cc   
  XMLWriterClass2.cc
  )
#  LAPWClass.cc     
#  PlaneWaveClass.cc

ADD_LIBRARY(qmctools ${TOOLS_SRCS})
link_libraries(qmctools common einspline)

ADD_EXECUTABLE(wfconvert  WFConv.cc)        
ADD_EXECUTABLE(ppconvert  NLPPClass.cc)        
#ADD_EXECUTABLE(fit_sk  FitSk.cc)
#ADD_EXECUTABLE(test_tile  TestTile.cc)
#ADD_EXECUTABLE(test_spline  TestSpline.cc)
#
