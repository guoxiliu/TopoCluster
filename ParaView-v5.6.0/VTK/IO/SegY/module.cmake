vtk_module(vtkIOSegY
  GROUPS
    StandAlone
  TEST_DEPENDS
    vtkInteractionStyle
    vtkRenderingOpenGL2
    vtkTestingCore
    vtkTestingRendering
  KIT
    vtkIO
  DEPENDS
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkIOImage
  PRIVATE_DEPENDS
    vtkCommonCore
  )