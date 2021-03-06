vtk_module(VisItLib
  EXCLUDE_FROM_WRAPPING
  DEPENDS
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkCommonMisc
    vtkFiltersFlowPaths
    vtkIOImage
    vtkIOLegacy
  PRIVATE_DEPENDS
    vtkCommonMath
    vtkCommonTransforms
    vtkFiltersCore
    vtkFiltersExtraction
    vtkFiltersGeneral
    vtkFiltersGeometry
    vtkFiltersSources
    vtkImagingHybrid
    vtkRenderingCore
    vtkRenderingVolume${VTK_RENDERING_BACKEND}
    vtkpng
    vtkzlib
    # old modules seem not used. Will remove them in the future.
    vtkFiltersAMR
    vtkIOGeometry
    vtkParallelCore
  )
