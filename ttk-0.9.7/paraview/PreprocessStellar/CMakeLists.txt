ttk_add_paraview_plugin(ttkPreprocessStellar
	SOURCES ${VTKWRAPPER_DIR}/ttkPreprocessStellar/ttkPreprocessStellar.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PreprocessStellar.xml
	LINK preprocessStellar)

