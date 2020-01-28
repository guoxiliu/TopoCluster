ttk_add_paraview_plugin(ttkTestStellar
	SOURCES ${VTKWRAPPER_DIR}/ttkTestStellar/ttkTestStellar.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TestStellar.xml
	LINK testStellar)

