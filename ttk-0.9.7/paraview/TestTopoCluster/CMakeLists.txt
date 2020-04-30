ttk_add_paraview_plugin(ttkTestTopoCluster
	SOURCES ${VTKWRAPPER_DIR}/ttkTestTopoCluster/ttkTestTopoCluster.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TestTopoCluster.xml
	LINK testTopoCluster)

