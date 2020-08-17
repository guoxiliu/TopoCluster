ttk_add_paraview_plugin(ttkFG_Segmentation
	SOURCES ${VTKWRAPPER_DIR}/ttkFG_Segmentation/ttkFG_Segmentation.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML FG_Segmentation.xml
	LINK fG_Segmentation)

