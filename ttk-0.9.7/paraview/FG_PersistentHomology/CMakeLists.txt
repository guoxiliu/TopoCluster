ttk_add_paraview_plugin(ttkFG_PersistentHomology
	SOURCES ${VTKWRAPPER_DIR}/ttkFG_PersistentHomology/ttkFG_PersistentHomology.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML FG_PersistentHomology.xml
	LINK fG_PersistentHomology)

