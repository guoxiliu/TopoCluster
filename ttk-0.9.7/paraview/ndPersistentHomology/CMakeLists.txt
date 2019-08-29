ttk_add_paraview_plugin(ttkndPersistentHomology
	SOURCES ${VTKWRAPPER_DIR}/ttkndPersistentHomology/ttkndPersistentHomology.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ndPersistentHomology.xml
	LINK ndPersistentHomology)

