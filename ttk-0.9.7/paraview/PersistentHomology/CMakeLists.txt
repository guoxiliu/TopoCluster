ttk_add_paraview_plugin(ttkPersistentHomology
	SOURCES ${VTKWRAPPER_DIR}/ttkPersistentHomology/ttkPersistentHomology.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PersistentHomology.xml
	LINK persistentHomology)

