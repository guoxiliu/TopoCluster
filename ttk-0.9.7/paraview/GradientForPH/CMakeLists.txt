ttk_add_paraview_plugin(ttkGradientForPH
	SOURCES ${VTKWRAPPER_DIR}/ttkGradientForPH/ttkGradientForPH.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML GradientForPH.xml
	LINK gradientForPH)

