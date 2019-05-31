ttk_add_paraview_plugin(ttkMultivariateGradient
	SOURCES ${VTKWRAPPER_DIR}/ttkMultivariateGradient/ttkMultivariateGradient.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML MultivariateGradient.xml
	LINK multivariateGradient)

