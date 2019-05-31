ttk_add_base_library(multivariateGradient
	SOURCES MultivariateGradient.cpp outputCells.cpp
	HEADERS MultivariateGradient.h CriticalClusterGraph.h multivariategradient_template.h
	LINK triangulation)

