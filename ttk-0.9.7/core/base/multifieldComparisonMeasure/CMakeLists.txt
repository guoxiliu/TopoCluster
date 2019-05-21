find_package (Eigen3 3.3 REQUIRED NO_MODULE)

ttk_add_base_library(multifieldComparisonMeasure
	SOURCES MultifieldComparisonMeasure.cpp
	HEADERS MultifieldComparisonMeasure.h
	LINK triangulation Eigen3::Eigen)

