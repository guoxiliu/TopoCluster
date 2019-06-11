ttk_add_base_library(persistentHomology
	SOURCES PersistentHomology.cpp
	HEADERS PersistentHomology.h persistenthomology_template.h
	LINK triangulation boundarymatrix)

