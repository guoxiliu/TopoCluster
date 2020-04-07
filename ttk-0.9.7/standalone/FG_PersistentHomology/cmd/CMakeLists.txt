cmake_minimum_required(VERSION 3.2)

project(fG_PersistentHomologyCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(fG_PersistentHomologyCmd main.cpp)
target_link_libraries(fG_PersistentHomologyCmd ${VTK_LIBRARIES}
	ttkFG_PersistentHomology ttkProgramBase)

install(TARGETS fG_PersistentHomologyCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
