cmake_minimum_required(VERSION 3.2)

project(testTopoClusterCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(testTopoClusterCmd main.cpp)
target_link_libraries(testTopoClusterCmd ${VTK_LIBRARIES}
	ttkTestTopoCluster ttkProgramBase)

install(TARGETS testTopoClusterCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
