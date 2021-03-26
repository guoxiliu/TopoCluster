# TopoCluster

TopoCluster is a new localized data structure for tetrahedral meshes, which provides efficient computation of theconnectivity of the mesh elements with a low memory footprint. 


## Installation

TopoCluster is developed under Topology Toolkit (TTK) framework. To use TopoCluster, TTK and ParaView need to be installed.

The following installation steps are based on [TTK offical installation guide](https://topology-tool-kit.github.io/installation-0.9.7.html) under Linux platform (based on Ubuntu Linux distribution).


### 1. Downloads
Use `git clone` to clone the repository.

### 2. Installing the dependencies
Please enter the following commands (omit the `$` character) in a terminal to install dependencies.

```
$ sudo apt-get install cmake-qt-gui
$ sudo apt-get install libvtk7-dev
$ sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev
```
### 3. Configuring, building and installing ParaView
Use `cd` command to change the current working directory to the project root when executing the following commands! 

#### Configuration 
To enter the configuration menu of ParaView's build, enter the following commands:

```
$ cd ./ParaView-v5.6.0/
$ mkdir build
$ cd build
$ cmake-gui ../
```

Click on the "Configure" button to proceed. Once the configuration is finished, please tick the "Advanced" check box and set the following variables as follows (required for TTK's installation):
- CMAKE_BUILD_TYPE=Release
- PARAVIEW_ENABLE_PYTHON=ON
- PARAVIEW_INSTALL_DEVELOPMENT_FILES=ON
- VTK_PYTHON_VERSION=3

Next, click on the "Generate" button and close the configuration window when the generation is completed.

#### Build 
Use `make -jN` command to start the compilation process, when `N` is the number of available cores on your system (this will take a **LONG** time).

#### Installation
Use `sudo make install` to install the build of ParaView on your system. 

### 4. Configuring, building and installing TTK
Use `cd` command to change the current working directory to the project root when executing the following commands! 

#### Configuration
To enter the configuration menu of ParaView's build, enter the following commands:

```
$ cd ./ttk-0.9.7/
$ mkdir build
$ cd build
$ cmake-gui ../
```

Click on the "Configure" button to proceed. Then, click on the "Generate" button. Once the generation is completed, close the configuration window.

Note: You can choose to build Explicit TopoCluster by setting `ENABLE_IMPLICIT_TOPOCLUSTER` to `OFF` (default) and to build Implicit TopoCluster by setting `ENABLE_IMPLICIT_TOPOCLUSTER` to `ON`.

#### Build 
Use `make -jN` command to start the compilation process, when `N` is the number of available cores on your system (this will take a **LONG** time).

#### Installation
Once the build is finished, use `sudo make install` to install your build of TTK on your system. After this, you have successfully installed TopoCluster.

