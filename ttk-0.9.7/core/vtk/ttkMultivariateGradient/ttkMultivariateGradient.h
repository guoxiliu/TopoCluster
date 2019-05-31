/// \ingroup vtk
/// \class ttkMultivariateGradient
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the multivariateGradient processing package.
///
/// VTK wrapping code for the @MultivariateGradient package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MultivariateGradient
#pragma once

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <MultivariateGradient.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkMultivariateGradient
#else
class ttkMultivariateGradient
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkMultivariateGradient* New();
    vtkTypeMacro(ttkMultivariateGradient, vtkDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    vtkSetMacro(extractInfluence, bool);
    vtkSetMacro(extractAscending, bool);
    vtkSetMacro(extractDescending, bool);
    vtkSetMacro(extractGradient, bool);
    vtkSetMacro(clusterInfluence, int);

    void SetSelectedFields(std::string name){
        if(multivariateGradient_ != NULL){
          multivariateGradient_ = NULL;
          SelectedFields.clear();
        }
        SelectedFields.push_back(name);
        Modified();
    }

    void SetMaxPathLenght(int path){
        Path = path;
        Modified();
    }

    void SetMaxClusterSize(int cl_size){
      CL_size = cl_size;
      Modified();
    }


    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters


    // TODO-4
    // set-getters macros to define from each variable you want to access from
    // the outside (in particular from paraview) - to adapt.
    // Note that the XML file for the ParaView plug-in specification needs to be
    // edited accordingly.

//    vtkSetMacro(ScalarField, std::string);
//    vtkGetMacro(ScalarField, std::string);
    // end of TODO-4

    // TODO-2
    // Over-ride the input types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the input types, on a per input basis.
    // In this example, the first input type is forced to vtkUnstructuredGrid.
    // The second input type is forced to vtkImageData.
     int FillInputPortInformation(int port, vtkInformation *info) override
     {
         switch(port){
           case 0:
             info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
             break;
         }

         return 1;
     }
    // end of TODO-2

    // TODO-3
    // Over-ride the output types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the output types, on a per output basis.
    // In this example, the first output type is forced to vtkUnstructuredGrid.
    // The second output type is forced to vtkImageData.
     int FillOutputPortInformation(int port, vtkInformation *info) override {
         switch(port){
           case 0: //critical cells
           case 1: //gradient pairs
           case 2: //critical graph
           case 3: //cluster involved
           case 4: //region of influence
           case 5: //region of influence
              info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
             break;
         }

         return 1;
     }
    // end of TODO-3

    void outputGradient(std::vector<vtkDataSet *> &outputs);
    void outputCriticalCells(std::vector<vtkDataSet *> &outputs);
    void outputClustersGraph(std::vector<vtkDataSet *> &outputs);
    void outputOneCluster(std::vector<vtkDataSet *> &outputs);
    void outputAscendingRegionOfInfluence(std::vector<vtkDataSet *> &outputs);
    void outputDescendingRegionOfInfluence(std::vector<vtkDataSet *> &outputs);


  protected:

    ttkMultivariateGradient(){


      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;

      extractInfluence = false;
      extractAscending = false;
      extractDescending= false;
      clusterInfluence = 0;

      Path = 0;
      CL_size = 0;

      multivariateGradient_=NULL;

      // TODO-1
      // Specify the number of input and output ports.
      // By default, this filter has one input and one output.
      // In this example, we define 2 inputs and 2 outputs.
       SetNumberOfInputPorts(1);
       SetNumberOfOutputPorts(6);
      // end of TODO-1
    }

    ~ttkMultivariateGradient(){};

    TTK_SETUP();


  private:

    std::vector<std::string>   SelectedFields;
    int Path;
    int CL_size;

    bool extractInfluence;
    bool extractAscending;
    bool extractDescending;
    bool extractGradient;
    int clusterInfluence;

    ttk::MultivariateGradient*            multivariateGradient_;

};
