/// \ingroup vtk
/// \class ttkTestTopoCluster
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Jan. 2020.
///
/// \brief TTK VTK-filter that wraps the testTopoCluster processing package.
///
/// VTK wrapping code for the @TestTopoCluster package.
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
/// \sa ttk::TestTopoCluster
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
#include                  <TestTopoCluster.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTestTopoCluster
#else
class ttkTestTopoCluster
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkTestTopoCluster* New();
    vtkTypeMacro(ttkTestTopoCluster, vtkDataSetAlgorithm)
    
    // default ttk setters
    vtkSetMacro(debugLevel_, int);
    
    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
    
    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);

    vtkSetMacro(CacheRatio, float);
    vtkGetMacro(CacheRatio, float);
    
    // TODO-2
    // Over-ride the input types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the input types, on a per input basis.
    // In this example, the first input type is forced to vtkUnstructuredGrid.
    // The second input type is forced to vtkImageData.
//     int FillInputPortInformation(int port, vtkInformation *info) override {
//       
//       switch(port){
//         case 0:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
//           break;
//         case 1:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData"); 
//           break;
//         default:
//           break;
//       }
//       
//       return 1;
//     }
    // end of TODO-2
    
    // TODO-3
    // Over-ride the output types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the output types, on a per output basis.
    // In this example, the first output type is forced to vtkUnstructuredGrid.
    // The second output type is forced to vtkImageData.
//     int FillOutputPortInformation(int port, vtkInformation *info) override {
//       
//       switch(port){
//         case 0:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
//           break;
//         case 1:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData"); 
//           break;
//         default:
//           break;
//       }
//       
//       return 1;
//     }
    // end of TODO-3
    
    
  protected:
   
    ttkTestTopoCluster(){
      
      // init
      outputScalarField_ = NULL;
      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;
      CacheRatio = 0.1;
    }
    
    ~ttkTestTopoCluster(){};
    
    TTK_SETUP();
    
    
  private:
    float                   CacheRatio;
    std::string           ScalarField;
    vtkDataArray          *outputScalarField_;
    ttk::TestTopoCluster      testTopoCluster_;
};
