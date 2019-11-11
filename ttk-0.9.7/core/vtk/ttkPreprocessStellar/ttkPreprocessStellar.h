/// \ingroup vtk
/// \class ttkPreprocessStellar
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the preprocessStellar processing package.
///
/// VTK wrapping code for the @PreprocessStellar package.
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
/// \sa ttk::PreprocessStellar
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
#include                  <PreprocessStellar.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPreprocessStellar
#else
class ttkPreprocessStellar
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkPreprocessStellar* New();
    vtkTypeMacro(ttkPreprocessStellar, vtkDataSetAlgorithm)
    
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
    
        
    // set-getters macros to define from each variable you want to access from 
    // the outside (in particular from paraview) - to adapt.
    vtkSetMacro(Threshold, int);
    vtkGetMacro(Threshold, int);
   
    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);

    // Over-ride the input types.
    // int FillInputPortInformation(int port, vtkInformation *info) override {
      
    //   switch(port){
    //     case 0:
    //       info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
    //       break;
    //     case 1:
    //       info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData"); 
    //       break;
    //     default:
    //       break;
    //   }
      
    //   return 1;
    // }
    
    // Over-ride the output types.
    // int FillOutputPortInformation(int port, vtkInformation *info) override {
      
    //   switch(port){
    //     case 0:
    //       info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
    //       break;
    //     case 1:
    //       info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData"); 
    //       break;
    //     default:
    //       break;
    //   }
      
    //   return 1;
    // }
    
    
  protected:
   
    ttkPreprocessStellar(){
      
        // init
      Threshold = 1000;
      outputScalarField_ = NULL;
      
      UseAllCores = true;
      ThreadNumber = 1;
      debugLevel_ = 3;
      
    }
    
    ~ttkPreprocessStellar(){};
    
    TTK_SETUP();
    
    
  private:
    
    int                   Threshold;
    std::string           ScalarField;
    vtkDataArray          *outputScalarField_;
    vector<SimplexId>     *vertexArray, *nodeArray, *cellArray;
    ttk::PreprocessStellar            preprocessStellar_;
    
};
