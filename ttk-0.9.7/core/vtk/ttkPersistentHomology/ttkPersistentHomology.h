/// \ingroup vtk
/// \class ttkPersistentHomology
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the persistentHomology processing package.
///
/// VTK wrapping code for the @PersistentHomology package.
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
/// \sa ttk::PersistentHomology
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
#include                  <PersistentHomology.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistentHomology
#else
class ttkPersistentHomology
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkPersistentHomology* New();
    vtkTypeMacro(ttkPersistentHomology, vtkDataSetAlgorithm)
    
    
    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);
    

  
    int FillInputPortInformation(int port, vtkInformation *info) override {
      
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
          break;
        default:
          break;
      }
      
      return 1;
    }
    
    int FillOutputPortInformation(int port, vtkInformation *info) override {
      
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
          break;
        default:
          break;
      }
      
      return 1;
    }
    
    
  protected:
   
    ttkPersistentHomology(){
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }
    
    ~ttkPersistentHomology(){};
    
    TTK_SETUP();
    
    
  private:
    
    std::string           ScalarField;
    ttk::PersistentHomology            persistentHomology_;
    
};
