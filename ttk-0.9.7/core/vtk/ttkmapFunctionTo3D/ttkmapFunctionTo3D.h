/// \ingroup vtk
/// \class ttkmapFunctionTo3D
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the mapFunctionTo3D processing package.
///
/// VTK wrapping code for the @mapFunctionTo3D package.
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
/// \sa ttk::mapFunctionTo3D
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
#include                  <mapFunctionTo3D.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkmapFunctionTo3D
#else
class ttkmapFunctionTo3D
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkmapFunctionTo3D* New();
    vtkTypeMacro(ttkmapFunctionTo3D, vtkDataSetAlgorithm)
  
    
    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);

    vtkSetMacro(Mult,float);
    vtkGetMacro(Mult,float);

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
    // end of TODO-2
    
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
   
    ttkmapFunctionTo3D(){
      
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
      
    }
    
    ~ttkmapFunctionTo3D(){};
    
    TTK_SETUP();
    
    
  private:
    
    std::string           ScalarField;
    float                 Mult;
    ttk::mapFunctionTo3D  mapFunctionTo3D_;
    
};
