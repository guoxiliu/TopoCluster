/// \ingroup vtk
/// \class ttkcrossDissolvePersistenceDiagrams
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the crossDissolvePersistenceDiagrams processing package.
///
/// VTK wrapping code for the @crossDissolvePersistenceDiagrams package.
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
/// \sa ttk::crossDissolvePersistenceDiagrams
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
#include                  <crossDissolvePersistenceDiagrams.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkcrossDissolvePersistenceDiagrams
#else
class ttkcrossDissolvePersistenceDiagrams
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkcrossDissolvePersistenceDiagrams* New();
    vtkTypeMacro(ttkcrossDissolvePersistenceDiagrams, vtkDataSetAlgorithm)
    
    vtkSetMacro(Alpha, double);
    vtkGetMacro(Alpha, double);

    void SetSelectedFields(std::string name){
      if(crossDissolvePersistenceDiagrams_ != NULL ){
        crossDissolvePersistenceDiagrams_ = NULL;
        SelectedFields.clear();
      }
        SelectedFields.push_back(name);
        Modified();
    }

    void outHomology(PersistentHomology* ph, vtkDataSet* output);
    void outPersistencePairs(PersistentHomology* ph, vtkDataSet* output);

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
        case 1:
        case 2:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
          break;
        default:
          break;
      }
      
      return 1;
    }
    
    
  protected:
   
    ttkcrossDissolvePersistenceDiagrams(){
      
      Alpha = 1;
      
      crossDissolvePersistenceDiagrams_ = NULL;
      
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(3);
    }
    
    ~ttkcrossDissolvePersistenceDiagrams(){};
    
    TTK_SETUP();
    
    
  private:
    
    double                Alpha;

    vector<string>           SelectedFields;
    ttk::crossDissolvePersistenceDiagrams*            crossDissolvePersistenceDiagrams_;
    
};
