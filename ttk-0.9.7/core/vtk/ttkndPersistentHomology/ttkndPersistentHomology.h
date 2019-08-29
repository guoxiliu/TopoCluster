/// \ingroup vtk
/// \class ttkndPersistentHomology
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ndPersistentHomology processing package.
///
/// VTK wrapping code for the @ndPersistentHomology package.
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
/// \sa ttk::ndPersistentHomology
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
#include                  <ndPersistentHomology.h>
#include                  <PersistentHomology.h>
#include                  <ttkWrapper.h>

using namespace std;

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkndPersistentHomology
#else
class ttkndPersistentHomology
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkndPersistentHomology* New();
    vtkTypeMacro(ttkndPersistentHomology, vtkDataSetAlgorithm)
   
        
    // TODO-4
    // set-getters macros to define from each variable you want to access from 
    // the outside (in particular from paraview) - to adapt.
    // Note that the XML file for the ParaView plug-in specification needs to be
    // edited accordingly.
    vtkSetMacro(Point, double);
    vtkGetMacro(Point, double);
   
    vtkSetMacro(Angle, double);
    vtkGetMacro(Angle, double);

    
    void SetSelectedFields(std::string name){
        // if(ndPersistentHomology_ != NULL){
        //   ndPersistentHomology_ = NULL;
        //   SelectedFields.clear();
        // }
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
   
    ttkndPersistentHomology(){
      
      Angle = 1;
      Point = 1;

      ndPersistentHomology_ = new ttk::ndPersistentHomology();
      
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(3);

    }
    
    ~ttkndPersistentHomology(){};
    
    TTK_SETUP();
    
    
  private:
    
    double                Point;
    double                Angle;

    vector<string>           SelectedFields;
    ttk::ndPersistentHomology*            ndPersistentHomology_;
    
};
