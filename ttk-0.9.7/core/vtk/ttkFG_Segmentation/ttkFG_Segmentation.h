/// \ingroup vtk
/// \class ttkFG_Segmentation
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the fG_Segmentation processing package.
///
/// VTK wrapping code for the @FG_Segmentation package.
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
/// \sa ttk::FG_Segmentation
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
#include                  <vtkCellData.h>

// ttk code includes
#include                  <FG_Segmentation.h>
#include                  <FG_Segmentation_template.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFG_Segmentation
#else
class ttkFG_Segmentation
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:
    
    static ttkFG_Segmentation* New();
    vtkTypeMacro(ttkFG_Segmentation, vtkDataSetAlgorithm)
    
    vtkSetMacro(ScalarField, std::string);
    vtkGetMacro(ScalarField, std::string);

    vtkSetMacro(CacheSize, int);
    vtkGetMacro(CacheSize, int);

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
   
    ttkFG_Segmentation(){
    
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(3);

      fG_Segmentation_ = new ttk::FG_Segmentation();
    }

    ~ttkFG_Segmentation(){};

    void outputCriticalCells(vector<vtkDataSet*> output, Triangulation* triangulation);
    

    
    TTK_SETUP();
    
    
  private:
    
    std::string ScalarField;
    ttk::FG_Segmentation* fG_Segmentation_;
    int CacheSize;
    
};
