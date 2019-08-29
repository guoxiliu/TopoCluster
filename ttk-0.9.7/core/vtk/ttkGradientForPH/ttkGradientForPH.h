/// \ingroup vtk
/// \class ttkMultivariateGradient
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the gradientForPH processing package.
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
#include                  <GradientForPH.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkGradientForPH
#else
class ttkGradientForPH
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkGradientForPH* New();
    vtkTypeMacro(ttkGradientForPH, vtkDataSetAlgorithm)

    // default ttk setters
    vtkSetMacro(debugLevel_, int);

    vtkSetMacro(extractGradient, bool);

    void SetSelectedField(std::string name){
        SelectedField = name;
    }

     int FillInputPortInformation(int port, vtkInformation *info) override
     {
         switch(port){
           case 0:
             info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
             break;
         }

         return 1;
     }

     int FillOutputPortInformation(int port, vtkInformation *info) override {
         switch(port){
           case 0: //critical cells
           case 1: //gradient pairs
              info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
             break;
         }

         return 1;
     }
    // end of TODO-3

    void outputGradient(std::vector<vtkDataSet *> &outputs);
    void outputCriticalCells(std::vector<vtkDataSet *> &outputs);

  protected:

    ttkGradientForPH(){

      gradientForPH_=NULL;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(2);
    }

    ~ttkGradientForPH(){};

    TTK_SETUP();


  private:

    std::string   SelectedField;
    bool extractGradient;


    ttk::GradientForPH*            gradientForPH_;

};
