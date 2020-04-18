/// \ingroup vtk
/// \class ttkFG_PersistentHomology
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the fG_PersistentHomology processing package.
///
/// VTK wrapping code for the @FG_PersistentHomology package.
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
/// \sa ttk::FG_PersistentHomology
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
#include                  <FG_PersistentHomology.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFG_PersistentHomology
#else
class ttkFG_PersistentHomology
#endif
  : public vtkDataSetAlgorithm, public ttk::Wrapper{

  public:

    static ttkFG_PersistentHomology* New();
    vtkTypeMacro(ttkFG_PersistentHomology, vtkDataSetAlgorithm)

    // vtkSetMacro(ScalarField, std::string);
    // vtkGetMacro(ScalarField, std::string);

    void SetScalarField(string field){
      ScalarField = field;
      delete fG_PersistentHomology_;
      fG_PersistentHomology_ = NULL;
      Modified();
    }

    vtkSetMacro(MinPers, double);
    vtkGetMacro(MinPers, double);

    vtkSetMacro(MaxPers, double);
    vtkGetMacro(MaxPers, double);

    vtkSetMacro(cycles1, bool);
    vtkGetMacro(cycles1, bool);

    vtkSetMacro(cycles2, bool);
    vtkGetMacro(cycles2, bool);

    vtkSetMacro(formanCycles, bool);
    vtkGetMacro(formanCycles, bool);

    vtkSetMacro(Gradient, bool);
    vtkGetMacro(Gradient, bool);

    vtkSetMacro(Pairs, bool);
    vtkGetMacro(Pairs, bool);


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
        case 3:
        case 4:
        case 5:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
          break;
        default:
          break;
      }

      return 1;
    }


  protected:

    ttkFG_PersistentHomology(){
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(6);
      fG_PersistentHomology_ = NULL;
      realMaxPers=realMinPers=maxFunctValue=0;
    }

    ~ttkFG_PersistentHomology(){};

    void outHomology(vtkDataSet* output);
    void outPersistencePairs(vtkDataArray *inputScalarField, vtkDataSet* output1, vtkDataSet* output2);
    void out1Cycles(vtkDataArray *inputScalarField,vtkDataSet* output);
    void out2Cycles(vtkDataArray *inputScalarField,vtkDataSet* output);
    void outGradient(vtkDataSet* output);

    TTK_SETUP();


  private:

    std::string           ScalarField;
    bool                      cycles1;
    bool                      cycles2;
    bool                 formanCycles;
    bool                     Gradient;
    bool                        Pairs;

    double                 MaxPers;
    double                 MinPers;

    double                 realMaxPers;
    double                 realMinPers;

    double                  maxFunctValue;
    double                  minFunctionValue;

    ttk::FG_PersistentHomology*            fG_PersistentHomology_;

};
