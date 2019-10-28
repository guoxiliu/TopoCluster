/// \ingroup base
/// \class ttk::StellarTriangulation 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %stellarTriangulation processing package.
///
/// %StellarTriangulation is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkStellarTriangulation.cpp %for a usage example.

#pragma once

// base code includes
#include                  <AbstractTriangulation.h>


namespace ttk{
  
  class StellarTriangulation : public AbstractTriangulation{

    public:

        StellarTriangulation();

        ~StellarTriangulation();


    inline int setInputCells(const SimplexId &cellNumber,
                           const LongSimplexId *cellArray){

      if(cellNumber_)
          clear();

      cellNumber_ = cellNumber;
      cellArray_ = cellArray;

      return 0;
    }

    inline int setInputPoints(const SimplexId &pointNumber, const void *pointSet,
                            const bool &doublePrecision = false){

      if(vertexNumber_)
          clear();

      vertexNumber_ = pointNumber;
      pointSet_ = pointSet;
      doublePrecision_ = doublePrecision;
      return 0;
    }



  protected:

      int clear();

      bool                doublePrecision_;
      SimplexId           cellNumber_, vertexNumber_;
      const void          *pointSet_;
      const LongSimplexId *cellArray_;

  };
}
