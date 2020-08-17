/// \ingroup base
/// \class ttk::FG_Segmentation 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %fG_Segmentation processing package.
///
/// %FG_Segmentation is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkFG_Segmentation.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

#include  <FormanGradient.h>

#include <queue>
#include <unordered_set>


namespace ttk{
  
  class FG_Segmentation : public FormanGradient{

    public:

        FG_Segmentation();

        ~FG_Segmentation();

        void setInputDataPointer(void* data){inputData_ = data;};

        template <class dataType> void computeIndexing();

        void readCriticalPoints(vector<float>& points, vector<char>& dimension, vector<vector<int> >& indexes);
        void extractDescendingCell(Simplex simpl, list<Simplex>& descendingCell, set<SimplexId>& vertices);
        void extractAscendingCell(Simplex simpl, list<Simplex>& ascendingCell, set<SimplexId>& vertices);

        void computeMorseSmale();
        void markDescendingCell(Simplex simpl, boost::dynamic_bitset<>& descending3Cells);
        void markAscendingCell(Simplex simpl, boost::dynamic_bitset<>& descending3Cells);

        const void* inputData_;

  };
}


