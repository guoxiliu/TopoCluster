/// \ingroup base
/// \class ttk::MultivariateGradient
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %gradientForPH processing package.
///
/// %MultivariateGradient is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkMultivariateGradient.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/functional/hash.hpp>
#include <boost/optional/optional.hpp>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;

namespace ttk{
  typedef std::pair<short int,SimplexId> Simplex;
  typedef std::vector<std::vector<SimplexId> > SimplexesVec;
  typedef std::set<Simplex, boost::function<bool(const Simplex &, const Simplex &)>> SimplexesSet;

  typedef std::unordered_map<SimplexId, int> CMap;

  class GradientForPH : public Debug{

    public:

      GradientForPH();

      ~GradientForPH();

      template <typename dataType>
      int computeGradient();

      inline int setInputDataPointer(void* data){
        inputData_ = data;
        return 0;
      }

      void readOutputGradientPairs(std::vector<float>& pairedPoints);
      void readOutputCriticalPoints(int val, std::vector<float>& points, std::vector<char>& dimensions);

      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
        dimensionality_ = triangulation_->getCellVertexNumber(0) - 1;

        if(triangulation_){
            triangulation_->preprocessVertexStars();
            triangulation_->preprocessVertexEdges();
            triangulation_->preprocessEdges();

            if (dimensionality_ >= 2) {
              triangulation_->preprocessTriangles();
              triangulation_->preprocessEdgeTriangles();
              triangulation_->preprocessTriangleEdges();
            }

            if (dimensionality_ == 3) {
              triangulation_->preprocessCellTriangles();
              triangulation_->preprocessTriangleStars();
              triangulation_->preprocessVertexTriangles();
            }

        }

        return 0;
      }

    protected:

    template <typename dataType>
    dataType getFieldValue(Simplex simplex);

    template <typename dataType>
    void computeIndexing();

    template <typename dataType>
    int assignGradient(SimplexesVec& vec);

    template <typename dataType>
    int numPairableFaces(Simplex simplex,
                         SimplexesSet& pairable,
                         Simplex& chosen);


    SimplexId getIndexing(Simplex);
    void getVertexIndexingOfVertices(Simplex simplex, std::vector<SimplexId>& vec);
    void simplexToVertices(Simplex, std::vector<SimplexId>&);
    void computeLowerStar(SimplexId, SimplexesVec&);
    bool cmpSimplexesLevelSets(const Simplex& s1, const Simplex& s2);

    
    int extractBoundary(Simplex simpl,int dimension_b, int index_b);
    int extractCoboundary(Simplex simpl, int dimension_cb, int index_cb);
    int sizeCofacets(Simplex simpl, int dimension_cb);
    
    void setPair(const Simplex& tail, const Simplex& head);
    bool getPair(const Simplex& simpl, Simplex& paired);
    bool isCritical(const Simplex& simpl);
    bool isTopCritical(const Simplex& simpl);

    void computeBarycenter(const Simplex& simplex, std::vector<float>& coords);
    float computeDistance(const vector<float>& point1, const vector<float>& point2);



    protected:
      std::vector<SimplexId>*       vertexIndexing;
      std::vector<std::vector<std::vector<SimplexId> > >*   gradient;

      int                    dimensionality_;
      int                    numberCriticalSimplices;

      //input
      void*        inputData_;
      Triangulation         *triangulation_;

  };

  // if the package is a pure template class, uncomment the following line
  #include                  <multivariategradient_template.h>
}
