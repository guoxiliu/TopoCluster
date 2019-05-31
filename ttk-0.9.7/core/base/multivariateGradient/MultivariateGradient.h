/// \ingroup base
/// \class ttk::MultivariateGradient
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %multivariateGradient processing package.
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

#include <CriticalClusterGraph.h>


namespace ttk{

  typedef std::vector<std::vector<SimplexId> > SimplexesVec;
  typedef std::set<Simplex, boost::function<bool(const Simplex &, const Simplex &)>> SimplexesSet;

  typedef std::unordered_map<SimplexId, int> CMap;

  class MultivariateGradient : public Debug{

    public:

      MultivariateGradient();

      ~MultivariateGradient();

      template <typename dataType>
      int computeGradient();

      int computeCriticalClusters();
      int computeCriticalClustersAlt();


      void simplifyGraph(int maxLenght);
      float simplifyClusters(int maxClusterSize);

      inline int setInputDataPointer(std::vector<void *>* data){
        inputData_ = data;
        return 0;
      }

      inline void restoreOriginalGraph(){
        *cclusters = *backup;
      }

      inline void sortClustersBySize(){cclusters->sort();}

      void readOutputGradientPairs(std::vector<float>& pairedPoints);
      void readOutputCriticalPoints(int val, std::vector<float>& points, std::vector<char>& dimensions, std::vector<char>& clusterType, std::vector<int>& label);
      void readOutputCriticalClusters(std::list<float>& points, std::list<std::vector<int> >& cells);
      void getOutputClusterGraph(std::list<float>& points, std::list<std::pair<int,int> >& cells,std::vector<int>& csize,std::vector<char>& type,std::vector<int>& labels);

      void readOutputDescendingCell(int clIndex, std::list<float>& points, std::list<std::vector<int> >& cells, std::vector<int>& nCells);
      void readOutputAscendingCell(int clIndex, std::list<float>& points, std::list<std::vector<int> >& cells, std::vector<int>& nCells);

      void fullSetOfClusters(int cluster, unordered_set<int>& all_clusters);

      template <typename dataType>
      void readOutputCriticalCells(std::list<float>& points, std::list<char>& dimensions);

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

      inline void setThreadNumber(bool useAll, int num){
        #ifdef TTK_ENABLE_OPENMP
          if(useAll){
            threadNumber = omp_get_num_procs();
          }
          else
          {
            threadNumber = num;
          }
        #endif
      }

    protected:

    template <typename dataType>
    std::vector<dataType> getFieldValue(Simplex simplex);

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

    void mergeRegions(int region_s, int region_b,unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void unionFindEdgeVertices(unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void unionFindTriangleTetras(unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void unionFindTetraEdges(unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void unionFindTriangleVertices(unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void unionFindTetraVertices(unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void unionFindTriangleEdges(unordered_map<Simplex, int, boost::hash<Simplex> >& simplex_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl);
    void collectVertices(list<Simplex>* star, unordered_set<int>& adjacent_v);

    int extractBoundary(Simplex simpl,int dimension_b, int index_b);
    int extractCoboundary(Simplex simpl, int dimension_cb, int index_cb);
    int sizeCofacets(Simplex simpl, int dimension_cb);
    void extractOutgoingPath(Simplex simpl, CMap& critical_reached, bool saveGeometry, unordered_set<int>& geometry);
    void extractIngoingPath(Simplex simpl, unordered_set<int>& critical_reached, bool saveGeomettry, unordered_set<int>& geometry);

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
      int                    threadNumber;

      CriticalClusterGraph* cclusters;
      CriticalClusterGraph* backup;

      //input
      std::vector<void *>*        inputData_;
      Triangulation         *triangulation_;

  };

  // if the package is a pure template class, uncomment the following line
  #include                  <multivariategradient_template.h>
}
