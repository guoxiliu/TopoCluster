/// \ingroup base
/// \class ttk::StellarTriangulation 
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Nov. 2019.
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
#include                  <algorithm>
#include                  <list>
#include                  <boost/unordered_map.hpp>
#include                  <boost/unordered_set.hpp>
#include                  <AbstractTriangulation.h>
#include                  <boost/functional/hash.hpp>

#define VERTEX_ID 0
#define EDGE_ID 1
#define TRIANGLE_ID 2
#define CELL_ID 3
#define pair_int std::pair<ttk::SimplexId, ttk::SimplexId>

using namespace std;
namespace ttk{

  class ExpandedNode
  {
  private:
    /* components */
    SimplexId nid;
    vector<pair_int> *internalEdgeList_;
    vector<vector<SimplexId>> *internalTriangleList_;
    boost::unordered_map<pair_int, SimplexId> *internalEdgeMap_;
    boost::unordered_map<pair_int, SimplexId> *externalEdgeMap_;
    boost::unordered_map<vector<SimplexId>, SimplexId> *internalTriangleMap_;
    boost::unordered_map<vector<SimplexId>, SimplexId> *externalTriangleMap_;
    /* vertex relationships */
    vector<vector<SimplexId>> *vertexEdges_;
    vector<vector<SimplexId>> *vertexTriangles_;
    vector<vector<SimplexId>> *vertexStars_;
    vector<vector<SimplexId>> *vertexNeighbors_;
    /* edge relationships */
    // edgeVertex relation can be extracted from internal edge list
    vector<vector<SimplexId>> *edgeTriangles_;
    vector<vector<SimplexId>> *edgeStars_;
    /* triangle relationships */
    // triangleVertex relation can be extracted from internal triangle list
    vector<vector<SimplexId>> *triangleEdges_;
    vector<vector<SimplexId>> *triangleStars_;
    /* cell relationships */
    vector<vector<SimplexId>> *cellEdges_;
    vector<vector<SimplexId>> *cellTriangles_;

  public:
    ExpandedNode(SimplexId id){
      /* components */
      nid = id;
      internalEdgeList_ = nullptr;
      internalTriangleList_ = nullptr;
      internalEdgeMap_ = nullptr;
      externalEdgeMap_ = nullptr;
      internalTriangleMap_ = nullptr;
      externalTriangleMap_ = nullptr;
      /* vertex relationships */
      vertexEdges_ = nullptr;
      vertexTriangles_ = nullptr;
      vertexStars_ = nullptr;
      /* edge relationships */
      edgeTriangles_ = nullptr;
      edgeStars_ = nullptr;
      /* triangle relationships */
      triangleEdges_ = nullptr;
      triangleStars_ = nullptr;
      /* cell relationships */
      cellEdges_ = nullptr;
      cellTriangles_ = nullptr;
    }
    ~ExpandedNode(){
      delete internalEdgeList_;
      delete internalTriangleList_;
      delete internalEdgeMap_;
      delete externalEdgeMap_;
      delete internalTriangleMap_;
      delete externalTriangleMap_;
      delete vertexEdges_;
      delete vertexTriangles_;
      delete vertexStars_;
      delete edgeTriangles_;
      delete edgeStars_;
      delete triangleEdges_;
      delete triangleStars_;
      delete cellEdges_;
      delete cellTriangles_;
    }

    friend class StellarTriangulation;
  };

  
  class StellarTriangulation : public AbstractTriangulation{

    public:

      StellarTriangulation();

      ~StellarTriangulation();


      int setInputCells(const SimplexId &cellNumber,
                            const LongSimplexId *cellArray){

        if(cellNumber_)
          clear();

        cellNumber_ = cellNumber;
        cellArray_ = cellArray;

        if((!cellArray_)||(!cellNumber_))
          return -1;
        
        // initialize the array of cell intervals
        cellIntervals_.resize(nodeNumber_+1);
        externalCells_.resize(nodeNumber_+1);
        cellIntervals_[0] = -1;

        SimplexId nodeNum = 1, cid = 0;
        vector<SimplexId> cell;

        while(nodeNum <= nodeNumber_){
          SimplexId startPos = (cellArray[0]+1)*cid+1;
          cell = vector<SimplexId>(cellArray+startPos, cellArray+startPos+cellArray[0]);

          if(cell[0] > vertexIntervals_[nodeNum] || cid >= cellNumber){
            cellIntervals_[nodeNum++] = cid - 1;
            continue;
          }
          // create external cell list 
          for(size_t i = 1; i < cell.size(); i++){
            if(cell[i] > vertexIntervals_[nodeNum]){
            SimplexId nid = findNodeIndex(cell[i], VERTEX_ID);
              if(externalCells_[nid].empty()){
                externalCells_[nid].push_back(cid);
              }else if(externalCells_[nid].back() != cid){
                externalCells_[nid].push_back(cid);
              }
            }
          }
          if(cid < cellNumber)
            cid++;
        }

        return 0;
      }

      int setInputPoints(const SimplexId &pointNumber, const void *pointSet, 
        const int *indexArray, const bool &doublePrecision = false){

        if(vertexNumber_)
            clear();

        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
        vertexIndices_ = indexArray;
        doublePrecision_ = doublePrecision;

        // initialize the array of vertex intervals
        // add a dummy node at the beginning
        vertexIntervals_.push_back(-1);

        SimplexId vid = 1;
        for(; vid < pointNumber; vid++){
          if(indexArray[vid] != indexArray[vid-1]){
            vertexIntervals_.push_back(vid-1);
          }
        }
        vertexIntervals_.push_back(vid-1);
        nodeNumber_ = vertexIntervals_.size()-1;
        initCache();

        cout << "[StellarTriangulation] node num: " << nodeNumber_ << ", cache size: " << cacheSize_ << endl;

        return 0;
      }

      int getCellEdge(const SimplexId &cellId, 
        const int &localEdgeId, SimplexId &edgeId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
          if((localEdgeId < 0))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(cellId, CELL_ID);
        SimplexId localCellId = cellId - cellIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        
        if(exnode->cellEdges_ == nullptr){
          exnode->cellEdges_ = new vector<vector<SimplexId>>();
          getCellEdges(exnode);
        }
        if(localEdgeId >= (SimplexId) exnode->cellEdges_->at(localCellId).size())
          return -2;
        edgeId = exnode->cellEdges_->at(localCellId).at(localEdgeId);
        return 0;
      }
        
      inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const{

       #ifndef TTK_ENABLE_KAMIKAZE
        if((cellId < 0)||(cellId >= cellNumber_))
          return -1;
        #endif

        return cellArray_[0] * (cellArray_[0]-1) / 2;
      }
      
      const vector<vector<SimplexId> > *getCellEdges(){
        cellEdgeList_.reserve(cellNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->cellEdges_ == nullptr){
            exnode->cellEdges_ = new vector<vector<SimplexId>>();
            getCellEdges(exnode);
          }
          cellEdgeList_.insert(cellEdgeList_.end(), exnode->cellEdges_->begin(), exnode->cellEdges_->end());
        }
        return &cellEdgeList_;
      }
      
      int getCellNeighbor(const SimplexId &cellId,
        const int &localNeighborId, SimplexId &neighborId) const{
        return 0;
      }
        
      SimplexId getCellNeighborNumber(const SimplexId &cellId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getCellNeighbors(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getCellTriangle(const SimplexId &cellId, 
        const int &localTriangleId, SimplexId &triangleId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((cellId < 0)||(cellId >= cellNumber_))
          return -1;
        if((localTriangleId < 0)||(localTriangleId >= getCellTriangleNumber(cellId)))
          return -2;
        #endif

        SimplexId nid = findNodeIndex(cellId, CELL_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->cellTriangles_ == nullptr){
          exnode->cellTriangles_ = new vector<vector<SimplexId>>();
          getCellTriangles(exnode);
        }
        triangleId = exnode->cellTriangles_->at(cellId-cellIntervals_[nid-1]-1)[localTriangleId];
        return 0;
      }
        
      SimplexId getCellTriangleNumber(const SimplexId &cellId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((cellId < 0)||(cellId >= cellNumber_))
          return -1;
        #endif

        return cellArray_[0] * (cellArray_[0]-1) * (cellArray_[0]-2) / 6;
      }
        
      const vector<vector<SimplexId> > *getCellTriangles(){
        cellTriangleList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->cellTriangles_ == nullptr){
            exnode->cellTriangles_ = new vector<vector<SimplexId>>();
            getCellTriangles(exnode);
          }
          cellTriangleList_.insert(cellTriangleList_.end(), exnode->cellTriangles_->begin(), exnode->cellTriangles_->end());
        }
        return &cellTriangleList_;
      }
      
      inline int getCellVertex(const SimplexId &cellId,
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
          if((localVertexId < 0)
            ||(localVertexId >= cellArray_[0]))
            return -2;
        #endif

        vertexId = cellArray_[(cellArray_[0] + 1)*cellId + localVertexId + 1];
        return 0;
      }
    
      inline SimplexId getCellVertexNumber(const SimplexId &cellId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
        #endif
        
        return cellArray_[0];
      }
        
      inline int getDimensionality() const{
        if((cellArray_)&&(cellNumber_)){
          return cellArray_[0] - 1;
        }
        return -1;
      }
      
      const vector<pair_int > *getEdges(){
        edgeList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          vector<pair_int> localInternalEdges;
          boost::unordered_map<pair_int, SimplexId> localInternalEdgeMap;
          buildInternalEdgeList(nid, &localInternalEdges, &localInternalEdgeMap, nullptr);
          edgeList_.insert(edgeList_.end(), localInternalEdges.begin(), localInternalEdges.end());
        }
        return &edgeList_;
      }
      
      int getEdgeLink(const SimplexId &edgeId, 
        const int &localLinkId, SimplexId &linkId) const{
        return 0;
      }
        
      SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getEdgeLinks(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getEdgeStar(const SimplexId &edgeId, 
        const int &localStarId, SimplexId &starId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localEdgeId = edgeId - edgeIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->edgeStars_ == nullptr){
          exnode->edgeStars_ = new vector<vector<SimplexId>>();
          buildInternalEdgeList(nid, nullptr, nullptr, exnode->edgeStars_);
        }
        if(localStarId >= (SimplexId) (*(exnode->edgeStars_))[localEdgeId].size())
          return -2;
        starId = (*(exnode->edgeStars_))[localEdgeId][localStarId];
        return 0;
      }
        
      SimplexId getEdgeStarNumber(const SimplexId &edgeId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        ExpandedNode *exnode = searchCache(nid);
        SimplexId localEdgeId = edgeId - edgeIntervals_[nid-1] - 1;
        if(exnode->edgeStars_ == nullptr){
          exnode->edgeStars_ = new vector<vector<SimplexId>>();
          buildInternalEdgeList(nid, nullptr, nullptr,  exnode->edgeStars_);
        }
        return (*(exnode->edgeStars_))[localEdgeId].size();
      }
      
      const vector<vector<SimplexId> > *getEdgeStars(){
        edgeStarList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->edgeStars_ == nullptr){
            exnode->edgeStars_ = new vector<vector<SimplexId>>();
            buildInternalEdgeList(nid, nullptr, nullptr, exnode->edgeStars_);
          }
          edgeStarList_.insert(edgeStarList_.end(), exnode->edgeStars_->begin(), exnode->edgeStars_->end());
        }
        return &edgeStarList_;
      }
     
      int getEdgeTriangle(const SimplexId &edgeId, 
        const int &localTriangleId, SimplexId &triangleId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.back()))
          return -1;
        if(localTriangleId < 0)
          return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->edgeTriangles_ == nullptr){
          exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
          getEdgeTriangles(exnode);
        }
        if(localTriangleId >= (SimplexId) (*(exnode->edgeTriangles_))[localEdgeId].size())
          return -2;
        triangleId = (*(exnode->edgeTriangles_))[localEdgeId][localTriangleId];
        return 0;
      }
      
      SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.back()))
          return -1;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->edgeTriangles_ == nullptr){
          exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
          getEdgeTriangles(exnode);
          
        }
        return exnode->edgeTriangles_->at(edgeId-edgeIntervals_[nid-1]-1).size();
      }
      
      const vector<vector<SimplexId> > *getEdgeTriangles(){
        edgeTriangleList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->edgeTriangles_ == nullptr){
            exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
            getEdgeTriangles(exnode);

          }
          edgeTriangleList_.insert(edgeTriangleList_.end(), exnode->edgeTriangles_->begin(), exnode->edgeTriangles_->end());
        }
        return &edgeTriangleList_;
      }
      
      int getEdgeVertex(const SimplexId &edgeId, 
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.back()))
            return -1;
          if((localVertexId != 0)&&(localVertexId != 1))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->internalEdgeList_ == nullptr){
          exnode->internalEdgeList_ = new vector<pair_int>();
          exnode->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildInternalEdgeList(nid, exnode->internalEdgeList_, exnode->internalEdgeMap_, nullptr);
        }
        if(localVertexId){
          vertexId = exnode->internalEdgeList_->at(localEdgeId).second;
        }
        else{
          vertexId = exnode->internalEdgeList_->at(localEdgeId).first;
        }
        return 0;
      }
      
      inline SimplexId getNumberOfCells() const { return cellNumber_; }
      
      inline SimplexId getNumberOfEdges() const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(!edgeIntervals_.size())
            return -1;
        #endif

        return edgeIntervals_.back()+1;
      }

      SimplexId getNumberOfTriangles() const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(!triangleIntervals_.size())
            return -1;
        #endif

        return triangleIntervals_.back()+1;
      }
      
      inline SimplexId getNumberOfVertices() const { return vertexNumber_; }
      
      const vector<vector<SimplexId> > *getTriangles(){
        triangleList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->internalTriangleList_ == nullptr){
            exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
            exnode->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
            buildInternalTriangleList(nid, exnode->internalTriangleList_, exnode->internalTriangleMap_, nullptr);
          }
          triangleList_.insert(triangleList_.end(), exnode->internalTriangleList_->begin(), exnode->internalTriangleList_->end());
        }
        return &triangleList_;
      }
      
      int getTriangleEdge(const SimplexId &triangleId,
        const int &localEdgeId, SimplexId &edgeId) const{
                    
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if((localEdgeId < 0)||(localEdgeId > 2))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleEdges_ == nullptr){
          exnode->triangleEdges_ = new vector<vector<SimplexId>>();
          getTriangleEdges(exnode);
        }
        edgeId = (*(exnode->triangleEdges_))[triangleId-triangleIntervals_[nid-1]-1][localEdgeId];
        return 0;
      }
      
      SimplexId getTriangleEdgeNumber(const SimplexId &triangleId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleEdges_ == nullptr){
          exnode->triangleEdges_ = new vector<vector<SimplexId>>();
          getTriangleEdges(exnode);
        }
        return (*(exnode->triangleEdges_))[triangleId-triangleIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleEdges(){
        triangleEdgeList_.reserve(triangleIntervals_.size()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->triangleEdges_ == nullptr){
            exnode->triangleEdges_ = new vector<vector<SimplexId>>();
            getTriangleEdges(exnode);
          }
          triangleEdgeList_.insert(triangleEdgeList_.end(), exnode->triangleEdges_->begin(), exnode->triangleEdges_->end());
        }
        return &triangleEdgeList_;
      }
      
      int getTriangleLink(const SimplexId &triangleId, 
        const int &localLinkId, SimplexId &linkId) const{
          return 0;
        }
        
      SimplexId getTriangleLinkNumber(const SimplexId &triangleId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getTriangleLinks(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getTriangleStar(const SimplexId &triangleId,
        const int &localStarId, SimplexId &starId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleStars_ == nullptr){
          exnode->triangleStars_ = new vector<vector<SimplexId>>();
          buildInternalTriangleList(nid, nullptr, nullptr, exnode->triangleStars_);
        }
        if(localStarId >= (SimplexId) exnode->triangleStars_->at(localTriangleId).size())
          return -2;
        
        starId = exnode->triangleStars_->at(localTriangleId)[localStarId];
        return 0;
      }
        
      SimplexId getTriangleStarNumber(const SimplexId &triangleId) const{

         #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleStars_ == nullptr){
          exnode->triangleStars_ = new vector<vector<SimplexId>>();
          buildInternalTriangleList(nid, nullptr, nullptr, exnode->triangleStars_);
        }
        
        return (*(exnode->triangleStars_))[triangleId-triangleIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleStars(){
        triangleStarList_.reserve(triangleIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->triangleStars_ == nullptr){
            exnode->triangleStars_ = new vector<vector<SimplexId>>();
            buildInternalTriangleList(nid, nullptr, nullptr, exnode->triangleStars_);
          }
          triangleStarList_.insert(triangleStarList_.end(), exnode->triangleStars_->begin(), exnode->triangleStars_->end());
        }
        return &triangleStarList_;
      }
      
      int getTriangleVertex(const SimplexId &triangleId,
        const int &localVertexId, SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if((localVertexId < 0)||(localVertexId > 2))
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->internalTriangleList_ == nullptr){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          exnode->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildInternalTriangleList(nid, exnode->internalTriangleList_, exnode->internalTriangleMap_, nullptr);
        }
        vertexId = exnode->internalTriangleList_->at(triangleId-triangleIntervals_[nid-1]-1)[localVertexId];
        return 0;
      }
      
      int getVertexEdge(const SimplexId &vertexId, 
        const int &localEdgeId, SimplexId &edgeId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localEdgeId < 0)
            return -2;
        #endif
        
        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexEdges_ == nullptr){
          exnode->vertexEdges_ = new vector<vector<SimplexId>>();
          getVertexEdges(exnode);
        }
        if(localEdgeId >= (SimplexId) (*exnode->vertexEdges_)[localVertexId].size())
            return -2;
        edgeId = exnode->vertexEdges_->at(localVertexId)[localEdgeId];
        return 0;
      }
      
      SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexEdges_ == nullptr){
          exnode->vertexEdges_ = new vector<vector<SimplexId>>();
          getVertexEdges(exnode);
        }
        return (*(exnode->vertexEdges_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexEdges(){
        vertexEdgeList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->vertexEdges_ == nullptr){
            exnode->vertexEdges_ = new vector<vector<SimplexId>>();
            getVertexEdges(exnode);
          }
          vertexEdgeList_.insert(vertexEdgeList_.end(), exnode->vertexEdges_->begin(), exnode->vertexEdges_->end());
        }

        return &vertexEdgeList_;
      }
      
      int getVertexLink(const SimplexId &vertexId, 
        const int &localLinkId, SimplexId &linkId) const{
        return 0;
      }
      
      SimplexId getVertexLinkNumber(const SimplexId &vertexId) const{
        return 0;
      }
      
      const vector<vector<SimplexId> > *getVertexLinks(){
        vector<vector<SimplexId>> *dummyPointer= new vector<vector<SimplexId>>();
        return dummyPointer;
      }
      
      int getVertexNeighbor(const SimplexId &vertexId, 
        const int &localNeighborId, SimplexId &neighborId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localNeighborId < 0)
            return -2;
        #endif
        
        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexNeighbors_ == nullptr){
          exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
          getVertexNeighbors(exnode);
        }
        if(localNeighborId < (SimplexId) (*(exnode->vertexNeighbors_))[localVertexId].size())
          return -2;  
        neighborId = (*(exnode->vertexNeighbors_))[localVertexId][localNeighborId];
        return 0;
      }
      
      SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexNeighbors_ == nullptr){
          exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
          getVertexNeighbors(exnode);
        }
        return (*(exnode->vertexNeighbors_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexNeighbors(){
        vertexNeighborList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->vertexNeighbors_ == nullptr){
            exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
            getVertexNeighbors(exnode);
          }
          vertexNeighborList_.insert(vertexNeighborList_.end(), exnode->vertexNeighbors_->begin(), exnode->vertexNeighbors_->end());
        }
        return &vertexNeighborList_;
      }
      
      int getVertexPoint(const SimplexId &vertexId,
        float &x, float &y, float &z) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif
        
        if(doublePrecision_){
          x = ((double *) pointSet_)[3*vertexId];
          y = ((double *) pointSet_)[3*vertexId + 1];
          z = ((double *) pointSet_)[3*vertexId + 2];
        }
        else{
          x = ((float *) pointSet_)[3*vertexId];
          y = ((float *) pointSet_)[3*vertexId + 1];
          z = ((float *) pointSet_)[3*vertexId + 2];
        }
      
        return 0;
      }
      
      int getVertexStar(const SimplexId &vertexId, const int &localStarId,
        SimplexId &starId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localStarId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexStars_ == nullptr){
          exnode->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(exnode);
        }
        if(localStarId >= (SimplexId) (*exnode->vertexStars_)[localVertexId].size())
          return -2;
        starId = (*exnode->vertexStars_)[localVertexId][localStarId];
        return 0;
      }
      
      SimplexId getVertexStarNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexStars_ == nullptr){
          exnode->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(exnode);
        }
        return (*(exnode->vertexStars_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexStars(){
        vertexStarList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->vertexStars_ == nullptr){
            exnode->vertexStars_ = new vector<vector<SimplexId>>();
            getVertexStars(exnode);
          }
          vertexStarList_.insert(vertexStarList_.end(), exnode->vertexStars_->begin(), exnode->vertexStars_->end());
        }
        return &vertexStarList_;
      }

      int getVertexTriangle(const SimplexId &vertexId, 
        const int &localTriangleId, SimplexId &triangleId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localTriangleId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        SimplexId localVertexId = vertexId-vertexIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexTriangles_ == nullptr){
          exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
          getVertexTriangles(exnode);
        }
        if(localTriangleId >= (SimplexId) (*(exnode->vertexTriangles_))[localVertexId].size())
          return -2;
        triangleId = (*(exnode->vertexTriangles_))[localVertexId][localTriangleId];

        return 0;
      }
      
      SimplexId getVertexTriangleNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(vertexId, VERTEX_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexTriangles_ == nullptr){
          exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
          getVertexTriangles(exnode);
        }
        return (*(exnode->vertexTriangles_))[vertexId-vertexIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getVertexTriangles(){
        vertexTriangleList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->vertexTriangles_ == nullptr){
            exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
            getVertexTriangles(exnode);
          }
          vertexTriangleList_.insert(vertexTriangleList_.end(), exnode->vertexTriangles_->begin(), exnode->vertexTriangles_->end());
        }
        return &vertexTriangleList_;
      }
        
      inline bool hasPreprocessedBoundaryEdges() const{
        return hasPreprocessedBoundaryEdges_;
      }
      
      inline bool hasPreprocessedBoundaryTriangles() const{
        return hasPreprocessedBoundaryTriangles_;
      }

      inline bool hasPreprocessedBoundaryVertices() const{
        return hasPreprocessedBoundaryVertices_;
      }
     
      inline bool hasPreprocessedCellEdges() const{
        return hasPreprocessedCellEdges_;
      }
      
      inline bool hasPreprocessedCellNeighbors() const{
        return hasPreprocessedCellNeighbors_;
      }
      
      inline bool hasPreprocessedCellTriangles() const{
        return hasPreprocessedCellTriangles_;
      }
       
      inline bool hasPreprocessedEdgeLinks() const{
        return hasPreprocessedEdgeLinks_;
      }
       
      inline bool hasPreprocessedEdgeStars() const{
        return hasPreprocessedEdgeStars_;
      }
      
      inline bool hasPreprocessedEdgeTriangles() const{
        return hasPreprocessedEdgeTriangles_;
      }
        
      inline bool hasPreprocessedEdges() const{
        return hasPreprocessedEdges_;
      }
      
      inline bool hasPreprocessedTriangles() const{
        return hasPreprocessedTriangles_;
      }
      
      inline bool hasPreprocessedTriangleEdges() const{
        return hasPreprocessedTriangleEdges_;
      }
      
      inline bool hasPreprocessedTriangleLinks() const{
        return hasPreprocessedTriangleLinks_;
      }
      
      inline bool hasPreprocessedTriangleStars() const{
        return hasPreprocessedTriangleStars_;
      }
        
      inline bool hasPreprocessedVertexEdges() const{
        return hasPreprocessedVertexEdges_;
      }
      
      inline bool hasPreprocessedVertexLinks() const{
        return hasPreprocessedVertexLinks_;
      }
      
      inline bool hasPreprocessedVertexNeighbors() const{
        return hasPreprocessedVertexNeighbors_;
      }
      
      inline bool hasPreprocessedVertexStars() const{
        return hasPreprocessedVertexStars_;
      }
      
      inline bool hasPreprocessedVertexTriangles() const{
        return hasPreprocessedVertexTriangles_;
      }
      
      bool isEdgeOnBoundary(const SimplexId &edgeId) const{
        return false;
      }
        
      bool isEmpty() const{
        return false;
      }
      
      bool isTriangleOnBoundary(const SimplexId &triangleId) const{
        return false;
      }
      
      bool isVertexOnBoundary(const SimplexId &vertexId) const{
        return false;
      }

      int preprocessBoundaryEdges(){
        preprocessEdges();
        hasPreprocessedBoundaryEdges_ = true;
        return 0;
      }
      
      int preprocessBoundaryTriangles(){
        preprocessTriangles();
        hasPreprocessedBoundaryTriangles_ = true;
        return 0;
      }
      
      int preprocessBoundaryVertices(){
        hasPreprocessedBoundaryVertices_ = true;
        return 0;
      }
      
      int preprocessCellEdges(){
        preprocessEdges();
        hasPreprocessedCellEdges_ = true;
        return 0;
      }
      
      int preprocessCellNeighbors(){
        hasPreprocessedCellNeighbors_ = true;
        return 0;
      }
      
      int preprocessCellTriangles(){
        preprocessTriangles();
        hasPreprocessedCellTriangles_ = true;
        return 0;
      }
      
      int preprocessEdges(){

        #ifndef TTK_ENABLE_KAMIKAZE
          if(vertexNumber_ <= 0)
            return -1;
          if(cellNumber_ <= 0)
            return -2;
          if(!cellArray_)
            return -3;
          if(nodeNumber_ <= 0)
            return -4;
        #endif

        if(!hasPreprocessedEdges_){
          edgeIntervals_.resize(nodeNumber_+1);
          edgeIntervals_[0] = -1;
          vector<SimplexId> edgeCount(nodeNumber_+1);
          #ifdef TTK_ENABLE_OPENMP
          #pragma omp parallel for num_threads(threadNumber_)
          #endif
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            edgeCount[nid] = countInternalEdges(nid);
          }

          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            edgeIntervals_[nid] = edgeIntervals_[nid-1]+edgeCount[nid];
          }

          hasPreprocessedEdges_ = true;
        }

        return 0;
      }
      
      int preprocessEdgeLinks(){
        preprocessEdges();
        hasPreprocessedEdgeLinks_ = true;
        return 0;
      }
      
      int preprocessEdgeStars(){
        preprocessEdges();
        hasPreprocessedEdgeStars_ = true;
        return 0;
      }
      
      int preprocessEdgeTriangles(){
        preprocessEdges();
        preprocessTriangles();
        hasPreprocessedEdgeTriangles_ = true;
        return 0;
      }
      
      int preprocessTriangles(){

        #ifndef TTK_ENABLE_KAMIKAZE
          if(vertexNumber_ <= 0)
            return -1;
          if(cellNumber_ <= 0)
            return -2;
          if(!cellArray_)
            return -3;
          if(nodeNumber_ <= 0)
            return -4;
        #endif

        // build triangle interval list
        if(!hasPreprocessedTriangles_){
          triangleIntervals_.resize(nodeNumber_+1);
          triangleIntervals_[0] = -1;
          vector<SimplexId> triangleCount(nodeNumber_+1);
          #ifdef TTK_ENABLE_OPENMP
          #pragma omp parallel for num_threads(threadNumber_)
          #endif
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            triangleCount[nid] = countInternalTriangles(nid);
          }

          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            triangleIntervals_[nid] = triangleIntervals_[nid-1]+triangleCount[nid];
          }

          hasPreprocessedTriangles_ = true;
        }

        return 0;
      }
      
      int preprocessTriangleEdges(){
        preprocessEdges();
        preprocessTriangles();
        hasPreprocessedTriangleEdges_ = true;
        return 0;
      }
      
      int preprocessTriangleLinks(){
        preprocessTriangles();
        hasPreprocessedTriangleLinks_ = true;
        return 0;
      }
      
      int preprocessTriangleStars(){
        preprocessTriangles();
        hasPreprocessedTriangleStars_ = true;
        return 0;
      }
      
      int preprocessVertexEdges() { 
        preprocessEdges();
        hasPreprocessedVertexEdges_ = true;
        return 0;
      }
      
      int preprocessVertexLinks(){
        hasPreprocessedVertexLinks_ = true;
        return 0;
      }
      
      int preprocessVertexNeighbors(){
        preprocessEdges();
        hasPreprocessedVertexNeighbors_ = true;
        return 0;
      }
      
      int preprocessVertexStars(){

        #ifndef TTK_ENABLE_KAMIKAZE
          if(!cellArray_)
            return -1;
        #endif
        
        vertexStarList_.resize(vertexNumber_);
        for(SimplexId i = 0; i < cellNumber_; i++){
          for(SimplexId j = 0; j < cellArray_[0];j++){
            vertexStarList_[cellArray_[(cellArray_[0] + 1)*i + 1 + j]].push_back(i);
          }
        }

        hasPreprocessedVertexStars_ = true;
        return 0;
      }
      
      int preprocessVertexTriangles(){
        preprocessTriangles();
        hasPreprocessedVertexTriangles_ = true;
        return 0;
      }


      /**
       * Initialize the cache.
       */
      void initCache(const size_t size=100){
        cacheSize_ = size;
      }

    protected:

      int clear();

      /**
       * Find the corresponding node index given the id.
       */ 
      SimplexId findNodeIndex(SimplexId id, int idType) const{
        if(idType == VERTEX_ID){
          return vertexIndices_[id]+1;
        }
        const vector<SimplexId> *intervals = nullptr;
        // determine which vector to search
        if(idType == EDGE_ID){
          intervals = &edgeIntervals_;
        }else if(idType == TRIANGLE_ID){
          intervals = &triangleIntervals_;
        }else if(idType == CELL_ID){
          intervals = &cellIntervals_;
        }else{
          return -1;
        }

        vector<SimplexId>::const_iterator low = lower_bound(intervals->begin(), intervals->end(), id);
        return (low-intervals->begin());
      }

      /**
       * Search the node in the cache.
       */
      ExpandedNode* searchCache(const SimplexId &nodeId, const SimplexId reservedId=0) const{
        // cannot find the expanded node in the cache
        if(cacheMap_.find(nodeId) == cacheMap_.end()){
          // missCount_++;
          if(cache_.size() >= cacheSize_){
            if(cache_.back()->nid == reservedId){
              return nullptr;
            }
            cacheMap_.erase(cache_.back()->nid);
            delete cache_.back();
            cache_.pop_back();
          }
          cache_.push_front(new ExpandedNode(nodeId));
          cacheMap_[nodeId] = cache_.begin();
        }
        return (*cacheMap_[nodeId]);
      }

      /**
       * Insert a vector of nodes in the cache.
       */
      void insertCache(const vector<ExpandedNode*> nodeVec, const SimplexId reservedId=0) const{
        for(ExpandedNode *exnode : nodeVec){
          if(cacheMap_.find(exnode->nid) != cacheMap_.end())
            continue;
          if(cache_.size() >= cacheSize_){
            if(cache_.back()->nid == reservedId)
              break;
            cacheMap_.erase(cache_.back()->nid);
            delete cache_.back();
            cache_.pop_back();
          }
          cache_.push_front(exnode);
          cacheMap_[exnode->nid] = cache_.begin();
        }
      }

      /** 
       * Build the internal edge list in the node.
       */ 
      int buildInternalEdgeList(SimplexId nodeId, 
        vector<pair_int> * const internalEdgeList,
        boost::unordered_map<pair_int, SimplexId> * const internalEdgeMap,
        vector<vector<SimplexId>> * const edgeStars) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId edgeCount = 0, verticesPerCell = cellArray_[0];
        boost::unordered_map<pair_int, SimplexId> *edgeMap;
        // loop through all the internal cells first
        if(edgeStars){
          edgeStars->clear();
          edgeStars->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(internalEdgeMap){
          internalEdgeMap->clear();
          edgeMap = internalEdgeMap;
        }else{
          edgeMap = new boost::unordered_map<pair_int, SimplexId>();
        }

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          pair_int edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            edgeIds.first = cellArray_[cellId + j + 1];
            // the edge does not belong to the current node
            if(edgeIds.first > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // not found in the edge map - assign new edge id
              if(edgeMap->find(edgeIds) == edgeMap->end()){
                edgeCount++;
                (*edgeMap)[edgeIds] = edgeCount+edgeIntervals_[nodeId-1];
                if(edgeStars)
                  edgeStars->push_back(vector<SimplexId>{cid});
              }
              else if(edgeStars){
                edgeStars->at(edgeMap->at(edgeIds)-edgeIntervals_[nodeId-1]-1).push_back(cid);
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          pair_int edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // the edge is in the current node
              if(edgeIds.first > vertexIntervals_[nodeId-1] && edgeIds.first <= vertexIntervals_[nodeId]){
                if(edgeMap->find(edgeIds) == edgeMap->end()){
                  edgeCount++;
                  (*edgeMap)[edgeIds] = edgeCount+edgeIntervals_[nodeId-1];
                  if(edgeStars)
                    edgeStars->push_back(vector<SimplexId>{cid});
                }
                else if(edgeStars){
                  edgeStars->at(edgeMap->at(edgeIds)-edgeIntervals_[nodeId-1]-1).push_back(cid);
                }
              }
            }
          }
        }

        if(internalEdgeList){
          internalEdgeList->clear();
          internalEdgeList->resize(edgeMap->size());
          int k = 0;
          for(auto iter = edgeMap->begin(); iter != edgeMap->end(); iter++){
            (*internalEdgeList)[k++] = iter->first;
          }
        }

        if(!internalEdgeMap){
          delete edgeMap;
        }
        return 0;
      }

      /** 
       * Build the external edge list in the node.
       */ 
      int buildExternalEdgeMap(ExpandedNode * const nodePtr) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        SimplexId verticesPerCell = cellArray_[0];
        boost::unordered_map<SimplexId, vector<pair_int>> edgeNodes;
        
        // loop through the external cell list
        for(size_t i = 0; i < externalCells_[nodePtr->nid].size(); i++){
          pair_int edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * externalCells_[nodePtr->nid][i];

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // check if the edge is an external edge
              if(edgeIds.first <= vertexIntervals_[nodePtr->nid-1] && edgeIds.second > vertexIntervals_[nodePtr->nid-1] && edgeIds.second <= vertexIntervals_[nodePtr->nid]){
                SimplexId nid = findNodeIndex(edgeIds.first, VERTEX_ID);
                edgeNodes[nid].push_back(edgeIds);
              }
            }
          }
        }

        boost::unordered_map<SimplexId, vector<pair_int>>::iterator iter;
        for(iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++){
          ExpandedNode *exnode = searchCache(iter->first, nodePtr->nid);
          if(!exnode){
            boost::unordered_map<pair_int, SimplexId> localInternalEdgeMap;
            buildInternalEdgeList(iter->first, nullptr, &localInternalEdgeMap, nullptr);
            for(pair_int edgePair : iter->second){
              (*(nodePtr->externalEdgeMap_))[edgePair] = localInternalEdgeMap.at(edgePair);
            }
          }
          else{
            if(exnode->internalEdgeMap_ == nullptr){
              exnode->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
              buildInternalEdgeList(iter->first, nullptr, exnode->internalEdgeMap_, nullptr);
            }
            for(pair_int edgePair : iter->second){
              (*(nodePtr->externalEdgeMap_))[edgePair] = exnode->internalEdgeMap_->at(edgePair);
            }
          }
        }

        return 0;
      }

      /**
       * Build the internal triangle list in the node.
       */
      int buildInternalTriangleList(SimplexId nodeId, 
        vector<vector<SimplexId>> * const internalTriangleList,
        boost::unordered_map<vector<SimplexId>, SimplexId> * const internalTriangleMap,
        vector<vector<SimplexId>> * const triangleStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId triangleCount = 0, verticesPerCell = cellArray_[0];
        SimplexId localVertexNum = vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1];
        boost::unordered_map<vector<SimplexId>, SimplexId> *triangleMap;
        set<vector<SimplexId>> externalTriangleSet;

        if(triangleStars){
          triangleStars->clear();
          triangleStars->reserve(9*localVertexNum);
        }
        if(internalTriangleMap){
          internalTriangleMap->clear();
          triangleMap = internalTriangleMap;
        }else{
          triangleMap = new boost::unordered_map<vector<SimplexId>, SimplexId>();
        }

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            // the triangle does not belong to the current node
            if(triangleIds[0] > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell-1; k++){
              for(SimplexId l = k+1; l < verticesPerCell; l++){
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];
                
                if(triangleMap->find(triangleIds) == triangleMap->end()){
                  triangleCount++;
                  (*triangleMap)[triangleIds] = triangleCount + triangleIntervals_[nodeId-1];
                  if(triangleStars){
                    triangleStars->push_back(vector<SimplexId>{cid});
                  }
                }
                else if(triangleStars){
                  triangleStars->at(triangleMap->at(triangleIds)-triangleIntervals_[nodeId-1]-1).push_back(cid);
                }
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            if(triangleIds[0] > vertexIntervals_[nodeId-1] && triangleIds[0] <= vertexIntervals_[nodeId]){
              for(SimplexId k = j+1; k < verticesPerCell-1; k++){
                for(SimplexId l = k+1; l < verticesPerCell; l++){
                  triangleIds[1] = cellArray_[cellId + k + 1];
                  triangleIds[2] = cellArray_[cellId + l + 1];

                  if(triangleMap->find(triangleIds) == triangleMap->end()){
                      triangleCount++;
                      (*triangleMap)[triangleIds] = triangleCount + triangleIntervals_[nodeId-1];
                      if(internalTriangleList){
                        internalTriangleList->push_back(triangleIds);
                      }
                      if(triangleStars){
                        triangleStars->push_back(vector<SimplexId>{cid});
                      }
                  }
                  else if(triangleStars){
                    triangleStars->at(triangleMap->at(triangleIds)-triangleIntervals_[nodeId-1]-1).push_back(cid);
                  }
                }
              }
            }
          }
        }

        if(internalTriangleList){
          internalTriangleList->clear();
          internalTriangleList->resize(triangleMap->size());
          int k = 0;
          for(auto iter = triangleMap->begin(); iter != triangleMap->end(); iter++){
            (*internalTriangleList)[k++] = iter->first;
          }
        }

        if(!internalTriangleMap){
          delete triangleMap;
        }

        return 0;
      }

      /**
       * Build the external triangle list in the node.
       */
      int buildExternalTriangleMap(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        SimplexId verticesPerCell = cellArray_[0];
        boost::unordered_map<SimplexId, vector<vector<SimplexId>>> triangleNodes;

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            if(triangleIds[0] <= vertexIntervals_[nodePtr->nid-1]){
              for(SimplexId k = j+1; k < verticesPerCell-1; k++){
                for(SimplexId l = k+1; l < verticesPerCell; l++){
                  triangleIds[1] = cellArray_[cellId + k + 1];
                  triangleIds[2] = cellArray_[cellId + l + 1];

                  if(triangleIds[1] > vertexIntervals_[nodePtr->nid-1] && triangleIds[1] <= vertexIntervals_[nodePtr->nid]){
                    triangleNodes[findNodeIndex(triangleIds[0], VERTEX_ID)].push_back(triangleIds);
                  }
                  else if(triangleIds[2] > vertexIntervals_[nodePtr->nid-1] && triangleIds[2] <= vertexIntervals_[nodePtr->nid]){
                    triangleNodes[findNodeIndex(triangleIds[0], VERTEX_ID)].push_back(triangleIds);
                  }
                }
              }
            }
          }
        }

        ExpandedNode *exnode;
        boost::unordered_map<SimplexId, vector<vector<SimplexId>>>::iterator iter;
        for(iter = triangleNodes.begin(); iter != triangleNodes.end(); iter++){
          exnode = searchCache(iter->first, nodePtr->nid);

          if(!exnode){
            boost::unordered_map<vector<SimplexId>, SimplexId> localInternalTriangleMap;
            buildInternalTriangleList(iter->first, nullptr, &localInternalTriangleMap, nullptr);
            for(vector<SimplexId> triangleVec : iter->second){
              (*(nodePtr->externalTriangleMap_))[triangleVec] = localInternalTriangleMap.at(triangleVec);
            }
          }
          else{
            if(exnode->internalTriangleMap_ == nullptr){
              exnode->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
              buildInternalTriangleList(iter->first, nullptr, exnode->internalTriangleMap_, nullptr);
            }
            for(vector<SimplexId> triangleVec : iter->second){
              (*(nodePtr->externalTriangleMap_))[triangleVec] = exnode->internalTriangleMap_->at(triangleVec);
            }
          }
        }

        return 0;
      }

      /** 
       * Get the number of internal edges in the node.
       */ 
      SimplexId countInternalEdges(SimplexId nodeId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId edgeCount = 0, verticesPerCell = cellArray_[0];
        boost::unordered_set<pair_int, boost::hash<pair_int>> edgeSet;

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          pair_int edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            edgeIds.first = cellArray_[cellId + j + 1];
            // the edge does not belong to the current node
            if(edgeIds.first > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // not found in the edge map - assign new edge id
              if(edgeSet.find(edgeIds) == edgeSet.end()){
                edgeCount++;
                edgeSet.insert(edgeIds);
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          pair_int edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // the edge is in the current node
              if(edgeIds.first > vertexIntervals_[nodeId-1] && edgeIds.first <= vertexIntervals_[nodeId]){
                if(edgeSet.find(edgeIds) == edgeSet.end()){
                  edgeCount++;
                  edgeSet.insert(edgeIds);
                }
              }
            }
          }
        }

        return edgeCount;
      }

      /**
       * Get the number of internal triangles in the node.
       */
      int countInternalTriangles(SimplexId nodeId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId triangleCount = 0, verticesPerCell = cellArray_[0];
        boost::unordered_set<vector<SimplexId>, boost::hash<vector<SimplexId>>> triangleSet;

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            // the triangle does not belong to the current node
            if(triangleIds[0] > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell-1; k++){
              for(SimplexId l = k+1; l < verticesPerCell; l++){
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];
                
                if(triangleSet.find(triangleIds) == triangleSet.end()){
                  triangleCount++;
                  triangleSet.insert(triangleIds);
                }
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            if(triangleIds[0] > vertexIntervals_[nodeId-1] && triangleIds[0] <= vertexIntervals_[nodeId]){
              for(SimplexId k = j+1; k < verticesPerCell-1; k++){
                for(SimplexId l = k+1; l < verticesPerCell; l++){
                  triangleIds[1] = cellArray_[cellId + k + 1];
                  triangleIds[2] = cellArray_[cellId + l + 1];

                  if(triangleSet.find(triangleIds) == triangleSet.end()){
                    triangleCount++;
                    triangleSet.insert(triangleIds);
                  }
                }
              }
            }
          }
        }

        return triangleCount;
      }

      /**
       * Get the cell edges for all cells in a given node. 
       */
      int getCellEdges(ExpandedNode * const nodePtr) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        SimplexId edgesPerCell = cellArray_[0]*(cellArray_[0]-1)/2;
        nodePtr->cellEdges_->clear();
        nodePtr->cellEdges_->resize(cellIntervals_[nodePtr->nid]-cellIntervals_[nodePtr->nid-1], vector<SimplexId>(edgesPerCell));
        boost::unordered_map<SimplexId, vector<vector<SimplexId>>> edgeNodes;

        if(nodePtr->internalEdgeMap_ == nullptr){
          nodePtr->internalEdgeList_ = new vector<pair_int>();
          nodePtr->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildInternalEdgeList(nodePtr->nid, nodePtr->internalEdgeList_, nodePtr->internalEdgeMap_, nullptr);
        }

        for(SimplexId i = cellIntervals_[nodePtr->nid-1]+1; i <= cellIntervals_[nodePtr->nid]; i++){
          SimplexId cellId = (cellArray_[0]+1)*i;
          // get the internal edge id from the map
          for(SimplexId k = 1; k < cellArray_[0]; k++){
            pair_int edgePair(cellArray_[cellId+1], cellArray_[cellId+k+1]);
            (*nodePtr->cellEdges_)[i-cellIntervals_[nodePtr->nid-1]-1][k-1] = nodePtr->internalEdgeMap_->at(edgePair);
          }
          for(SimplexId j = 1; j < cellArray_[0]-1; j++){
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              pair_int edgePair(cellArray_[cellId+j+1], cellArray_[cellId+k+1]);
              if(edgePair.first <= vertexIntervals_[nodePtr->nid]){
                (*nodePtr->cellEdges_)[i-cellIntervals_[nodePtr->nid-1]-1][j+k] = nodePtr->internalEdgeMap_->at(edgePair);
              }
              // group the external edges by node id
              else{
                vector<SimplexId> edgeTuple{i, j+k, edgePair.first, edgePair.second};
                edgeNodes[findNodeIndex(edgePair.first, VERTEX_ID)].push_back(edgeTuple);
              }
            }
          }
        }

        ExpandedNode *exnode;
        for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++){
          exnode = searchCache(iter->first, nodePtr->nid);
          if(!exnode){
            boost::unordered_map<pair_int, SimplexId> localInternalEdgeMap;
            buildInternalEdgeList(iter->first, nullptr, &localInternalEdgeMap, nullptr);
            for(vector<SimplexId> edgeTuple : iter->second){
              pair_int edgePair(edgeTuple[2], edgeTuple[3]);
              (*nodePtr->cellEdges_)[edgeTuple[0]-cellIntervals_[nodePtr->nid-1]-1][edgeTuple[1]] = localInternalEdgeMap.at(edgePair);
            }
          }
          else{
            if(exnode->internalEdgeMap_ == nullptr){
              exnode->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
              buildInternalEdgeList(iter->first, nullptr, exnode->internalEdgeMap_, nullptr);
            }
            for(vector<SimplexId> edgeTuple : iter->second){
              pair_int edgePair(edgeTuple[2], edgeTuple[3]);
              (*nodePtr->cellEdges_)[edgeTuple[0]-cellIntervals_[nodePtr->nid-1]-1][edgeTuple[1]] = exnode->internalEdgeMap_->at(edgePair);
            }
          }
        }

        return 0;
      }

      /**
       * Get the cell triangles for all cells in a given node. 
       */
      int getCellTriangles(ExpandedNode * const nodePtr) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        SimplexId trianglesPerCell = cellArray_[0] * (cellArray_[0]-1) * (cellArray_[0]-2) / 6;
        boost::unordered_map<SimplexId, vector<vector<SimplexId>>> triangleNodes;

        nodePtr->cellTriangles_->clear();
        nodePtr->cellTriangles_->resize(cellIntervals_[nodePtr->nid]-cellIntervals_[nodePtr->nid-1], vector<SimplexId>(trianglesPerCell));

        if(nodePtr->internalTriangleMap_ == nullptr){
          nodePtr->internalTriangleList_ = new vector<vector<SimplexId>>();
          nodePtr->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildInternalTriangleList(nodePtr->nid, nodePtr->internalTriangleList_, nodePtr->internalTriangleMap_, nullptr);
        }

        for(SimplexId i = cellIntervals_[nodePtr->nid-1]+1; i <= cellIntervals_[nodePtr->nid]; i++){
          SimplexId cellId = (cellArray_[0]+1)*i;
          vector<SimplexId> triangleVec(3);
          // get the internal triangle from the map
          triangleVec[0] = cellArray_[cellId+1]; 
          for(SimplexId k = 1; k < cellArray_[0]-1; k++){
            triangleVec[1] = cellArray_[cellId+1+k];
            for(SimplexId l = k+1; l < cellArray_[0]; l++){
              triangleVec[2] = cellArray_[cellId+1+l];
              (*(nodePtr->cellTriangles_))[i-cellIntervals_[nodePtr->nid-1]-1][k+l-3] = nodePtr->internalTriangleMap_->at(triangleVec);
            }
          }
          // group the external triangles by node id
          triangleVec[0] = cellArray_[cellId+2];
          triangleVec[1] = cellArray_[cellId+3];
          triangleVec[2] = cellArray_[cellId+4];
          if(triangleVec[0] <= vertexIntervals_[nodePtr->nid]){
            (*(nodePtr->cellTriangles_))[i-cellIntervals_[nodePtr->nid-1]-1].back() = nodePtr->internalTriangleMap_->at(triangleVec);
          }
          else{
            vector<SimplexId> triangleTuple{i, triangleVec[0], triangleVec[1], triangleVec[2]};
            triangleNodes[findNodeIndex(triangleVec[0], VERTEX_ID)].push_back(triangleTuple);
          }
        }

        ExpandedNode *exnode;
        for(auto iter = triangleNodes.begin(); iter != triangleNodes.end(); iter++){
          exnode = searchCache(iter->first, nodePtr->nid);
          if(!exnode){
            boost::unordered_map<vector<SimplexId>, SimplexId> localInternalTriangleMap;
            buildInternalTriangleList(iter->first, nullptr, &localInternalTriangleMap, nullptr);
            for(vector<SimplexId> triangleVec : iter->second){
              vector<SimplexId> triangle(triangleVec.begin()+1, triangleVec.end());
              (*(nodePtr->cellTriangles_))[triangleVec[0]-cellIntervals_[nodePtr->nid-1]-1].back() = localInternalTriangleMap.at(triangle);
            }
          }
          else{
            if(exnode->internalTriangleMap_ == nullptr){
              exnode->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
              buildInternalTriangleList(iter->first, nullptr, exnode->internalTriangleMap_, nullptr);
            }
            for(vector<SimplexId> triangleVec : iter->second){
              vector<SimplexId> triangle(triangleVec.begin()+1, triangleVec.end());
              (*(nodePtr->cellTriangles_))[triangleVec[0]-cellIntervals_[nodePtr->nid-1]-1].back() = exnode->internalTriangleMap_->at(triangle);
            }
          }
        }

        return 0;
      }

      /**
       * Get the edge triangles for all the edges in a given node.
       */
      int getEdgeTriangles(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
              return -1;
        #endif

        nodePtr->edgeTriangles_->clear();
        nodePtr->edgeTriangles_->resize(edgeIntervals_[nodePtr->nid] - edgeIntervals_[nodePtr->nid-1]);

        if(nodePtr->internalEdgeMap_ == nullptr){
          nodePtr->internalEdgeList_ = new vector<pair_int>();
          nodePtr->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildInternalEdgeList(nodePtr->nid, nodePtr->internalEdgeList_, nodePtr->internalEdgeMap_, nullptr);
        }
        if(nodePtr->internalTriangleList_ == nullptr){
          nodePtr->internalTriangleList_ = new vector<vector<SimplexId>>();
          nodePtr->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildInternalTriangleList(nodePtr->nid, nodePtr->internalTriangleList_, nodePtr->internalTriangleMap_, nullptr);
        }
        if(nodePtr->externalTriangleMap_ == nullptr){
          nodePtr->externalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildExternalTriangleMap(nodePtr);
        }

        // for internal triangles
        for(SimplexId tid = 0; tid < (SimplexId) nodePtr->internalTriangleList_->size(); tid++){
          pair_int edge1 = pair_int((*(nodePtr->internalTriangleList_))[tid][0], (*(nodePtr->internalTriangleList_))[tid][1]);
          pair_int edge2 = pair_int((*(nodePtr->internalTriangleList_))[tid][0], (*(nodePtr->internalTriangleList_))[tid][2]);

          (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge1)-edgeIntervals_[nodePtr->nid-1]-1].push_back(tid + triangleIntervals_[nodePtr->nid-1] + 1);

          (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge2)-edgeIntervals_[nodePtr->nid-1]-1].push_back(tid + triangleIntervals_[nodePtr->nid-1] + 1);

          if((*(nodePtr->internalTriangleList_))[tid][1] <= vertexIntervals_[nodePtr->nid]){
            edge1 = pair_int((*(nodePtr->internalTriangleList_))[tid][1], (*(nodePtr->internalTriangleList_))[tid][2]);
            (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge1)-edgeIntervals_[nodePtr->nid-1]-1].push_back(tid + triangleIntervals_[nodePtr->nid-1] + 1);
          }
        }

        // for external triangles
        boost::unordered_map<vector<SimplexId>, SimplexId>::iterator iter;
        // loop through each edge of the cell
        for(iter = nodePtr->externalTriangleMap_->begin(); iter != nodePtr->externalTriangleMap_->end(); iter++){
          pair_int edge = pair<SimplexId, SimplexId >(iter->first.at(1), iter->first.at(2));
          if(edge.first > vertexIntervals_[nodePtr->nid-1] && edge.first <= vertexIntervals_[nodePtr->nid]){
            (*nodePtr->edgeTriangles_)[nodePtr->internalEdgeMap_->at(edge)-edgeIntervals_[nodePtr->nid-1]-1].push_back(iter->second);
          }
        }

        return 0;
      }

      /**
       * Get the triangle edges for all the triangles in a given node.
       */
      int getTriangleEdges(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->triangleEdges_->clear();
        nodePtr->triangleEdges_->resize(triangleIntervals_[nodePtr->nid]-triangleIntervals_[nodePtr->nid-1], vector<SimplexId>(3));
        boost::unordered_map<SimplexId, vector<vector<SimplexId>>> edgeNodes;

        if(nodePtr->internalEdgeMap_ == nullptr){
          nodePtr->internalEdgeList_ = new vector<pair_int>();
          nodePtr->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildInternalEdgeList(nodePtr->nid, nodePtr->internalEdgeList_, nodePtr->internalEdgeMap_, nullptr);
        }
        if(nodePtr->internalTriangleList_ == nullptr){
          nodePtr->internalTriangleList_ = new vector<vector<SimplexId>>();
          nodePtr->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildInternalTriangleList(nodePtr->nid, nodePtr->internalTriangleList_, nodePtr->internalTriangleMap_, nullptr);
        }

        for(SimplexId tid = 0; tid < (SimplexId) nodePtr->internalTriangleList_->size(); tid++){
          // since the first vertex of the triangle is in the node ...
          pair_int edgePair((*(nodePtr->internalTriangleList_))[tid][0], (*(nodePtr->internalTriangleList_))[tid][1]);
          (*nodePtr->triangleEdges_)[tid][0] = nodePtr->internalEdgeMap_->at(edgePair);
          edgePair.second = (*(nodePtr->internalTriangleList_))[tid][2];
          (*nodePtr->triangleEdges_)[tid][1] = nodePtr->internalEdgeMap_->at(edgePair);
          edgePair.first = (*(nodePtr->internalTriangleList_))[tid][1];
          if(edgePair.first > vertexIntervals_[nodePtr->nid-1] && edgePair.first <= vertexIntervals_[nodePtr->nid]){
            (*nodePtr->triangleEdges_)[tid][2] = nodePtr->internalEdgeMap_->at(edgePair);
          }else{
            vector<SimplexId> edgeTuple{tid, edgePair.first, edgePair.second};
            SimplexId nid = findNodeIndex(edgePair.first, VERTEX_ID);
            edgeNodes[nid].push_back(edgeTuple);
          }
        }

        ExpandedNode *exnode;
        for(auto iter = edgeNodes.begin(); iter != edgeNodes.end(); iter++){
          exnode = searchCache(iter->first, nodePtr->nid);
          if(!exnode){
            boost::unordered_map<pair_int, SimplexId> localInternalEdgeMap;
            buildInternalEdgeList(iter->first, nullptr, &localInternalEdgeMap, nullptr);
            for(vector<SimplexId> edgeTuple : iter->second){
              pair_int edgePair(edgeTuple[1], edgeTuple[2]);
              (*nodePtr->triangleEdges_)[edgeTuple[0]][2] = localInternalEdgeMap.at(edgePair);
            }
          }
          else{
            if(exnode->internalEdgeMap_ == nullptr){
              exnode->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
              buildInternalEdgeList(iter->first, nullptr, exnode->internalEdgeMap_, nullptr);
            }
            for(vector<SimplexId> edgeTuple : iter->second){
              pair_int edgePair(edgeTuple[1], edgeTuple[2]);
              (*nodePtr->triangleEdges_)[edgeTuple[0]][2] = exnode->internalEdgeMap_->at(edgePair);
            }
          }
        }

        return 0;
      }

      /**
       * Get the vertex edges for all the vertices in a given node.
       */
      int getVertexEdges(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->vertexEdges_->clear();
        nodePtr->vertexEdges_->resize(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1]);

        if(nodePtr->internalEdgeList_ == nullptr){
          nodePtr->internalEdgeList_ = new vector<pair_int>();
          nodePtr->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildInternalEdgeList(nodePtr->nid, nodePtr->internalEdgeList_, nodePtr->internalEdgeMap_, nullptr);
        }
        if(nodePtr->externalEdgeMap_ == nullptr){
          nodePtr->externalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildExternalEdgeMap(nodePtr);
        }

        for(SimplexId i = 0; i < (SimplexId) nodePtr->internalEdgeList_->size(); i++){
          (*(nodePtr->vertexEdges_))[nodePtr->internalEdgeList_->at(i).first-vertexIntervals_[nodePtr->nid-1]-1].push_back(edgeIntervals_[nodePtr->nid-1]+i+1);
          // the second vertex id of the edge must be greater than the first one
          if(nodePtr->internalEdgeList_->at(i).second <= vertexIntervals_[nodePtr->nid]){
            (*(nodePtr->vertexEdges_))[nodePtr->internalEdgeList_->at(i).second-vertexIntervals_[nodePtr->nid-1]-1].push_back(edgeIntervals_[nodePtr->nid-1]+i+1);
          }
        }

        boost::unordered_map<pair_int, SimplexId>::iterator iter;
        for(iter = nodePtr->externalEdgeMap_->begin(); iter != nodePtr->externalEdgeMap_->end(); iter++){
          (*(nodePtr->vertexEdges_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second);
        }

        return 0;
      }

      /**
       * Get the vertex neighbors for all the vertices in a given node.
       */
      int getVertexNeighbors(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->vertexNeighbors_->clear();
        nodePtr->vertexNeighbors_->resize(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1]);

        if(nodePtr->internalEdgeList_ == nullptr){
          nodePtr->internalEdgeList_ = new vector<pair_int>();
          nodePtr->internalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildInternalEdgeList(nodePtr->nid, nodePtr->internalEdgeList_, nodePtr->internalEdgeMap_, nullptr);
        }
        if(nodePtr->externalEdgeMap_ == nullptr){
          nodePtr->externalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildExternalEdgeMap(nodePtr);
        }

        for(SimplexId i = 0; i < (SimplexId) nodePtr->internalEdgeList_->size(); i++){
          (*(nodePtr->vertexNeighbors_))[nodePtr->internalEdgeList_->at(i).first-vertexIntervals_[nodePtr->nid-1]-1].push_back(nodePtr->internalEdgeList_->at(i).second);
          if(nodePtr->internalEdgeList_->at(i).second <= vertexIntervals_[nodePtr->nid])
            (*(nodePtr->vertexNeighbors_))[nodePtr->internalEdgeList_->at(i).second-vertexIntervals_[nodePtr->nid-1]-1].push_back(i+edgeIntervals_[nodePtr->nid]+1);
        }

        boost::unordered_map<pair_int, SimplexId>::iterator iter;
        for(iter = nodePtr->externalEdgeMap_->begin(); iter != nodePtr->externalEdgeMap_->end(); iter++){
          (*(nodePtr->vertexNeighbors_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second);
        }

        return 0;
      }

      /** 
       * Get the vertex stars for all the vertices in a given node. 
       * The function is similar as getVertexCells().
       */ 
      int getVertexStars(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->vertexStars_->clear();
        nodePtr->vertexStars_->resize(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1]);

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
          SimplexId cellId = (cellArray_[0]+1) * cid;
          for(SimplexId j = 0; j < cellArray_[0]; j++){
            // see if it is in the current node
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+j+1] <= vertexIntervals_[nodePtr->nid])
              (*(nodePtr->vertexStars_))[cellArray_[cellId+j+1]-vertexIntervals_[nodePtr->nid-1]-1].push_back(cid);
          }
        }

        // and also external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]){
          SimplexId cellId = (cellArray_[0]+1) * cid;
          for(SimplexId j = 0; j < cellArray_[0]; j++){
            // see if it is in the current node
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+j+1] <= vertexIntervals_[nodePtr->nid])
              (*(nodePtr->vertexStars_))[cellArray_[cellId+j+1]-vertexIntervals_[nodePtr->nid-1]-1].push_back(cid);
          }
        }

        return 0;
      }

      /** 
       * Get the vertex triangles for all the vertices in a given node. 
       */ 
      int getVertexTriangles(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->vertexTriangles_->clear();
        nodePtr->vertexTriangles_->resize(vertexIntervals_[nodePtr->nid] - vertexIntervals_[nodePtr->nid-1]);

        if(nodePtr->internalTriangleList_ == nullptr){
          nodePtr->internalTriangleList_ = new vector<vector<SimplexId>>();
          nodePtr->internalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildInternalTriangleList(nodePtr->nid, nodePtr->internalTriangleList_, nodePtr->internalTriangleMap_, nullptr);
        }

        if(nodePtr->externalTriangleMap_ == nullptr){
          nodePtr->externalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildExternalTriangleMap(nodePtr);
        }

        SimplexId triangleId = triangleIntervals_[nodePtr->nid-1]+1;
        for(SimplexId i = 0; i < (SimplexId) nodePtr->internalTriangleList_->size(); i++){
          for(SimplexId j = 0; j < 3; j++){
            if(nodePtr->internalTriangleList_->at(i)[j] > vertexIntervals_[nodePtr->nid-1] && nodePtr->internalTriangleList_->at(i)[j] <= vertexIntervals_[nodePtr->nid])
              (*(nodePtr->vertexTriangles_))[nodePtr->internalTriangleList_->at(i)[j]-vertexIntervals_[nodePtr->nid-1]-1].push_back(triangleId);
          }
          triangleId++;
        }

        boost::unordered_map<vector<SimplexId>, SimplexId>::iterator iter;
        for(iter = nodePtr->externalTriangleMap_->begin(); iter != nodePtr->externalTriangleMap_->end(); iter++){
          for(SimplexId j = 0; j < 3; j++){
            if(iter->first.at(j) > vertexIntervals_[nodePtr->nid-1] && iter->first.at(j) <= vertexIntervals_[nodePtr->nid])
              (*(nodePtr->vertexTriangles_))[iter->first.at(j)-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second);
          }
        }
        
        return 0;
      }

      
      /**
       * Protected class variables.
       */ 
      bool                doublePrecision_;
      SimplexId           cellNumber_, vertexNumber_, nodeNumber_;
      const void          *pointSet_;
      const int           *vertexIndices_;
      vector<SimplexId>   vertexIntervals_;
      vector<SimplexId>   edgeIntervals_;
      vector<SimplexId>   triangleIntervals_;
      vector<SimplexId>   cellIntervals_;
      const LongSimplexId *cellArray_;
      vector<vector<SimplexId>> externalCells_;

      // LRU cache
      size_t cacheSize_;
      mutable list<ExpandedNode*> cache_;
      mutable boost::unordered_map<SimplexId, list<ExpandedNode*>::iterator> cacheMap_;

      friend class TestStellar;
  };
}
