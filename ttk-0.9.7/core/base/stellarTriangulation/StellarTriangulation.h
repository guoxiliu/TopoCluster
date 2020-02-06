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
#include                  <set>
#include                  <map>
#include                  <list>
#include                  <unordered_map>
#include                  <AbstractTriangulation.h>

#define VERTEX_ID 0
#define EDGE_ID 1
#define TRIANGLE_ID 2
#define CELL_ID 3

using namespace std;
namespace ttk{

  class ExpandedNode
  {
  private:
    /* components */
    SimplexId nid;
    vector<pair<SimplexId, SimplexId>> *internalEdgeList_;
    vector<vector<SimplexId>> *internalTriangleList_;
    vector<pair<SimplexId, SimplexId>> *externalEdgeList_;
    vector<vector<SimplexId>> *externalTriangleList_;
    map<pair<SimplexId, SimplexId>, SimplexId> *internalEdgeMap_;
    map<vector<SimplexId>, SimplexId> *internalTriangleMap_;
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
      externalEdgeList_ = nullptr;
      externalTriangleList_ = nullptr;
      internalEdgeMap_ = nullptr;
      internalTriangleMap_ = nullptr;
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
      delete externalEdgeList_;
      delete externalTriangleList_;
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
      }

      int setInputPoints(const SimplexId &pointNumber, const void *pointSet, 
        const int *indexArray, const bool &doublePrecision = false){

        if(vertexNumber_)
            clear();

        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
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
        cacheSize_ = 100;

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
        
        if(!exnode->cellEdges_){
          exnode->cellEdges_ = new vector<vector<SimplexId>>();
          getCellEdges(nid, exnode->cellEdges_);
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
          if(!exnode->cellEdges_){
            exnode->cellEdges_ = new vector<vector<SimplexId>>();
            getCellEdges(nid, exnode->cellEdges_);
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
        if(!exnode->cellTriangles_){
          exnode->cellTriangles_ = new vector<vector<SimplexId>>();
          getCellTriangles(nid, exnode->cellTriangles_);
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
          if(!exnode->cellTriangles_){
            exnode->cellTriangles_ = new vector<vector<SimplexId>>();
            getCellTriangles(nid, exnode->cellTriangles_);
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
      
      const vector<pair<SimplexId, SimplexId> > *getEdges(){
        edgeList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->internalEdgeList_){
            exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
            exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
            buildEdgeList(nid, exnode->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
          }
          edgeList_.insert(edgeList_.end(), exnode->internalEdgeList_->begin(), exnode->internalEdgeList_->end());
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
        if(!exnode->edgeStars_){
          exnode->edgeStars_ = new vector<vector<SimplexId>>();
          buildEdgeList(nid, nullptr, nullptr, nullptr, exnode->edgeStars_);
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
        if(!exnode->edgeStars_){
          exnode->edgeStars_ = new vector<vector<SimplexId>>();
          buildEdgeList(nid, nullptr, nullptr, nullptr,  exnode->edgeStars_);
        }
        return (*(exnode->edgeStars_))[localEdgeId].size();
      }
      
      const vector<vector<SimplexId> > *getEdgeStars(){
        edgeStarList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->edgeStars_){
            exnode->edgeStars_ = new vector<vector<SimplexId>>();
            buildEdgeList(nid, nullptr, nullptr, nullptr, exnode->edgeStars_);
          }
          edgeStarList_.insert(edgeStarList_.end(), exnode->edgeStars_->begin(), exnode->edgeStars_->end());
        }
        return &edgeStarList_;
      }
     
      int getEdgeTriangle(const SimplexId &edgeId, 
        const int &localTriangleId, SimplexId &triangleId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.size()))
          return -1;
        if(localTriangleId < 0)
          return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(!exnode->edgeTriangles_){
          exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
          getEdgeTriangles(nid, exnode->edgeTriangles_);
        }
        if(localTriangleId >= (SimplexId) (*(exnode->edgeTriangles_))[localEdgeId].size())
          return -2;
        triangleId = (*(exnode->edgeTriangles_))[localEdgeId][localTriangleId];
        return 0;
      }
      
      SimplexId getEdgeTriangleNumber(const SimplexId &edgeId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((edgeId < 0)||(edgeId > (SimplexId) edgeIntervals_.size()))
          return -1;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        ExpandedNode *exnode = searchCache(nid);
        if(!exnode->edgeTriangles_){
          exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
          getEdgeTriangles(nid, exnode->edgeTriangles_);
          
        }
        return (*(exnode->edgeTriangles_))[edgeId-edgeIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getEdgeTriangles(){
        edgeTriangleList_.reserve(edgeIntervals_.size()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->edgeTriangles_){
            exnode->edgeTriangles_ = new vector<vector<SimplexId>>();
            getEdgeTriangles(nid, exnode->edgeTriangles_);

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
        if(!exnode->internalEdgeList_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nid, exnode->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
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
          if(!exnode->internalTriangleList_){
            exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
            exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
            buildTriangleList(nid, exnode->internalTriangleList_, nullptr, exnode->internalTriangleMap_, nullptr);
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
        if(!exnode->triangleEdges_){
          exnode->triangleEdges_ = new vector<vector<SimplexId>>();
          getTriangleEdges(nid, exnode->triangleEdges_);
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
        if(!exnode->triangleEdges_){
          exnode->triangleEdges_ = new vector<vector<SimplexId>>();
          getTriangleEdges(nid, exnode->triangleEdges_);
        }
        return (*(exnode->triangleEdges_))[triangleId-triangleIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleEdges(){
        triangleEdgeList_.reserve(triangleIntervals_.size()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->triangleEdges_){
            exnode->triangleEdges_ = new vector<vector<SimplexId>>();
            getTriangleEdges(nid, exnode->triangleEdges_);
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
        if(!exnode->triangleStars_){
          exnode->triangleStars_ = new vector<vector<SimplexId>>();
          buildTriangleList(nid, nullptr, nullptr, nullptr, exnode->triangleStars_);
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
        if(!exnode->triangleStars_){
          exnode->triangleStars_ = new vector<vector<SimplexId>>();
          buildTriangleList(nid, nullptr, nullptr, nullptr, exnode->triangleStars_);
        }
        
        return (*(exnode->triangleStars_))[triangleId-triangleIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleStars(){
        triangleStarList_.reserve(triangleIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->triangleStars_){
            exnode->triangleStars_ = new vector<vector<SimplexId>>();
            buildTriangleList(nid, nullptr, nullptr, nullptr, exnode->triangleStars_);
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
        if(!exnode->internalTriangleList_){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
          buildTriangleList(nid, exnode->internalTriangleList_, nullptr, exnode->internalTriangleMap_, nullptr);
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
        if(!exnode->vertexEdges_){
          exnode->vertexEdges_ = new vector<vector<SimplexId>>();
          getVertexEdges(nid, exnode->vertexEdges_);
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
        if(!exnode->vertexEdges_){
          exnode->vertexEdges_ = new vector<vector<SimplexId>>();
          getVertexEdges(nid, exnode->vertexEdges_);
        }
        return (*(exnode->vertexEdges_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexEdges(){
        vertexEdgeList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->vertexEdges_){
            exnode->vertexEdges_ = new vector<vector<SimplexId>>();
            getVertexEdges(nid, exnode->vertexEdges_);
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
        if(!exnode->vertexNeighbors_){
          exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
          getVertexNeighbors(nid, exnode->vertexNeighbors_);
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
        if(!exnode->vertexNeighbors_){
          exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
          getVertexNeighbors(nid, exnode->vertexNeighbors_);
        }
        return (*(exnode->vertexNeighbors_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexNeighbors(){
        vertexNeighborList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->vertexNeighbors_){
            exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
            getVertexNeighbors(nid, exnode->vertexNeighbors_);
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
        if(!exnode->vertexStars_){
          exnode->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(nid, exnode->vertexStars_);
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
        if(!exnode->vertexStars_){
          exnode->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(nid, exnode->vertexStars_);
        }
        return (*(exnode->vertexStars_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexStars(){
        vertexStarList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->vertexStars_){
            exnode->vertexStars_ = new vector<vector<SimplexId>>();
            getVertexStars(nid, exnode->vertexStars_);
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
        if(!exnode->vertexTriangles_){
          exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
          getVertexTriangles(nid, exnode->vertexTriangles_);
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
        if(!exnode->vertexTriangles_){
          exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
          getVertexTriangles(nid, exnode->vertexTriangles_);
        }
        return (*(exnode->vertexTriangles_))[vertexId-vertexIntervals_[nid-1]-1].size();
      }
      
      const vector<vector<SimplexId> > *getVertexTriangles(){
        vertexTriangleList_.reserve(vertexNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(!exnode->vertexTriangles_){
            exnode->vertexTriangles_ = new vector<vector<SimplexId>>();
            getVertexTriangles(nid, exnode->vertexTriangles_);
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
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            ExpandedNode *exnode = searchCache(nid);
            exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
            exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
            SimplexId edgeNum = buildEdgeList(nid, exnode ->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
            edgeIntervals_[nid] = edgeIntervals_[nid-1] + edgeNum;
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
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            ExpandedNode *exnode = searchCache(nid);
            exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
            exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
            SimplexId triangleNum = buildTriangleList(nid, exnode->internalTriangleList_, nullptr, exnode->internalTriangleMap_, nullptr);
            triangleIntervals_[nid] = triangleIntervals_[nid-1] + triangleNum;
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

      
    protected:

      int clear();

      /**
       * Find the corresponding node index given the id.
       */ 
      SimplexId findNodeIndex(SimplexId id, int idType) const{
        const vector<SimplexId> *intervals = nullptr;
        // determine which vector to search
        if(idType == VERTEX_ID){
          intervals = &vertexIntervals_;
        }else if(idType == EDGE_ID){
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
      ExpandedNode* searchCache(const SimplexId &nodeId, bool insertion=true) const{
        // cannot find the expanded node in the cache
        if(cacheMap_.find(nodeId) == cacheMap_.end()){
          if(!insertion){
            return new ExpandedNode(nodeId);
          }
          if(cache_.size() >= cacheSize_){
            cacheMap_.erase(cache_.back()->nid);
            delete cache_.back();
            cache_.pop_back();
          }
          cache_.push_front(new ExpandedNode(nodeId));
          cacheMap_[nodeId] = cache_.begin();
        }
        // find the expanded node in the cache
        else{
          cache_.splice(cache_.begin(), cache_, cacheMap_[nodeId]);
          cacheMap_[nodeId] = cache_.begin();
        }
        return cache_.front();
      }

      /** 
       * Build the whole edge list and return the number of edges in the node.
       */ 
      SimplexId buildEdgeList(SimplexId nodeId, 
        vector<pair<SimplexId, SimplexId>> * const internalEdgeList,
        vector<pair<SimplexId, SimplexId>> * const externalEdgeList,
        map<pair<SimplexId, SimplexId>, SimplexId> * const internalEdgeMap,
        vector<vector<SimplexId>> * const edgeStars) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId edgeCount = 0, verticesPerCell = cellArray_[0];
        map<pair<SimplexId, SimplexId>, SimplexId> *edgeMap = new map<pair<SimplexId, SimplexId>, SimplexId>();
        // loop through all the internal cells first
        if(internalEdgeList){
          internalEdgeList->clear();
          internalEdgeList->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(externalEdgeList){
          externalEdgeList->clear();
          externalEdgeList->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(edgeStars){
          edgeStars->clear();
          edgeStars->reserve(6*(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]));
        }
        if(internalEdgeMap){
          edgeMap = internalEdgeMap;
        }

        set<pair<SimplexId, SimplexId>> externalEdgeSet;

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          pair<SimplexId, SimplexId> edgeIds;
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
                if(internalEdgeList)
                  internalEdgeList->push_back(edgeIds);
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
          pair<SimplexId, SimplexId> edgeIds;
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
                  if(internalEdgeList)
                    internalEdgeList->push_back(edgeIds);
                  if(edgeStars)
                    edgeStars->push_back(vector<SimplexId>{cid});
                }
                else if(edgeStars){
                  edgeStars->at(edgeMap->at(edgeIds)-edgeIntervals_[nodeId-1]-1).push_back(cid);
                }
              }
              // the edge is a bridge between two nodes
              else if(externalEdgeList && edgeIds.second > vertexIntervals_[nodeId-1] && edgeIds.second <= vertexIntervals_[nodeId]){
                externalEdgeSet.insert(edgeIds);
              }
            }
          }
        }

        if(externalEdgeList)
          externalEdgeList->insert(externalEdgeList->end(), externalEdgeSet.begin(), externalEdgeSet.end());
        
        return edgeCount;
      }

      /**
       * Build the whole triangle list and return the number of triangles in the node.
       */
      int buildTriangleList(SimplexId nodeId, 
        vector<vector<SimplexId>> * const internalTriangleList,
        vector<vector<SimplexId>> * const externalTriangleList, 
        map<vector<SimplexId>,SimplexId> * const internalTriangleMap,
        vector<vector<SimplexId>> * const triangleStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId triangleCount = 0, verticesPerCell = cellArray_[0];
        SimplexId localVertexNum = vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1];
        map<vector<SimplexId>,SimplexId> *triangleMap = new map<vector<SimplexId>,SimplexId>();
        set<vector<SimplexId>> externalTriangleSet;

        if(internalTriangleList){
          internalTriangleList->clear();
          internalTriangleList->reserve(9*localVertexNum);
        }
        if(externalTriangleList){
          externalTriangleList->clear();
          externalTriangleList->reserve(9*localVertexNum);
        }
        if(triangleStars){
          triangleStars->clear();
          triangleStars->reserve(9*localVertexNum);
        }
        if(internalTriangleMap){
          triangleMap = internalTriangleMap;
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

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            for(SimplexId k = j+1; k < verticesPerCell-1; k++){
              for(SimplexId l = k+1; l < verticesPerCell; l++){
                triangleIds[0] = cellArray_[cellId + j + 1];
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];

                for(SimplexId v = 0; v < 3; v++){
                  if(triangleIds[v] > vertexIntervals_[nodeId-1] && triangleIds[v] <= vertexIntervals_[nodeId]){
                    if(triangleMap->find(triangleIds) == triangleMap->end()){
                      if(v == 0){
                        triangleCount++;
                        (*triangleMap)[triangleIds] = triangleCount + triangleIntervals_[nodeId-1];
                        if(internalTriangleList){
                          internalTriangleList->push_back(triangleIds);
                        }
                        if(triangleStars){
                          triangleStars->push_back(vector<SimplexId>{cid});
                        }
                      }
                      else if(externalTriangleList){
                        externalTriangleSet.insert(triangleIds);
                      }
                    }
                    else{
                      if(v == 0 && triangleStars){
                        triangleStars->at(triangleMap->at(triangleIds)-triangleIntervals_[nodeId-1]-1).push_back(cid);
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if(externalTriangleList){
          externalTriangleList->insert(externalTriangleList->end(), externalTriangleSet.begin(), externalTriangleSet.end());
        }

        return triangleCount;
      }

      /**
       * Get the cell edges for all cells in a given node. 
       */
      int getCellEdges(const SimplexId &nodeId, vector<vector<SimplexId>> * const cellEdges) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        cellEdges->clear();
        cellEdges->resize(cellIntervals_[nodeId]-cellIntervals_[nodeId-1]);

        ExpandedNode *exnode = searchCache(nodeId);
        if(!exnode->internalEdgeMap_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nodeId, exnode->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
        }

        for(SimplexId i = cellIntervals_[nodeId-1]+1; i <= cellIntervals_[nodeId]; i++){
          SimplexId cellId = (cellArray_[0]+1)*i;
          for(SimplexId j = 0; j < cellArray_[0]-1; j++){
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              pair<SimplexId, SimplexId> edgePair(cellArray_[cellId+j+1], cellArray_[cellId+k+1]);
              if(edgePair.first <= vertexIntervals_[nodeId]){
                cellEdges->at(i-cellIntervals_[nodeId-1]-1).push_back(exnode->internalEdgeMap_->at(edgePair));
              }
              else{
                cellEdges->at(i-cellIntervals_[nodeId-1]-1).push_back(getEdgeId(edgePair));
              }
            }
          }
        }
        return 0;
      }

      /**
       * Get the cell triangles for all cells in a given node. 
       */
      int getCellTriangles(const SimplexId &nodeId, vector<vector<SimplexId>> * const cellTriangles) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        cellTriangles->clear();
        cellTriangles->resize(cellIntervals_[nodeId]-cellIntervals_[nodeId-1]);

        for(SimplexId i = cellIntervals_[nodeId-1]+1; i <= cellIntervals_[nodeId]; i++){
          SimplexId cellId = (cellArray_[0]+1)*i;
          vector<SimplexId> triangleVec(3);
          for(SimplexId j = 0; j < cellArray_[0]-2; j++){
            triangleVec[0] = cellArray_[cellId+1+j]; 
            for(SimplexId k = j+1; k < cellArray_[0]-1; k++){
              triangleVec[1] = cellArray_[cellId+1+k];
              for(SimplexId l = k+1; l < cellArray_[0]; l++){
                triangleVec[2] = cellArray_[cellId+1+l];
                cellTriangles->at(i-cellIntervals_[nodeId-1]-1).push_back(getTriangleId(triangleVec));
              }
            }
          }
        }
        return 0;
      }

      /**
       * Get the edge id given a pair of vertex ids.
       * Note: make sure the passed pair is already sorted.
       */
      SimplexId getEdgeId(pair<SimplexId, SimplexId> &edgePair) const{
        SimplexId nid = findNodeIndex(edgePair.first, VERTEX_ID);
        ExpandedNode *exnode = searchCache(nid, false);
        if(!exnode->internalEdgeMap_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nid, exnode->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
        }

        if(exnode->internalEdgeMap_->find(edgePair) == exnode->internalEdgeMap_->end()){
          return -1;
        }

        return exnode->internalEdgeMap_->at(edgePair);
      }

      /**
       * Get the edge triangles for all the edges in a given node.
       */
      int getEdgeTriangles(const SimplexId& nodeId, vector<vector<SimplexId>> * const edgeTriangles) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
              return -1;
        #endif

        edgeTriangles->clear();
        edgeTriangles->resize(edgeIntervals_[nodeId] - edgeIntervals_[nodeId-1]);

        ExpandedNode *exnode = searchCache(nodeId);
        if(!exnode->internalEdgeMap_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nodeId, exnode->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
        }
        if(!exnode->internalTriangleList_){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
          exnode->externalTriangleList_ = new vector<vector<SimplexId>>();
          buildTriangleList(nodeId, exnode->internalTriangleList_, exnode->externalTriangleList_, exnode->internalTriangleMap_, nullptr);
        }
        if(!exnode->externalTriangleList_){
          exnode->externalTriangleList_ = new vector<vector<SimplexId>>();
          buildTriangleList(nodeId, exnode->internalTriangleList_, exnode->externalTriangleList_, exnode->internalTriangleMap_, nullptr);
        }


        // for internal triangles
        pair<SimplexId, SimplexId> edgeIds;
        for(SimplexId tid = 0; tid < (SimplexId) exnode->internalTriangleList_->size(); tid++){
          for(SimplexId j = 0; j < 2; j++){
            if((*(exnode->internalTriangleList_))[tid][j] > vertexIntervals_[nodeId]){
              break;
            }
            for(SimplexId k = j+1; k < 3; k++){
              edgeIds = pair<SimplexId,SimplexId>((*(exnode->internalTriangleList_))[tid][j], (*(exnode->internalTriangleList_))[tid][k]);
              (*edgeTriangles)[(*(exnode->internalEdgeMap_))[edgeIds]-edgeIntervals_[nodeId-1]-1].push_back(
                tid + triangleIntervals_[nodeId-1] + 1);
            }
          }
        }

        // for external triangles
        for(vector<SimplexId> triangle : (*(exnode->externalTriangleList_))){          // loop through each edge of the cell
          SimplexId triangleId = getTriangleId(triangle);
          for(SimplexId j = 0; j < 2; j++){
            for(SimplexId k = j+1; k < 3; k++){
              edgeIds.first = triangle[j];
              edgeIds.second = triangle[k];
              
              // the edge is in the current node
              if(edgeIds.first > vertexIntervals_[nodeId-1] && edgeIds.first <= vertexIntervals_[nodeId]){
                (*edgeTriangles)[(*(exnode->internalEdgeMap_))[edgeIds]-edgeIntervals_[nodeId-1]-1].push_back(triangleId);
              }
            }
          }
        }

        return 0;
      }

      /**
       * Get the edge id given a pair of vertex ids.
       * Note: make sure the passed vector is already sorted.
       */
      SimplexId getTriangleId(vector<SimplexId> &triangle) const{
        SimplexId nid = findNodeIndex(triangle.front(), VERTEX_ID);
        ExpandedNode *exnode = searchCache(nid, false);
        if(!exnode->internalTriangleMap_){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
          buildTriangleList(nid, exnode->internalTriangleList_, nullptr, exnode->internalTriangleMap_, nullptr);
        }

        if(exnode->internalTriangleMap_->find(triangle) == exnode->internalTriangleMap_->end())
          return -1;

        return exnode->internalTriangleMap_->at(triangle);
      }

      /**
       * Get the triangle edges for all the triangles in a given node.
       */
      int getTriangleEdges(const SimplexId &nodeId, vector<vector<SimplexId>> *triangleEdges) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        triangleEdges->clear();
        triangleEdges->resize(triangleIntervals_[nodeId]-triangleIntervals_[nodeId-1], vector<SimplexId>(3));

        ExpandedNode *exnode = searchCache(nodeId);
        if(!exnode->internalEdgeMap_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nodeId, exnode->internalEdgeList_, nullptr, exnode->internalEdgeMap_, nullptr);
        }
        if(!exnode->internalTriangleList_){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
          buildTriangleList(nodeId, exnode->internalTriangleList_, nullptr, exnode->internalTriangleMap_, nullptr);
        }

        for(SimplexId tid = 0; tid < (SimplexId) exnode->internalTriangleList_->size(); tid++){
          // since the first vertex of the triangle is in the node ...
          pair<SimplexId,SimplexId> edgePair(exnode->internalTriangleList_->at(tid)[0], exnode->internalTriangleList_->at(tid)[1]);
          (*triangleEdges)[tid][0] = exnode->internalEdgeMap_->at(edgePair);
          edgePair.second = exnode->internalTriangleList_->at(tid)[2];
          (*triangleEdges)[tid][1] = exnode->internalEdgeMap_->at(edgePair);
          edgePair.first = exnode->internalTriangleList_->at(tid)[1];
          if(edgePair.first > vertexIntervals_[nodeId-1] && edgePair.second <= vertexIntervals_[nodeId]){
            (*triangleEdges)[tid][2] = exnode->internalEdgeMap_->at(edgePair);
          }else{
            (*triangleEdges)[tid][2] = getEdgeId(edgePair);
          }
        }

        return 0;
      }

      /**
       * Get the vertex edges for all the vertices in a given node.
       */
      int getVertexEdges(const SimplexId &nodeId, vector<vector<SimplexId>> * const vertexEdges) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexEdges->clear();
        vertexEdges->resize(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);

        ExpandedNode *exnode = searchCache(nodeId);
        if(!exnode->internalEdgeList_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->externalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nodeId, exnode->internalEdgeList_, exnode->externalEdgeList_, exnode->internalEdgeMap_, nullptr);
        }
        if(!exnode->externalEdgeList_){
          exnode->externalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          buildEdgeList(nodeId, nullptr, exnode->externalEdgeList_, nullptr, nullptr);
        }

        for(SimplexId i = 0; i < (SimplexId) exnode->internalEdgeList_->size(); i++){
          (*vertexEdges)[exnode->internalEdgeList_->at(i).first-vertexIntervals_[nodeId-1]-1].push_back(edgeIntervals_[nodeId-1]+i+1);
          // the second vertex id of the edge must be greater than the first one
          if(exnode->internalEdgeList_->at(i).second <= vertexIntervals_[nodeId])
            (*vertexEdges)[exnode->internalEdgeList_->at(i).second-vertexIntervals_[nodeId-1]-1].push_back(edgeIntervals_[nodeId-1]+i+1);
        }

        for(SimplexId i = 0; i < (SimplexId) exnode->externalEdgeList_->size(); i++){
          (*vertexEdges)[exnode->externalEdgeList_->at(i).second-vertexIntervals_[nodeId-1]-1].push_back(getEdgeId(exnode->externalEdgeList_->at(i)));
        }

        return 0;
      }

      /**
       * Get the vertex neighbors for all the vertices in a given node.
       */
      int getVertexNeighbors(const SimplexId &nodeId, vector<vector<SimplexId>> * const vertexNeighbors) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexNeighbors->clear();
        vertexNeighbors->resize(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);

        ExpandedNode *exnode = searchCache(nodeId);
        if(!exnode->internalEdgeList_){
          exnode->internalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->externalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          exnode->internalEdgeMap_ = new map<pair<SimplexId, SimplexId>, SimplexId>();
          buildEdgeList(nodeId, exnode->internalEdgeList_, exnode->externalEdgeList_, exnode->internalEdgeMap_, nullptr);
        }
        if(!exnode->externalEdgeList_){
          exnode->externalEdgeList_ = new vector<pair<SimplexId, SimplexId>>();
          buildEdgeList(nodeId, nullptr, exnode->externalEdgeList_, nullptr, nullptr);
        }

        for(SimplexId i = 0; i < (SimplexId) exnode->internalEdgeList_->size(); i++){
          (*vertexNeighbors)[exnode->internalEdgeList_->at(i).first-vertexIntervals_[nodeId-1]-1].push_back(exnode->internalEdgeList_->at(i).second);
          if(exnode->internalEdgeList_->at(i).second <= vertexIntervals_[nodeId])
            (*vertexNeighbors)[exnode->internalEdgeList_->at(i).second-vertexIntervals_[nodeId-1]-1].push_back(i+edgeIntervals_[nodeId]+1);
        }

        for(SimplexId i = 0; i < (SimplexId) exnode->externalEdgeList_->size(); i++){
          (*vertexNeighbors)[exnode->externalEdgeList_->at(i).second-vertexIntervals_[nodeId-1]-1].push_back(exnode->externalEdgeList_->at(i).first);
        }

        return 0;
      }

      /** 
       * Get the vertex stars for all the vertices in a given node. 
       * The function is similar as getVertexCells().
       */ 
      int getVertexStars(const SimplexId &nodeId, vector<vector<SimplexId>> * const vertexStars) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexStars->clear();
        vertexStars->resize(vertexIntervals_[nodeId]-vertexIntervals_[nodeId-1]);
        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodeId-1]+1; cid <= cellIntervals_[nodeId]; cid++){
          SimplexId cellId = (cellArray_[0]+1) * cid;
          for(SimplexId j = 0; j < cellArray_[0]; j++){
            // see if it is in the current node
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodeId-1] && cellArray_[cellId+j+1] <= vertexIntervals_[nodeId])
              (*vertexStars)[cellArray_[cellId+j+1]-vertexIntervals_[nodeId-1]-1].push_back(cid);
          }
        }

        // and also external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          SimplexId cellId = (cellArray_[0]+1) * cid;
          for(SimplexId j = 0; j < cellArray_[0]; j++){
            // see if it is in the current node
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodeId-1] && cellArray_[cellId+j+1] <= vertexIntervals_[nodeId])
              (*vertexStars)[cellArray_[cellId+j+1]-vertexIntervals_[nodeId-1]-1].push_back(cid);
          }
        }

        return 0;
      }

      /** 
       * Get the vertex triangles for all the vertices in a given node. 
       */ 
      int getVertexTriangles(const SimplexId &nodeId, vector<vector<SimplexId>> * const vertexTriangles) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        vertexTriangles->clear();
        vertexTriangles->resize(vertexIntervals_[nodeId] - vertexIntervals_[nodeId-1]);

        ExpandedNode *exnode = searchCache(nodeId);
        if(!exnode->internalTriangleList_){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          exnode->internalTriangleMap_ = new map<vector<SimplexId>, SimplexId>();
          exnode->externalTriangleList_ = new vector<vector<SimplexId>>();
          buildTriangleList(nodeId, exnode->internalTriangleList_, exnode->externalTriangleList_, exnode->internalTriangleMap_, nullptr);
        }
        if(!exnode->externalTriangleList_){
          exnode->externalTriangleList_ = new vector<vector<SimplexId>>();
          buildTriangleList(nodeId, nullptr, exnode->externalTriangleList_, nullptr, nullptr);
        }

        SimplexId triangleId = triangleIntervals_[nodeId-1]+1;
        for(SimplexId i = 0; i < (SimplexId) exnode->internalTriangleList_->size(); i++){
          for(SimplexId j = 0; j < 3; j++){
            if(exnode->internalTriangleList_->at(i)[j] > vertexIntervals_[nodeId-1] && exnode->internalTriangleList_->at(i)[j] <= vertexIntervals_[nodeId])
              (*vertexTriangles)[exnode->internalTriangleList_->at(i)[j]-vertexIntervals_[nodeId-1]-1].push_back(triangleId);
          }
          triangleId++;
        }
        for(SimplexId i = 0; i < (SimplexId) exnode->externalTriangleList_->size(); i++){
          triangleId = getTriangleId(exnode->externalTriangleList_->at(i));
          for(SimplexId j = 0; j < 3; j++){
            if(exnode->externalTriangleList_->at(i)[j] > vertexIntervals_[nodeId-1] && exnode->externalTriangleList_->at(i)[j] <= vertexIntervals_[nodeId])
              (*vertexTriangles)[exnode->externalTriangleList_->at(i)[j]-vertexIntervals_[nodeId-1]-1].push_back(triangleId);
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
      vector<SimplexId>   vertexIntervals_;
      vector<SimplexId>   edgeIntervals_;
      vector<SimplexId>   triangleIntervals_;
      vector<SimplexId>   cellIntervals_;
      const LongSimplexId *cellArray_;
      vector<vector<SimplexId>> externalCells_;

      // LRU cache
      size_t cacheSize_;
      mutable list<ExpandedNode*> cache_;
      mutable unordered_map<SimplexId, list<ExpandedNode*>::iterator> cacheMap_;
  };
}
