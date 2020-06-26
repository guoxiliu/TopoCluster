/// \ingroup base
/// \class ttk::ExplicitTopoCluster 
/// \author Guoxi Liu <guoxil@g.clemson.edu>
/// \date Nov. 2019.
///
/// \brief TTK %ExplicitTopoCluster processing package.
///
/// %ExplicitTopoCluster is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkExplicitTopoCluster.cpp %for a usage example.

#pragma once

// base code includes
#include                  <algorithm>
#include                  <list>
#include                  <boost/unordered_map.hpp>
#include                  <boost/unordered_set.hpp>
#include                  <AbstractTriangulation.h>
#include                  <boost/functional/hash.hpp>

#define EDGE_ID 1
#define TRIANGLE_ID 2
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
    boost::unordered_map<pair_int, SimplexId> *externalEdgeMap_;
    boost::unordered_map<vector<SimplexId>, SimplexId> *externalTriangleMap_;
    /* boundary cells */
    vector<bool> *boundaryVertices_;
    vector<bool> *boundaryEdges_;
    vector<bool> *boundaryTriangles_;
    /* vertex relationships */
    vector<vector<SimplexId>> *vertexEdges_;
    vector<vector<SimplexId>> *vertexLinks_;
    vector<vector<SimplexId>> *vertexNeighbors_;
    vector<vector<SimplexId>> *vertexStars_;
    vector<vector<SimplexId>> *vertexTriangles_;
    /* edge relationships */
    // edgeVertex relation can be extracted from internal edge list
    vector<vector<SimplexId>> *edgeLinks_;
    vector<vector<SimplexId>> *edgeStars_;
    vector<vector<SimplexId>> *edgeTriangles_;
    /* triangle relationships */
    // triangleVertex relation can be extracted from internal triangle list
    vector<vector<SimplexId>> *triangleEdges_;
    vector<vector<SimplexId>> *triangleLinks_;
    vector<vector<SimplexId>> *triangleStars_;
    /* cell relationships */
    vector<vector<SimplexId>> *cellEdges_;
    vector<vector<SimplexId>> *cellNeighbors_;
    vector<vector<SimplexId>> *cellTriangles_;

  public:
    ExpandedNode(SimplexId id){
      /* components */
      nid = id;
      internalEdgeList_ = nullptr;
      internalTriangleList_ = nullptr;
      externalEdgeMap_ = nullptr;
      externalTriangleMap_ = nullptr;
      /* boundary cells */
      boundaryEdges_ = nullptr;
      boundaryVertices_ = nullptr;
      boundaryTriangles_ = nullptr;
      /* vertex relationships */
      vertexEdges_ = nullptr;
      vertexLinks_ = nullptr;
      vertexNeighbors_ = nullptr;
      vertexStars_ = nullptr;
      vertexTriangles_ = nullptr;
      /* edge relationships */
      edgeLinks_ = nullptr;
      edgeStars_ = nullptr;
      edgeTriangles_ = nullptr;
      /* triangle relationships */
      triangleLinks_ = nullptr;
      triangleEdges_ = nullptr;
      triangleStars_ = nullptr;
      /* cell relationships */
      cellEdges_ = nullptr;
      cellNeighbors_ = nullptr;
      cellTriangles_ = nullptr;
    }
    ~ExpandedNode(){
      delete internalEdgeList_;
      delete internalTriangleList_;
      delete externalEdgeMap_;
      delete externalTriangleMap_;
      delete boundaryEdges_;
      delete boundaryVertices_;
      delete boundaryTriangles_;
      delete vertexEdges_;
      delete vertexLinks_;
      delete vertexNeighbors_;
      delete vertexStars_;
      delete vertexTriangles_;
      delete edgeLinks_;
      delete edgeStars_;
      delete edgeTriangles_;
      delete triangleEdges_;
      delete triangleLinks_;
      delete triangleStars_;
      delete cellEdges_;
      delete cellNeighbors_;
      delete cellTriangles_;
    }

    friend class ExplicitTopoCluster;
  };

  
  class ExplicitTopoCluster : public AbstractTriangulation{

    public:

      ExplicitTopoCluster();

      ~ExplicitTopoCluster();

      /**
       * Set up vertices from the input.
       */ 
      int setInputPoints(const SimplexId &pointNumber, const void *pointSet, 
        const int *indexArray, const bool &doublePrecision = false){

        if(vertexNumber_)
            clear();

        vertexNumber_ = pointNumber;
        pointSet_ = pointSet;
        vertexIndices_ = indexArray;
        doublePrecision_ = doublePrecision;

        return 0;
      }

      /**
       * Set up cells from the input.
       */ 
      int setInputCells(const SimplexId &cellNumber,
        const LongSimplexId *cellArray){

        if(cellNumber_)
          clear();

        if((!cellArray)||(!cellNumber))
          return -1;
        
        cellNumber_ = cellNumber;
        cellArray_ = cellArray;

        vector<SimplexId> vertexMap(vertexNumber_);
        reorderVertices(vertexMap);
        reorderCells(vertexMap);

        return 0;
      }
      
      /**
       * Reorder the input vertices.
       */
      int reorderVertices(vector<SimplexId>& vertexMap){
        // get the number of nodes (the max value in the array)
        for(SimplexId vid = 0; vid < vertexNumber_; vid++){
          if(vertexIndices_[vid] > nodeNumber_){
            nodeNumber_ = vertexIndices_[vid];
          }
        }
        nodeNumber_++;  // since the index starts from 0
        vector<vector<SimplexId>> nodeVertices(nodeNumber_);
        for(SimplexId vid = 0; vid < vertexNumber_; vid++){
          nodeVertices[vertexIndices_[vid]].push_back(vid);
        }

        // update the vertex intervals
        vertexIntervals_.resize(nodeNumber_+1);
        vertexIntervals_[0] = -1;
        SimplexId vertexCount = 0;
        for(SimplexId nid = 0; nid < nodeNumber_; nid++){
          for(SimplexId vid : nodeVertices[nid]){
            vertexMap[vid] = vertexCount++;
          }
          vertexIntervals_[nid+1] = vertexCount - 1;
        }

        // rearange the vertex coordinate values
        if(doublePrecision_){
          double *newPointSet = new double[3*vertexNumber_];
          for(SimplexId vid = 0; vid < vertexNumber_; vid++){
            for(int j = 0; j < 3; j++){
              newPointSet[3*vertexMap[vid]+j] = ((double*) pointSet_)[3*vid+j];
            }
          }
          for(SimplexId vid = 0; vid < vertexNumber_; vid++){
            for(int j = 0; j < 3; j++){
              ((double*) pointSet_)[3*vid+j] = newPointSet[3*vid+j];
            }
          }
          delete[] newPointSet;
        }
        else{
          float *newPointSet = new float[3*vertexNumber_];
          for(SimplexId vid = 0; vid < vertexNumber_; vid++){
            for(int j = 0; j < 3; j++){
              newPointSet[3*vertexMap[vid]+j] = ((float*) pointSet_)[3*vid+j];
            }
          }
          for(SimplexId vid = 0; vid < vertexNumber_; vid++){
            for(int j = 0; j < 3; j++){
              ((float*) pointSet_)[3*vid+j] = newPointSet[3*vid+j];
            }
          }
          delete[] newPointSet;
        }

        // change the vertex indices
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          for(SimplexId vid = vertexIntervals_[nid-1]+1; vid <= vertexIntervals_[nid]; vid++){
            ((int*)vertexIndices_)[vid] = nid;
          }
        }

        return 0;
      }
      
      /**
       * Reorder the input cells.
       */
      int reorderCells(const vector<SimplexId>& vertexMap){
        // change the indices in cell array
        SimplexId cellCount = 0, verticesPerCell = cellArray_[0];
        vector<vector<SimplexId>> nodeCells(nodeNumber_+1);
        ttk::LongSimplexId *cellArr = ((ttk::LongSimplexId *)cellArray_);

        for(SimplexId cid = 0; cid < cellNumber_; cid++){
          SimplexId cellId = (verticesPerCell+1) * cid;
          for(int j = 1; j <= verticesPerCell; j++){
            cellArr[cellId+j] = vertexMap[cellArr[cellId+j]];
          }
          sort(cellArr+cellId+1, cellArr+cellId+1+verticesPerCell);
          nodeCells[vertexIndices_[cellArr[cellId+1]]].push_back(cid);
        }
        
        // rearange the cell array
        cellIntervals_.resize(nodeNumber_+1);
        externalCells_.resize(nodeNumber_+1);
        cellIntervals_[0] = -1;
        ttk::LongSimplexId *newCellArray = new ttk::LongSimplexId[(verticesPerCell+1)*cellNumber_];
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          for(SimplexId cid : nodeCells[nid]){
            SimplexId cellId = (verticesPerCell+1)*cid;
            SimplexId newCellId = (verticesPerCell+1)*cellCount;
            newCellArray[newCellId] = verticesPerCell;
            for(int j = 1; j <= verticesPerCell; j++){
              newCellArray[newCellId+j] = cellArray_[cellId+j];
              if(newCellArray[newCellId+j] > vertexIntervals_[nid]){
                SimplexId nodeNum = vertexIndices_[newCellArray[newCellId+j]];
                if(externalCells_[nodeNum].empty() || externalCells_[nodeNum].back() != cid){
                  externalCells_[nodeNum].push_back(cid);
                }
              }
            }
            cellCount++;
          }
          cellIntervals_[nid] = cellCount - 1;
        }

        // copy the new cell array back to original one
        for(SimplexId i = 0; i < (verticesPerCell+1)*cellNumber_; i++){
          ((ttk::LongSimplexId*)cellArray_)[i] = newCellArray[i];
        }
        delete[] newCellArray;

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

        SimplexId nid = vertexIndices_[cellArray_[(cellArray_[0]+1)*cellId+1]];
        SimplexId localCellId = cellId - cellIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->cellEdges_ == nullptr){
          exnode->cellEdges_ = new vector<vector<SimplexId>>();
          getCellEdges(exnode);
        }
        if(localEdgeId >= (SimplexId) (*(exnode->cellEdges_))[localCellId].size())
          return -2;
        
        edgeId = (*(exnode->cellEdges_))[localCellId][localEdgeId];
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
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
          if((localNeighborId < 0))
            return -2;
        #endif

        SimplexId nid = vertexIndices_[cellArray_[(cellArray_[0]+1)*cellId+1]];
        SimplexId localCellId = cellId - cellIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->cellNeighbors_ == nullptr){
          exnode->cellNeighbors_ = new vector<vector<SimplexId>>();
          getCellNeighbors(exnode);
        }
        if(localNeighborId >= (SimplexId) (*(exnode->cellNeighbors_))[localCellId].size())
          return -2;
        
        neighborId = (*(exnode->cellNeighbors_))[localCellId][localNeighborId];
        return 0;
      }
        
      SimplexId getCellNeighborNumber(const SimplexId &cellId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((cellId < 0)||(cellId >= cellNumber_))
            return -1;
        #endif

        SimplexId nid = vertexIndices_[cellArray_[(cellArray_[0]+1)*cellId+1]];
        SimplexId localCellId = cellId - cellIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->cellNeighbors_ == nullptr){
          exnode->cellNeighbors_ = new vector<vector<SimplexId>>();
          getCellNeighbors(exnode);
        }
        return (*(exnode->cellNeighbors_))[localCellId].size();
      }
      
      const vector<vector<SimplexId> > *getCellNeighbors(){
        cellNeighborList_.reserve(cellNumber_);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->cellNeighbors_ == nullptr){
            exnode->cellNeighbors_ = new vector<vector<SimplexId>>();
            getCellNeighbors(exnode);
          }
          cellNeighborList_.insert(cellNeighborList_.end(), exnode->cellNeighbors_->begin(), exnode->cellNeighbors_->end());
        }
        return &cellNeighborList_;
      }
      
      int getCellTriangle(const SimplexId &cellId, 
        const int &localTriangleId, SimplexId &triangleId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
        if((cellId < 0)||(cellId >= cellNumber_))
          return -1;
        if((localTriangleId < 0)||(localTriangleId >= getCellTriangleNumber(cellId)))
          return -2;
        #endif

        SimplexId nid = vertexIndices_[cellArray_[(cellArray_[0]+1)*cellId+1]];
        SimplexId localCellId = cellId-cellIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->cellTriangles_ == nullptr){
          exnode->cellTriangles_ = new vector<vector<SimplexId>>();
          getCellTriangles(exnode);
        }
        triangleId = (*(exnode->cellTriangles_))[localCellId][localTriangleId];
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
        cellTriangleList_.reserve(cellNumber_);
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
          buildInternalEdgeList(nid, &localInternalEdges);
          edgeList_.insert(edgeList_.end(), localInternalEdges.begin(), localInternalEdges.end());
        }
        return &edgeList_;
      }
      
      int getEdgeLink(const SimplexId &edgeId, 
        const int &localLinkId, SimplexId &linkId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return -1;
          if(localLinkId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->edgeLinks_ == nullptr){
          exnode->edgeLinks_ = new vector<vector<SimplexId>>();
          getEdgeLinks(exnode);
        }
        if(localLinkId >= (SimplexId) (*(exnode->edgeLinks_))[localEdgeId].size())
          return -2;
        
        linkId = (*(exnode->edgeLinks_))[localEdgeId][localLinkId];
        return 0;
      }
        
      SimplexId getEdgeLinkNumber(const SimplexId &edgeId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localEdgeId = edgeId-edgeIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->edgeLinks_ == nullptr){
          exnode->edgeLinks_ = new vector<vector<SimplexId>>();
          getEdgeLinks(exnode);
        }
        return (*(exnode->edgeLinks_))[localEdgeId].size();
      }
      
      const vector<vector<SimplexId> > *getEdgeLinks(){
        edgeLinkList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->edgeLinks_ == nullptr){
            exnode->edgeLinks_ = new vector<vector<SimplexId>>();
            getEdgeLinks(exnode);
          }
          edgeLinkList_.insert(edgeLinkList_.end(), exnode->edgeLinks_->begin(), exnode->edgeLinks_->end());
        }
        return &edgeLinkList_;
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
          getEdgeStars(exnode);
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
          getEdgeStars(exnode);
        }
        return (*(exnode->edgeStars_))[localEdgeId].size();
      }
      
      const vector<vector<SimplexId> > *getEdgeStars(){
        edgeStarList_.reserve(edgeIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->edgeStars_ == nullptr){
            exnode->edgeStars_ = new vector<vector<SimplexId>>();
            getEdgeStars(exnode);
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
        return (*(exnode->edgeTriangles_))[edgeId-edgeIntervals_[nid-1]-1].size();
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
          buildInternalEdgeList(nid, exnode->internalEdgeList_);
        }
        if(localVertexId){
          vertexId = (*(exnode->internalEdgeList_))[localEdgeId].second;
        }
        else{
          vertexId = (*(exnode->internalEdgeList_))[localEdgeId].first;
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
        // if it is a triangle mesh
        if(getDimensionality() == 2){
          triangleList_.resize(cellNumber_, vector<SimplexId>(3));
          for(SimplexId cid = 0; cid < cellNumber_; cid++){
            SimplexId cellId = (cellArray_[0]+1)*cid;
            triangleList_[cid][0] = cellArray_[cellId+1];
            triangleList_[cid][1] = cellArray_[cellId+2];
            triangleList_[cid][2] = cellArray_[cellId+3];
          }
        }
        else{
          triangleList_.reserve(triangleIntervals_.back() + 1);
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            ExpandedNode *exnode = searchCache(nid);
            if(exnode->internalTriangleList_ == nullptr){
              exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
              buildInternalTriangleList(nid, exnode->internalTriangleList_);
            }
            triangleList_.insert(triangleList_.end(), exnode->internalTriangleList_->begin(), exnode->internalTriangleList_->end());
          }
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
        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleEdges_ == nullptr){
          exnode->triangleEdges_ = new vector<vector<SimplexId>>();
          getTriangleEdges(exnode);
        }
        edgeId = (*(exnode->triangleEdges_))[localTriangleId][localEdgeId];
        return 0;
      }
      
      SimplexId getTriangleEdgeNumber(const SimplexId &triangleId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        return 3;
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
          
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
          if(localLinkId < 0)
            return -2;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleLinks_ == nullptr){
          exnode->triangleLinks_ = new vector<vector<SimplexId>>();
          getTriangleLinks(exnode);
        }
        if(localLinkId >= (SimplexId) (*(exnode->triangleLinks_))[localTriangleId].size())
          return -2;
        
        linkId = (*(exnode->triangleLinks_))[localTriangleId][localLinkId];
        return 0;
      }
      
      SimplexId getTriangleLinkNumber(const SimplexId &triangleId) const{
                 
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleLinks_ == nullptr){
          exnode->triangleLinks_ = new vector<vector<SimplexId>>();
          getTriangleLinks(exnode);
        }
        return (*(exnode->triangleLinks_))[localTriangleId].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleLinks(){
        triangleLinkList_.reserve(triangleIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->triangleLinks_ == nullptr){
            exnode->triangleLinks_ = new vector<vector<SimplexId>>();
            getTriangleLinks(exnode);
          }
          triangleLinkList_.insert(triangleLinkList_.end(), exnode->triangleLinks_->begin(), exnode->triangleLinks_->end());
        }
        return &triangleLinkList_;
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
          getTriangleStars(exnode);
        }
        if(localStarId >= (SimplexId) (*(exnode->triangleStars_))[localTriangleId].size())
          return -2;
        
        starId = (*(exnode->triangleStars_))[localTriangleId][localStarId];
        return 0;
      }
        
      SimplexId getTriangleStarNumber(const SimplexId &triangleId) const{

         #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->triangleStars_ == nullptr){
          exnode->triangleStars_ = new vector<vector<SimplexId>>();
          getTriangleStars(exnode);
        }
        
        return (*(exnode->triangleStars_))[localTriangleId].size();
      }
      
      const vector<vector<SimplexId> > *getTriangleStars(){
        triangleStarList_.reserve(triangleIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->triangleStars_ == nullptr){
            exnode->triangleStars_ = new vector<vector<SimplexId>>();
            getTriangleStars(exnode);
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
        SimplexId localTriangleId = triangleId-triangleIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->internalTriangleList_ == nullptr){
          exnode->internalTriangleList_ = new vector<vector<SimplexId>>();
          buildInternalTriangleList(nid, exnode->internalTriangleList_);
        }
        vertexId = (*(exnode->internalTriangleList_))[localTriangleId][localVertexId];
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
        
        SimplexId nid = vertexIndices_[vertexId];
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexEdges_ == nullptr){
          exnode->vertexEdges_ = new vector<vector<SimplexId>>();
          getVertexEdges(exnode);
        }
        if(localEdgeId >= (SimplexId) (*exnode->vertexEdges_)[localVertexId].size())
            return -2;
        edgeId = (*(exnode->vertexEdges_))[localVertexId][localEdgeId];
        return 0;
      }
      
      SimplexId getVertexEdgeNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = vertexIndices_[vertexId];
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
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId > vertexIntervals_.back()))
            return -1;
          if(localLinkId < 0)
            return -2;
        #endif

        SimplexId nid = vertexIndices_[vertexId];
        SimplexId localVertexId = vertexId-vertexIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexLinks_ == nullptr){
          exnode->vertexLinks_ = new vector<vector<SimplexId>>();
          getVertexLinks(exnode);
        }
        if(localLinkId >= (SimplexId) (*(exnode->vertexLinks_))[localVertexId].size())
          return -2;
        
        linkId = (*(exnode->vertexLinks_))[localVertexId][localLinkId];
        return 0;
      }
      
      SimplexId getVertexLinkNumber(const SimplexId &vertexId) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId > vertexIntervals_.back()))
            return -1;
        #endif

        SimplexId nid = vertexIndices_[vertexId];
        SimplexId localVertexId = vertexId-vertexIntervals_[nid-1]-1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexLinks_ == nullptr){
          exnode->vertexLinks_ = new vector<vector<SimplexId>>();
          getVertexLinks(exnode);
        }
        return (*(exnode->vertexLinks_))[localVertexId].size();
      }
      
      const vector<vector<SimplexId> > *getVertexLinks(){
        vertexLinkList_.reserve(vertexIntervals_.back()+1);
        for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
          ExpandedNode *exnode = searchCache(nid);
          if(exnode->vertexLinks_ == nullptr){
            exnode->vertexLinks_ = new vector<vector<SimplexId>>();
            getVertexLinks(exnode);
          }
          vertexLinkList_.insert(vertexLinkList_.end(), exnode->vertexLinks_->begin(), exnode->vertexLinks_->end());
        }
        return &vertexLinkList_;
      }
      
      int getVertexNeighbor(const SimplexId &vertexId, 
        const int &localNeighborId, SimplexId &neighborId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
          if(localNeighborId < 0)
            return -2;
        #endif
        
        SimplexId nid = vertexIndices_[vertexId];
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexNeighbors_ == nullptr){
          exnode->vertexNeighbors_ = new vector<vector<SimplexId>>();
          getVertexNeighbors(exnode);
        }
        if(localNeighborId >= (SimplexId) (*(exnode->vertexNeighbors_))[localVertexId].size())
          return -2;  
        neighborId = (*(exnode->vertexNeighbors_))[localVertexId][localNeighborId];
        return 0;
      }
      
      SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = vertexIndices_[vertexId];
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

        SimplexId nid = vertexIndices_[vertexId];
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        if(exnode->vertexStars_ == nullptr){
          exnode->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(exnode);
        }
        if(localStarId >= (SimplexId) (*exnode->vertexStars_)[localVertexId].size())
          return -2;
        starId = (*(exnode->vertexStars_))[localVertexId][localStarId];
        return 0;
      }
      
      SimplexId getVertexStarNumber(const SimplexId &vertexId) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return -1;
        #endif

        SimplexId nid = vertexIndices_[vertexId];
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

        SimplexId nid = vertexIndices_[vertexId];
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

        SimplexId nid = vertexIndices_[vertexId];
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
        #ifndef TTK_ENABLE_KAMIKAZE
          if((edgeId < 0)||(edgeId > edgeIntervals_.back()))
            return false;
        #endif
        SimplexId nid = findNodeIndex(edgeId, EDGE_ID);
        SimplexId localedgeId = edgeId - edgeIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        getBoundaryCells(exnode, 1);
        return (*(exnode->boundaryEdges_))[localedgeId];
      }
        
      bool isEmpty() const{
        return !vertexNumber_;
      }
      
      bool isTriangleOnBoundary(const SimplexId &triangleId) const{
        if(getDimensionality() == 2)
          return false;
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if((triangleId < 0)||(triangleId > triangleIntervals_.back()))
            return false;
        #endif
        SimplexId nid = findNodeIndex(triangleId, TRIANGLE_ID);
        SimplexId localtriangleId = triangleId - triangleIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        getBoundaryCells(exnode);

        return (*(exnode->boundaryTriangles_))[localtriangleId];
      }
      
      bool isVertexOnBoundary(const SimplexId &vertexId) const{
        #ifndef TTK_ENABLE_KAMIKAZE
          if((vertexId < 0)||(vertexId >= vertexNumber_))
            return false;
        #endif
        SimplexId nid = vertexIndices_[vertexId];
        SimplexId localVertexId = vertexId - vertexIntervals_[nid-1] - 1;
        ExpandedNode *exnode = searchCache(nid);
        getBoundaryCells(exnode, 0);
        return (*(exnode->boundaryVertices_))[localVertexId];
      }

      int preprocessBoundaryEdges(){
        if(getDimensionality() == 2 || getDimensionality() == 3){
          preprocessEdges();
          hasPreprocessedBoundaryEdges_ = true;
        }
        else{
          // unsupported dimension
          std::stringstream msg;
          msg << "[ExplicitTopoCluster] Unsupported dimension for boundary "
            << "preprocessing." << std::endl;
          dMsg(std::cerr, msg.str(), infoMsg);
          return -1;
        }
        return 0;
      }
      
      int preprocessBoundaryTriangles(){
        if(getDimensionality() == 2 || getDimensionality() == 3){
          preprocessTriangles();
          hasPreprocessedBoundaryTriangles_ = true;
        }
        else{
          // unsupported dimension
          std::stringstream msg;
          msg << "[ExplicitTopoCluster] Unsupported dimension for boundary "
            << "preprocessing." << std::endl;
          dMsg(std::cerr, msg.str(), infoMsg);
          return -1;
        }
        return 0;
      }
      
      int preprocessBoundaryVertices(){
        preprocessTriangles();
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
          Timer t;
          edgeIntervals_.resize(nodeNumber_+1);
          edgeIntervals_[0] = -1;
          internalEdgeMaps_.resize(nodeNumber_+1);
          vector<SimplexId> edgeCount(nodeNumber_+1);
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            edgeCount[nid] = buildInternalEdgeMap(nid, &internalEdgeMaps_[nid]);
          }

          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            edgeIntervals_[nid] = edgeIntervals_[nid-1]+edgeCount[nid];
          }

          hasPreprocessedEdges_ = true;

          cout << "[ExplicitTopoCluster] Edges processed in " << t.getElapsedTime() << " s.\n";
        }

        return 0;
      }
      
      int preprocessEdgeLinks(){
        if(getDimensionality() == 2 || getDimensionality() == 3){
          preprocessEdges();
          hasPreprocessedEdgeLinks_ = true;
        }
        else{
          // unsupported dimension
          std::stringstream msg;
          msg 
            << "[ExplicitTopoCluster] Unsupported dimension for edge link "
            << "preprocessing." << std::endl;
          dMsg(std::cerr, msg.str(), infoMsg);
          return -1;
        }
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
        if(getDimensionality() == 2){
          hasPreprocessedTriangles_ = true;
          return 0;
        }

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
          Timer t;
          triangleIntervals_.resize(nodeNumber_+1);
          triangleIntervals_[0] = -1;
          internalTriangleMaps_.resize(nodeNumber_+1);
          vector<SimplexId> triangleCount(nodeNumber_+1);
          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            triangleCount[nid] = buildInternalTriangleMap(nid, &internalTriangleMaps_[nid]);
          }

          for(SimplexId nid = 1; nid <= nodeNumber_; nid++){
            triangleIntervals_[nid] = triangleIntervals_[nid-1]+triangleCount[nid];
          }

          hasPreprocessedTriangles_ = true;

          cout << "[ExplicitTopoCluster] Triangles processed in " << t.getElapsedTime() << " s.\n";
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
        if(getDimensionality() == 2){
          preprocessEdges();
          hasPreprocessedVertexLinks_ = true;
        }
        else if(getDimensionality() == 3){
          preprocessTriangles();
          hasPreprocessedVertexLinks_ = true;
        }
        else{
          // unsupported dimension
          std::stringstream msg;
          msg << "[ExplicitTopoCluster] Unsupported dimension for vertex"
            << " link preprocessing." << std::endl;
          dMsg(std::cerr, msg.str(), infoMsg);
          return -1;
        }
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
      void initCache(const size_t size=10){
        cacheSize_ = size;
        cache_.clear();
        cacheMap_.clear();
      }

    protected:

      int clear();

      /**
       * Find the corresponding node index given the id.
       */ 
      SimplexId findNodeIndex(SimplexId id, int idType) const{
        const vector<SimplexId> *intervals = nullptr;
        // determine which vector to search
        if(idType == EDGE_ID){
          intervals = &edgeIntervals_;
        }else if(idType == TRIANGLE_ID){
          intervals = &triangleIntervals_;
        }else{
          return -1;
        }

        vector<SimplexId>::const_iterator low = lower_bound(intervals->begin(), intervals->end(), id);
        return (low-intervals->begin());
      }

      /**
       * Search the node in the cache.
       */
      ExpandedNode* searchCache(const SimplexId &nodeId) const{
        // cannot find the expanded node in the cache
        if(cacheMap_.find(nodeId) == cacheMap_.end()){
          // missCount_++;
          if(cache_.size() >= cacheSize_){
            cacheMap_.erase(cache_.back()->nid);
            delete cache_.back();
            cache_.pop_back();
          }
          cache_.push_front(new ExpandedNode(nodeId));
          cacheMap_[nodeId] = cache_.begin();
          return cache_.front();
        }
        return (*cacheMap_[nodeId]);
      }

      /** 
       * Build the internal edge list in the node.
       */
      int buildInternalEdgeList(SimplexId nodeId,
        vector<pair_int> * const internalEdgeList) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        internalEdgeList->clear();
        internalEdgeList->resize(internalEdgeMaps_[nodeId].size());
        for(auto iter = internalEdgeMaps_[nodeId].begin(); iter != internalEdgeMaps_[nodeId].end(); iter++){
          (*internalEdgeList)[iter->second - 1] = iter->first;
        }

        return 0;
      }

      /** 
       * Build the internal edge map in the node.
       */ 
      int buildInternalEdgeMap(SimplexId nodeId, 
        boost::unordered_map<pair_int, SimplexId> * const internalEdgeMap) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId edgeCount = 0, verticesPerCell = cellArray_[0];
        internalEdgeMap->clear();

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
              if(internalEdgeMap->find(edgeIds) == internalEdgeMap->end()){
                edgeCount++;
                (*internalEdgeMap)[edgeIds] = edgeCount;
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
                if(internalEdgeMap->find(edgeIds) == internalEdgeMap->end()){
                  edgeCount++;
                  (*internalEdgeMap)[edgeIds] = edgeCount;
                }
              }
            }
          }
        }
        return edgeCount;
      }

      /** 
       * Build the external edge list in the node.
       */ 
      int buildExternalEdgeMap(SimplexId nodeId, 
        boost::unordered_map<pair_int, SimplexId> * const externalEdgeMap) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
          if(!externalEdgeMap)
            return -2;
        #endif

        SimplexId verticesPerCell = cellArray_[0];
        
        // loop through the external cell list
        for(size_t i = 0; i < externalCells_[nodeId].size(); i++){
          pair_int edgeIds;
          SimplexId cellId = (verticesPerCell + 1) * externalCells_[nodeId][i];

          // loop through each edge of the cell
          for(SimplexId j = 0; j < verticesPerCell-1; j++){
            for(SimplexId k = j+1; k < verticesPerCell; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // check if the edge is an external edge
              if(edgeIds.first <= vertexIntervals_[nodeId-1] && edgeIds.second > vertexIntervals_[nodeId-1] && edgeIds.second <= vertexIntervals_[nodeId]){
                SimplexId nodeNum = vertexIndices_[edgeIds.first];
                (*externalEdgeMap)[edgeIds] = internalEdgeMaps_[nodeNum].at(edgeIds) + edgeIntervals_[nodeNum-1];
              }
            }
          }
        }

        return 0;
      }

      /**
       * Build the internal triangle list in the node.
       */
      int buildInternalTriangleList(SimplexId nodeId, 
        vector<vector<SimplexId>> * const internalTriangleList) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        internalTriangleList->clear();
        internalTriangleList->resize(internalTriangleMaps_[nodeId].size());
        for(auto iter = internalTriangleMaps_[nodeId].begin(); iter != internalTriangleMaps_[nodeId].end(); iter++){
          (*internalTriangleList)[iter->second - 1] = iter->first;
        }

        return 0;
      }

      /**
       * Build the internal triangle map in the node.
       */
      int buildInternalTriangleMap(SimplexId nodeId, 
        boost::unordered_map<vector<SimplexId>, SimplexId> * const internalTriangleMap) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId triangleCount = 0, verticesPerCell = cellArray_[0];
        internalTriangleMap->clear();

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
                
                if(internalTriangleMap->find(triangleIds) == internalTriangleMap->end()){
                  triangleCount++;
                  (*internalTriangleMap)[triangleIds] = triangleCount;
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

                  if(internalTriangleMap->find(triangleIds) == internalTriangleMap->end()){
                    triangleCount++;
                    (*internalTriangleMap)[triangleIds] = triangleCount;
                  }
                }
              }
            }
          }
        }

        return triangleCount;
      }

      /**
       * Build the external triangle list in the node.
       */
      int buildExternalTriangleMap(SimplexId nodeId, 
        boost::unordered_map<vector<SimplexId>, SimplexId> * const externalTriangleMap) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodeId <= 0 || nodeId > nodeNumber_)
            return -1;
        #endif

        SimplexId verticesPerCell = cellArray_[0];

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodeId]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            if(triangleIds[0] <= vertexIntervals_[nodeId-1]){
              for(SimplexId k = j+1; k < verticesPerCell-1; k++){
                for(SimplexId l = k+1; l < verticesPerCell; l++){
                  triangleIds[1] = cellArray_[cellId + k + 1];
                  triangleIds[2] = cellArray_[cellId + l + 1];

                  if(triangleIds[1] > vertexIntervals_[nodeId-1] && triangleIds[1] <= vertexIntervals_[nodeId]){
                    SimplexId nodeNum = vertexIndices_[triangleIds[0]];
                    (*externalTriangleMap)[triangleIds] = internalTriangleMaps_[nodeNum].at(triangleIds) + triangleIntervals_[nodeNum-1];
                  }
                  else if(triangleIds[2] > vertexIntervals_[nodeId-1] && triangleIds[2] <= vertexIntervals_[nodeId]){
                    SimplexId nodeNum = vertexIndices_[triangleIds[0]];
                    (*externalTriangleMap)[triangleIds] = internalTriangleMaps_[nodeNum].at(triangleIds) + triangleIntervals_[nodeNum-1];
                  }
                }
              }
            }
          }
        }

        return 0;
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

        for(SimplexId i = cellIntervals_[nodePtr->nid-1]+1; i <= cellIntervals_[nodePtr->nid]; i++){
          SimplexId cellId = (cellArray_[0]+1)*i;
          int cnt = 0;
          // get the internal edge id from the map
          for(SimplexId k = 1; k < cellArray_[0]; k++){
            pair_int edgePair(cellArray_[cellId+1], cellArray_[cellId+k+1]);
            (*(nodePtr->cellEdges_))[i-cellIntervals_[nodePtr->nid-1]-1][cnt++] = internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1];
          }
          for(SimplexId j = 1; j < cellArray_[0]-1; j++){
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              pair_int edgePair(cellArray_[cellId+j+1], cellArray_[cellId+k+1]);
              if(edgePair.first <= vertexIntervals_[nodePtr->nid]){
                (*(nodePtr->cellEdges_))[i-cellIntervals_[nodePtr->nid-1]-1][cnt++] = internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1];
              }
              else{
                SimplexId nodeNum = vertexIndices_[edgePair.first];
                (*(nodePtr->cellEdges_))[i-cellIntervals_[nodePtr->nid-1]-1][cnt++] = internalEdgeMaps_[nodeNum].at(edgePair) + edgeIntervals_[nodeNum-1];
              }
            }
          }
        }

        return 0;
      }

      /**
       * Get the cell triangles for all cells in a given node. 
       */
      int getCellNeighbors(ExpandedNode * const nodePtr) const{
        
        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->cellNeighbors_->clear();
        nodePtr->cellNeighbors_->resize(cellIntervals_[nodePtr->nid]-cellIntervals_[nodePtr->nid-1]);

        if(nodePtr->vertexStars_ == nullptr){
          nodePtr->vertexStars_ = new vector<vector<SimplexId>>();
          getVertexStars(nodePtr);
          for(int i = 0; i < (int) nodePtr->vertexStars_->size(); i++){
            sort((*(nodePtr->vertexStars_))[i].begin(), (*(nodePtr->vertexStars_))[i].end());
          }
        }

        boost::unordered_map<SimplexId, ExpandedNode*> nodeMaps;
        for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
          SimplexId cellId = (cellArray_[0] + 1) * cid;
          for(SimplexId j = 1; j < cellArray_[0]; j++){
            if(cellArray_[cellId+j+1] > vertexIntervals_[nodePtr->nid]){
              SimplexId nodeId = vertexIndices_[cellArray_[cellId+j+1]];
              if(nodeMaps.find(nodeId) == nodeMaps.end()){
                ExpandedNode *newNode = new ExpandedNode(nodeId);
                newNode->vertexStars_ = new vector<vector<SimplexId>>();
                getVertexStars(newNode);
                for(int i = 0; i < (int) newNode->vertexStars_->size(); i++){
                  sort((*(newNode->vertexStars_))[i].begin(), (*(newNode->vertexStars_))[i].end());
                }
                nodeMaps[nodeId] = newNode;
              }
            }
          }
        }

        if(getDimensionality() == 2){
          for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
            for(SimplexId j = 0; j < 3; j++){
              
              SimplexId v0 = cellArray_[4*cid + 1 + j];
              SimplexId v1 = cellArray_[4*cid + 1 + (j+1)%3];

              vector<SimplexId> star0, star1;
              if(v0 <= vertexIntervals_[nodePtr->nid]){
                star0 = (*(nodePtr->vertexStars_))[v0-vertexIntervals_[nodePtr->nid-1]-1];
              }
              else{
                SimplexId nid = vertexIndices_[v0];
                star0 = (*(nodeMaps[nid]->vertexStars_))[v0-vertexIntervals_[nid-1]-1];
              }
              if(v1 <= vertexIntervals_[nodePtr->nid]){
                star1 = (*(nodePtr->vertexStars_))[v1-vertexIntervals_[nodePtr->nid-1]-1];
              }
              else{
                SimplexId nid = vertexIndices_[v1];
                star1 = (*(nodeMaps[nid]->vertexStars_))[v1-vertexIntervals_[nid-1]-1];
              }
              
              // perform an intersection of the 2 sorted star lists
              SimplexId pos0 = 0, pos1 = 0;
              SimplexId intersection = -1;

              while((pos0 < (SimplexId) star0.size()) && (pos1 < (SimplexId) star1.size())){
                
                SimplexId biggest = star0[pos0];
                if(star1[pos1] > biggest){
                  biggest = star1[pos1];
                }
                
                for(SimplexId l = pos0; l < (SimplexId) star0.size(); l++){
                  if(star0[l] < biggest){
                    pos0++;
                  }
                  else{
                    break;
                  }
                }
                for(SimplexId l = pos1; l < (SimplexId) star1.size(); l++){
                  if(star1[l] < biggest){
                    pos1++;
                  }
                  else{
                    break;
                  }
                }
                
                if(pos0 >= (SimplexId) star0.size() || pos1 >= (SimplexId) star1.size())
                  break;

                if(star0[pos0] == star1[pos1]){
                  if(star0[pos0] != cid){
                    intersection = star0[pos0];
                    break;
                  }
                  pos0++;
                  pos1++;
                }
              }
              
              if(intersection != -1){
                (*(nodePtr->cellNeighbors_))[cid-cellIntervals_[nodePtr->nid-1]-1].push_back(intersection);
              }
            }
          }
        }

        else if(getDimensionality() == 3){
          for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
            // go triangle by triangle
            for(SimplexId j = 0; j < 4; j++){
            
              SimplexId v0 = cellArray_[5*cid + 1 + j%4];
              SimplexId v1 = cellArray_[5*cid + 1 + (j+1)%4];
              SimplexId v2 = cellArray_[5*cid + 1 + (j+2)%4];

              vector<SimplexId> star0, star1, star2;
              if(v0 <= vertexIntervals_[nodePtr->nid]){
                star0 = (*(nodePtr->vertexStars_))[v0-vertexIntervals_[nodePtr->nid-1]-1];
              }
              else{
                SimplexId nid = vertexIndices_[v0];
                star0 = (*(nodeMaps[nid]->vertexStars_))[v0-vertexIntervals_[nid-1]-1];
              }
              if(v1 <= vertexIntervals_[nodePtr->nid]){
                star1 = (*(nodePtr->vertexStars_))[v1-vertexIntervals_[nodePtr->nid-1]-1];
              }
              else{
                SimplexId nid = vertexIndices_[v1];
                star1 = (*(nodeMaps[nid]->vertexStars_))[v1-vertexIntervals_[nid-1]-1];
              }
              if(v2 <= vertexIntervals_[nodePtr->nid]){
                star2 = (*(nodePtr->vertexStars_))[v2-vertexIntervals_[nodePtr->nid-1]-1];
              }
              else{
                SimplexId nid = vertexIndices_[v2];
                star2 = (*(nodeMaps[nid]->vertexStars_))[v2-vertexIntervals_[nid-1]-1];
              }
              
              // perform an intersection of the 3 (sorted) star lists
              SimplexId pos0 = 0, pos1 = 0, pos2 = 0;
              SimplexId intersection = -1;
              
              while((pos0 < (SimplexId) star0.size())
                &&(pos1 < (SimplexId) star1.size())
                &&(pos2 < (SimplexId) star2.size())){
                
                SimplexId biggest = star0[pos0];
                if(star1[pos1] > biggest){
                  biggest = star1[pos1];
                }
                if(star2[pos2] > biggest){
                  biggest = star2[pos2];
                }
                
                for(SimplexId l = pos0; l < (SimplexId) star0.size(); l++){
                  if(star0[l] < biggest){
                    pos0++;
                  }
                  else{
                    break;
                  }
                }
                for(SimplexId l = pos1; l < (SimplexId) star1.size(); l++){
                  if(star1[l] < biggest){
                    pos1++;
                  }
                  else{
                    break;
                  }
                }
                for(SimplexId l = pos2; l < (SimplexId) star2.size(); l++){
                  if(star2[l] < biggest){
                    pos2++;
                  }
                  else{
                    break;
                  }
                }
                
                if(pos0 >= (SimplexId) star0.size() 
                  || pos1 >= (SimplexId) star1.size() 
                  || pos2 >= (SimplexId) star2.size())
                  break;  
                
                if((star0[pos0] == star1[pos1]) && (star0[pos0] == star2[pos2])){
                  if(star0[pos0] != cid){
                    intersection = star0[pos0];
                    break;
                  }
                  pos0++; pos1++; pos2++;
                }
              }

              if(intersection != -1){
                (*(nodePtr->cellNeighbors_))[cid-cellIntervals_[nodePtr->nid-1]-1].push_back(intersection);
              }
            }
          }
        }

        // release the memory
        for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++){
          delete iter->second;
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
        nodePtr->cellTriangles_->clear();
        nodePtr->cellTriangles_->resize(cellIntervals_[nodePtr->nid]-cellIntervals_[nodePtr->nid-1], vector<SimplexId>(trianglesPerCell));

        for(SimplexId i = cellIntervals_[nodePtr->nid-1]+1; i <= cellIntervals_[nodePtr->nid]; i++){
          SimplexId cellId = (cellArray_[0]+1)*i;
          vector<SimplexId> triangleVec(3);
          // get the internal triangle from the map
          triangleVec[0] = cellArray_[cellId+1]; 
          for(SimplexId k = 1; k < cellArray_[0]-1; k++){
            triangleVec[1] = cellArray_[cellId+1+k];
            for(SimplexId l = k+1; l < cellArray_[0]; l++){
              triangleVec[2] = cellArray_[cellId+1+l];
              (*(nodePtr->cellTriangles_))[i-cellIntervals_[nodePtr->nid-1]-1][k+l-3] = internalTriangleMaps_[nodePtr->nid].at(triangleVec) + triangleIntervals_[nodePtr->nid-1];
            }
          }
          // group the external triangles by node id
          triangleVec[0] = cellArray_[cellId+2];
          triangleVec[1] = cellArray_[cellId+3];
          triangleVec[2] = cellArray_[cellId+4];
          if(triangleVec[0] <= vertexIntervals_[nodePtr->nid]){
            (*(nodePtr->cellTriangles_))[i-cellIntervals_[nodePtr->nid-1]-1].back() = internalTriangleMaps_[nodePtr->nid].at(triangleVec) + triangleIntervals_[nodePtr->nid-1];
          }
          else{
            SimplexId nodeNum = vertexIndices_[triangleVec[0]];
            (*(nodePtr->cellTriangles_))[i-cellIntervals_[nodePtr->nid-1]-1].back() = internalTriangleMaps_[nodeNum].at(triangleVec) + triangleIntervals_[nodeNum-1];
          }
        }

        return 0;
      }

      /**
       * Get the triangle links for all the triangles in a given node.
       */
      int getEdgeLinks(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->edgeLinks_->clear();
        nodePtr->edgeLinks_->resize(edgeIntervals_[nodePtr->nid]-edgeIntervals_[nodePtr->nid-1]);
        if(getDimensionality() == 2){
          if(nodePtr->edgeStars_ == nullptr){
            nodePtr->edgeStars_ = new vector<vector<SimplexId>>();
            getEdgeStars(nodePtr);
          }
          boost::unordered_map<pair_int, SimplexId>::const_iterator iter;
          for(iter = internalEdgeMaps_[nodePtr->nid].begin(); iter != internalEdgeMaps_[nodePtr->nid].end(); iter++){
            for(SimplexId j = 0; j < (SimplexId) (*(nodePtr->edgeStars_))[iter->second-1].size(); j++){
              SimplexId vertexId = -1;
              for(int k = 0; k < 3; k++){
                if((cellArray_[(*(nodePtr->edgeStars_))[iter->second-1][j] + 1 + k] != iter->first.first) &&
                  (cellArray_[(*(nodePtr->edgeStars_))[iter->second-1][j] + 1 + k] != iter->first.second)){
                  vertexId = cellArray_[(*(nodePtr->edgeStars_))[iter->second-1][j] + 1 + k];
                  break;
                }
              }
              if(vertexId != -1){
                (*(nodePtr->edgeLinks_))[iter->second-1].push_back(vertexId);
              }
            }
          }
        }
        else if(getDimensionality() == 3){
          if(nodePtr->cellEdges_ == nullptr){
            nodePtr->cellEdges_ = new vector<vector<SimplexId>>();
            getCellEdges(nodePtr);
          }

          for(SimplexId cid = 0; cid < cellIntervals_[nodePtr->nid]-cellIntervals_[nodePtr->nid-1]; cid++){
            SimplexId cellId = (cid+cellIntervals_[nodePtr->nid-1]+1) * 5;
            pair_int edgePair; edgePair.first = cellArray_[cellId+1];
            for(SimplexId j = 1; j < 4; j++){
              edgePair.second = cellArray_[cellId+j+1];
              (*(nodePtr->edgeLinks_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1].push_back((*(nodePtr->cellEdges_))[cid][6-j]);
            }
            if(cellArray_[cellId+2] <= vertexIntervals_[nodePtr->nid]){
              edgePair.first = cellArray_[cellId+2];
              for(int j = 2; j < 4; j++){
                edgePair.second = cellArray_[cellId+j+1];
                (*(nodePtr->edgeLinks_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1].push_back((*(nodePtr->cellEdges_))[cid][4-j]);
              }
              if(cellArray_[cellId+3] <= vertexIntervals_[nodePtr->nid]){
                edgePair = pair_int(cellArray_[cellId+3], cellArray_[cellId+4]);
                (*(nodePtr->edgeLinks_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1].push_back((*(nodePtr->cellEdges_))[cid][0]);
              }
            }
          }

          // loop through the external cell list
          for(SimplexId cid : externalCells_[nodePtr->nid]){
            pair_int edgeIds;
            SimplexId cellId = 5 * cid;

            // loop through each edge of the cell
            for(SimplexId j = 1; j < 4; j++){
              for(SimplexId k = j+1; k < 5; k++){
                edgeIds.first = cellArray_[cellId+j];
                edgeIds.second = cellArray_[cellId+k];
                
                // the edge is in the current node
                if(edgeIds.first > vertexIntervals_[nodePtr->nid-1] && edgeIds.first <= vertexIntervals_[nodePtr->nid]){
                  pair_int otherEdge(-1, -1);
                  for(int i = 1; i < 5; i++){
                    if(cellArray_[cellId+i] != edgeIds.first && cellArray_[cellId+i] != edgeIds.second){
                      if(otherEdge.first == -1){
                        otherEdge.first = cellArray_[cellId+i];
                      }
                      else if(otherEdge.second == -1){
                        otherEdge.second = cellArray_[cellId+i];
                      }
                      else{
                        cerr << "[ExplicitTopoCluster] More than two other vertices are found in the edge!\n";
                      }
                    }
                  }
                  SimplexId nodeId = vertexIndices_[otherEdge.first];
                  (*(nodePtr->edgeLinks_))[internalEdgeMaps_[nodePtr->nid].at(edgeIds)-1].push_back(internalEdgeMaps_[nodeId].at(otherEdge));
                }
              }
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

        if(nodePtr->externalTriangleMap_ == nullptr){
          nodePtr->externalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildExternalTriangleMap(nodePtr->nid, nodePtr->externalTriangleMap_);
        }

        // for internal triangles'
        boost::unordered_map<vector<SimplexId>, SimplexId>::const_iterator iter;
        for(iter = internalTriangleMaps_[nodePtr->nid].begin(); iter != internalTriangleMaps_[nodePtr->nid].end(); iter++){
          pair_int edge1(iter->first[0], iter->first[1]);
          pair_int edge2(iter->first[0], iter->first[2]);

          (*(nodePtr->edgeTriangles_))[internalEdgeMaps_[nodePtr->nid].at(edge1)-1].push_back(iter->second + triangleIntervals_[nodePtr->nid-1]);
          (*(nodePtr->edgeTriangles_))[internalEdgeMaps_[nodePtr->nid].at(edge2)-1].push_back(iter->second + triangleIntervals_[nodePtr->nid-1]);

          if(iter->first[1] <= vertexIntervals_[nodePtr->nid]){
            edge2.first = iter->first[1];
            (*(nodePtr->edgeTriangles_))[internalEdgeMaps_[nodePtr->nid].at(edge2)-1].push_back(iter->second + triangleIntervals_[nodePtr->nid-1]);
          }
        }

        // for external triangles
        for(iter = nodePtr->externalTriangleMap_->begin(); iter != nodePtr->externalTriangleMap_->end(); iter++){
          if(iter->first[1] > vertexIntervals_[nodePtr->nid-1] && iter->first[1] <= vertexIntervals_[nodePtr->nid]){
            pair_int edge(iter->first[1], iter->first[2]);
            (*(nodePtr->edgeTriangles_))[internalEdgeMaps_[nodePtr->nid].at(edge)-1].push_back(iter->second);
          }
        }

        return 0;
      }

      /**
       * Get the edge stars for all the edges in a given node.
       */
      int getEdgeStars(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
              return -1;
        #endif

        nodePtr->edgeStars_->clear();
        nodePtr->edgeStars_->resize(edgeIntervals_[nodePtr->nid]-edgeIntervals_[nodePtr->nid-1]);

        // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
          pair_int edgeIds;
          SimplexId cellId = (cellArray_[0] + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < cellArray_[0]-1; j++){
            edgeIds.first = cellArray_[cellId + j + 1];
            // the edge does not belong to the current node
            if(edgeIds.first > vertexIntervals_[nodePtr->nid]){
              break;
            }
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              edgeIds.second = cellArray_[cellId + k + 1];
              (*(nodePtr->edgeStars_))[internalEdgeMaps_[nodePtr->nid].at(edgeIds)-1].push_back(cid);
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]){
          pair_int edgeIds;
          SimplexId cellId = (cellArray_[0] + 1) * cid;

          // loop through each edge of the cell
          for(SimplexId j = 0; j < cellArray_[0]-1; j++){
            for(SimplexId k = j+1; k < cellArray_[0]; k++){
              edgeIds.first = cellArray_[cellId + j + 1];
              edgeIds.second = cellArray_[cellId + k + 1];
              
              // the edge is in the current node
              if(edgeIds.first > vertexIntervals_[nodePtr->nid-1] && edgeIds.first <= vertexIntervals_[nodePtr->nid]){
                (*(nodePtr->edgeStars_))[internalEdgeMaps_[nodePtr->nid].at(edgeIds)-1].push_back(cid);
              }
            }
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

        boost::unordered_map<vector<SimplexId>, SimplexId>::const_iterator iter;
        for(iter = internalTriangleMaps_[nodePtr->nid].begin(); iter != internalTriangleMaps_[nodePtr->nid].end(); iter++){
          pair_int edgePair(iter->first[0], iter->first[1]);
          (*(nodePtr->triangleEdges_))[iter->second-1][0] = internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1];
          edgePair.second = iter->first[2];
          (*(nodePtr->triangleEdges_))[iter->second-1][1] = internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1];
          edgePair.first = iter->first[1];
          if(edgePair.first > vertexIntervals_[nodePtr->nid-1] && edgePair.first <= vertexIntervals_[nodePtr->nid]){
            (*(nodePtr->triangleEdges_))[iter->second-1][2] = internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1];
          }else{
            SimplexId nodeNum = vertexIndices_[edgePair.first];
            (*(nodePtr->triangleEdges_))[iter->second-1][2] = internalEdgeMaps_[nodeNum].at(edgePair) + edgeIntervals_[nodeNum-1];
          }
        }

        return 0;
      }
      
      /**
       * Get the triangle links for all the triangles in a given node.
       */
      int getTriangleLinks(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->triangleLinks_->clear();
        nodePtr->triangleLinks_->resize(triangleIntervals_[nodePtr->nid]-triangleIntervals_[nodePtr->nid-1]);

        if(nodePtr->triangleStars_ == nullptr){
          nodePtr->triangleStars_ = new vector<vector<SimplexId>>();
          getTriangleStars(nodePtr);
        }

        boost::unordered_map<vector<SimplexId>, SimplexId>::const_iterator iter;
        for(iter = internalTriangleMaps_[nodePtr->nid].begin(); iter != internalTriangleMaps_[nodePtr->nid].end(); iter++){
          for(SimplexId i = 0; i < (SimplexId) (*(nodePtr->triangleStars_))[iter->second-1].size(); i++){
            for(int j = 0; j < 4; j++){
              SimplexId vertexId = cellArray_[(*(nodePtr->triangleStars_))[iter->second-1][i]*5+j+1];
              if((vertexId != iter->first[0]) && (vertexId != iter->first[1]) && (vertexId != iter->first[2])){
                (*(nodePtr->triangleLinks_))[iter->second-1].push_back(vertexId);
                break;
              }
            }
          }
        }

        return 0;
      }

      /**
       * Get the triangle edges for all the triangles in a given node.
       */
      int getTriangleStars(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->triangleStars_->clear();
        nodePtr->triangleStars_->resize(triangleIntervals_[nodePtr->nid]-triangleIntervals_[nodePtr->nid-1]);

        SimplexId verticesPerCell = cellArray_[0];

         // loop through the internal cell list
        for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            // the triangle does not belong to the current node
            if(triangleIds[0] > vertexIntervals_[nodePtr->nid]){
              break;
            }
            for(SimplexId k = j+1; k < verticesPerCell-1; k++){
              for(SimplexId l = k+1; l < verticesPerCell; l++){
                triangleIds[1] = cellArray_[cellId + k + 1];
                triangleIds[2] = cellArray_[cellId + l + 1];
                (*(nodePtr->triangleStars_))[internalTriangleMaps_[nodePtr->nid].at(triangleIds)-1].push_back(cid);
              }
            }
          }
        }

        // loop through the external cell list
        for(SimplexId cid : externalCells_[nodePtr->nid]){
          vector<SimplexId> triangleIds(3);
          SimplexId cellId = (verticesPerCell + 1) * cid;

          // loop through each triangle of the cell
          for(SimplexId j = 0; j < verticesPerCell-2; j++){
            triangleIds[0] = cellArray_[cellId + j + 1];
            if(triangleIds[0] > vertexIntervals_[nodePtr->nid-1] && triangleIds[0] <= vertexIntervals_[nodePtr->nid]){
              for(SimplexId k = j+1; k < verticesPerCell-1; k++){
                for(SimplexId l = k+1; l < verticesPerCell; l++){
                  triangleIds[1] = cellArray_[cellId + k + 1];
                  triangleIds[2] = cellArray_[cellId + l + 1];
                  (*(nodePtr->triangleStars_))[internalTriangleMaps_[nodePtr->nid].at(triangleIds)-1].push_back(cid);
                }
              }
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

        (nodePtr->vertexEdges_)->clear();
        (nodePtr->vertexEdges_)->resize(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1]);

        if(nodePtr->externalEdgeMap_ == nullptr){
          nodePtr->externalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildExternalEdgeMap(nodePtr->nid, nodePtr->externalEdgeMap_);
        }

        boost::unordered_map<pair_int, SimplexId>::const_iterator iter;
        for(iter = internalEdgeMaps_[nodePtr->nid].begin(); iter != internalEdgeMaps_[nodePtr->nid].end(); iter++){
          (*(nodePtr->vertexEdges_))[iter->first.first-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second + edgeIntervals_[nodePtr->nid-1]);
          // the second vertex id of the edge must be greater than the first one
          if(iter->first.second <= vertexIntervals_[nodePtr->nid]){
            (*(nodePtr->vertexEdges_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second + edgeIntervals_[nodePtr->nid-1]);
          }
        }

        for(iter = nodePtr->externalEdgeMap_->begin(); iter != nodePtr->externalEdgeMap_->end(); iter++){
          (*(nodePtr->vertexEdges_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second);
        }

        return 0;
      }

      /**
       * Get the vertex links for all the vertices in a given node.
       */
      int getVertexLinks(ExpandedNode * const nodePtr) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif

        nodePtr->vertexLinks_->clear();
        nodePtr->vertexLinks_->resize(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1]);

        if(getDimensionality() == 2){
          for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
            SimplexId cellId = 4 * cid;

            // the first vertex of the cell must be in the cluster
            pair_int edgePair(cellArray_[cellId+2], cellArray_[cellId+3]);
            SimplexId nodeId = vertexIndices_[cellArray_[cellId+2]];
            (*(nodePtr->vertexLinks_))[cellArray_[cellId+1]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
              internalEdgeMaps_[nodeId].at(edgePair) + edgeIntervals_[nodeId-1]);
            if(cellArray_[cellId+2] <= vertexIntervals_[nodePtr->nid]){
              edgePair.first = cellArray_[cellId+1];
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+2]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
              internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1]);
              if(cellArray_[cellId+3] <= vertexIntervals_[nodePtr->nid]){
                edgePair.second = cellArray_[cellId+2];
                (*(nodePtr->vertexLinks_))[cellArray_[cellId+3]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalEdgeMaps_[nodePtr->nid].at(edgePair) + edgeIntervals_[nodePtr->nid-1]);
              }
            }
          }

          // loop through the external cell list
          for(SimplexId cid : externalCells_[nodePtr->nid]){
            SimplexId cellId = 4 * cid;
            pair_int edgePair(cellArray_[cellId+1], cellArray_[cellId+3]);
            SimplexId nodeId = vertexIndices_[edgePair.first];
            if(cellArray_[cellId+2] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+2] <= vertexIntervals_[nodePtr->nid]){
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+2]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalEdgeMaps_[nodeId].at(edgePair)+edgeIntervals_[nodeId-1]);
            }
            if(cellArray_[cellId+3] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+3] <= vertexIntervals_[nodePtr->nid]){
              edgePair.second = cellArray_[cellId+2];
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+3]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalEdgeMaps_[nodeId].at(edgePair)+edgeIntervals_[nodeId-1]);
            }
          }
        }
        else if(getDimensionality() == 3){
          for(SimplexId cid = cellIntervals_[nodePtr->nid-1]+1; cid <= cellIntervals_[nodePtr->nid]; cid++){
            SimplexId cellId = 5 * cid;

            // v1: (v2, v3, v4)
            vector<SimplexId> triangleVec(3, -1);
            triangleVec[0] = cellArray_[cellId+2];
            triangleVec[1] = cellArray_[cellId+3];
            triangleVec[2] = cellArray_[cellId+4];
            SimplexId nodeId = vertexIndices_[cellArray_[cellId+2]];
            (*(nodePtr->vertexLinks_))[cellArray_[cellId+1]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
              internalTriangleMaps_[nodeId].at(triangleVec) + triangleIntervals_[nodeId-1]);
            // v2: (v1, v3, v4)
            if(cellArray_[cellId+2] <= vertexIntervals_[nodePtr->nid]){
              triangleVec[0] = cellArray_[cellId+1];
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+2]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalTriangleMaps_[nodePtr->nid].at(triangleVec) + triangleIntervals_[nodePtr->nid-1]);
              // v3: (v1, v2, v4)
              if(cellArray_[cellId+3] <= vertexIntervals_[nodePtr->nid]){
                triangleVec[1] = cellArray_[cellId+2];
                (*(nodePtr->vertexLinks_))[cellArray_[cellId+3]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                  internalTriangleMaps_[nodePtr->nid].at(triangleVec) + triangleIntervals_[nodePtr->nid-1]);
              }
              // v4: (v1, v2, v3)
              if(cellArray_[cellId+4] <= vertexIntervals_[nodePtr->nid]){
                triangleVec[2] = cellArray_[cellId+3];
                (*(nodePtr->vertexLinks_))[cellArray_[cellId+4]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                  internalTriangleMaps_[nodePtr->nid].at(triangleVec) + triangleIntervals_[nodePtr->nid-1]);
              }
            }
          }

          // loop through the external cell list
          for(SimplexId cid : externalCells_[nodePtr->nid]){
            SimplexId cellId = 5 * cid;
            // start from v2
            vector<SimplexId> triangleVec(3, -1);
            triangleVec[0] = cellArray_[cellId+1];
            triangleVec[1] = cellArray_[cellId+3];
            triangleVec[2] = cellArray_[cellId+4];
            SimplexId nodeId = vertexIndices_[triangleVec[0]];
            if(cellArray_[cellId+2] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+2] <= vertexIntervals_[nodePtr->nid]){
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+2]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalTriangleMaps_[nodeId].at(triangleVec) + triangleIntervals_[nodeId-1]);
            }
            if(cellArray_[cellId+3] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+3] <= vertexIntervals_[nodePtr->nid]){
              triangleVec[1] = cellArray_[cellId+2];
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+3]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalTriangleMaps_[nodeId].at(triangleVec) + triangleIntervals_[nodeId-1]);
            }
            if(cellArray_[cellId+4] > vertexIntervals_[nodePtr->nid-1] && cellArray_[cellId+4] <= vertexIntervals_[nodePtr->nid]){
              triangleVec[1] = cellArray_[cellId+2];
              triangleVec[2] = cellArray_[cellId+3];
              (*(nodePtr->vertexLinks_))[cellArray_[cellId+4]-vertexIntervals_[nodePtr->nid-1]-1].push_back(
                internalTriangleMaps_[nodeId].at(triangleVec) + triangleIntervals_[nodeId-1]);
            }
          }
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

        if(nodePtr->externalEdgeMap_ == nullptr){
          nodePtr->externalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
          buildExternalEdgeMap(nodePtr->nid, nodePtr->externalEdgeMap_);
        }

        boost::unordered_map<pair_int, SimplexId>::const_iterator iter;
        for(iter = internalEdgeMaps_[nodePtr->nid].begin(); iter != internalEdgeMaps_[nodePtr->nid].end(); iter++){
          (*(nodePtr->vertexNeighbors_))[iter->first.first-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->first.second);
          if(iter->first.second <= vertexIntervals_[nodePtr->nid])
            (*(nodePtr->vertexNeighbors_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->first.first);
        }

        for(iter = nodePtr->externalEdgeMap_->begin(); iter != nodePtr->externalEdgeMap_->end(); iter++){
          (*(nodePtr->vertexNeighbors_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->first.first);
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

        if(nodePtr->externalTriangleMap_ == nullptr){
          nodePtr->externalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
          buildExternalTriangleMap(nodePtr->nid, nodePtr->externalTriangleMap_);
        }


        boost::unordered_map<vector<SimplexId>, SimplexId>::const_iterator iter;
        for(iter = internalTriangleMaps_[nodePtr->nid].begin(); iter != internalTriangleMaps_[nodePtr->nid].end(); iter++){
          for(SimplexId j = 0; j < 3; j++){
            if(iter->first[j] > vertexIntervals_[nodePtr->nid-1] && iter->first[j] <= vertexIntervals_[nodePtr->nid])
              (*(nodePtr->vertexTriangles_))[iter->first[j]-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second + triangleIntervals_[nodePtr->nid-1]);
          }
        }

        for(iter = nodePtr->externalTriangleMap_->begin(); iter != nodePtr->externalTriangleMap_->end(); iter++){
          for(SimplexId j = 0; j < 3; j++){
            if(iter->first.at(j) > vertexIntervals_[nodePtr->nid-1] && iter->first.at(j) <= vertexIntervals_[nodePtr->nid])
              (*(nodePtr->vertexTriangles_))[iter->first.at(j)-vertexIntervals_[nodePtr->nid-1]-1].push_back(iter->second);
          }
        }
        
        return 0;
      }

      /**
       * Get the boundary cells in a given node.
       */
      int getBoundaryCells(ExpandedNode * const nodePtr, const SimplexId dim = 2) const{

        #ifndef TTK_ENABLE_KAMIKAZE
          if(nodePtr->nid <= 0 || nodePtr->nid > nodeNumber_)
            return -1;
        #endif
        
        if(getDimensionality() == 2){
          SimplexId localEdgeNum = edgeIntervals_[nodePtr->nid]-edgeIntervals_[nodePtr->nid-1];
          if(nodePtr->boundaryEdges_ == nullptr){
            nodePtr->boundaryEdges_ = new vector<bool>(localEdgeNum, false);
            if(nodePtr->edgeStars_ == nullptr){
              nodePtr->edgeStars_ = new vector<vector<SimplexId>>();
              getEdgeStars(nodePtr);
            }
            for(SimplexId i = 0; i < localEdgeNum; i++){
              if((*(nodePtr->edgeStars_))[i].size() == 1){
                (*(nodePtr->boundaryEdges_))[i] = true;
              }
            }
          }
          // if boundary vertices are requested
          if(dim == 0 && nodePtr->boundaryVertices_ == nullptr){
            nodePtr->boundaryVertices_ = new vector<bool>(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1], false);
            if(nodePtr->externalEdgeMap_ == nullptr){
              nodePtr->externalEdgeMap_ = new boost::unordered_map<pair_int, SimplexId>();
              buildExternalEdgeMap(nodePtr->nid, nodePtr->externalEdgeMap_);
            }
            // internal edges
            for(auto iter = internalEdgeMaps_[nodePtr->nid].begin(); iter != internalEdgeMaps_[nodePtr->nid].end(); iter++){
              if((*(nodePtr->boundaryEdges_))[iter->second-1]){
                (*(nodePtr->boundaryVertices_))[iter->first.first-vertexIntervals_[nodePtr->nid-1]-1] = true;
                if(iter->first.second <= vertexIntervals_[nodePtr->nid]){
                  (*(nodePtr->boundaryVertices_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1] = true;
                }
              }
            }
            // external edges
            boost::unordered_map<SimplexId, ExpandedNode*> nodeMaps;
            for(auto iter = nodePtr->externalEdgeMap_->begin(); iter != nodePtr->externalEdgeMap_->end(); iter++){
              SimplexId nodeId = vertexIndices_[iter->first.first];
              if(nodeMaps.find(nodeId) == nodeMaps.end()){
                nodeMaps[nodeId] = new ExpandedNode(nodeId);
                getBoundaryCells(nodeMaps[nodeId]);
              }
              if((*(nodeMaps[nodeId]->boundaryEdges_))[iter->second-edgeIntervals_[nodeId-1]-1]){
                (*(nodePtr->boundaryVertices_))[iter->first.second-vertexIntervals_[nodePtr->nid-1]-1] = true;
              }
            }

            // release the memory
            for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++){
              delete iter->second;
            }
          }
        }
        else if(getDimensionality() == 3){
          // get the boundary triangles first
          SimplexId localTriangleNum = triangleIntervals_[nodePtr->nid]-triangleIntervals_[nodePtr->nid-1];
          if(nodePtr->boundaryTriangles_ == nullptr){
            nodePtr->boundaryTriangles_ = new vector<bool>(localTriangleNum, false);
            if(nodePtr->triangleStars_ == nullptr){
              nodePtr->triangleStars_ = new vector<vector<SimplexId>>();
              getTriangleStars(nodePtr);
            }
            for(SimplexId i = 0; i < localTriangleNum; i++){
              if((*(nodePtr->triangleStars_))[i].size() == 1){
                (*(nodePtr->boundaryTriangles_))[i] = true;
              }
            }
          }

          // if the boundary edges are requested
          if(dim == 1 && nodePtr->boundaryEdges_ == nullptr){
            nodePtr->boundaryEdges_ = new vector<bool>(edgeIntervals_[nodePtr->nid]-edgeIntervals_[nodePtr->nid-1], false);
            if(nodePtr->externalTriangleMap_ == nullptr){
              nodePtr->externalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
              buildExternalTriangleMap(nodePtr->nid, nodePtr->externalTriangleMap_);
            }
            // internal triangles
            for(auto iter = internalTriangleMaps_[nodePtr->nid].begin(); iter != internalTriangleMaps_[nodePtr->nid].end(); iter++){
              if((*(nodePtr->boundaryTriangles_))[iter->second-1]){
                pair_int edgePair(iter->first[0], iter->first[1]);
                (*(nodePtr->boundaryEdges_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1] = true;
                edgePair.second = iter->first[2];
                (*(nodePtr->boundaryEdges_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1] = true;
                if(iter->first[1] <= vertexIntervals_[nodePtr->nid]){
                  edgePair.first = iter->first[1];
                  (*(nodePtr->boundaryEdges_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1] = true;
                }
              }
            }
            // external triangles
            boost::unordered_map<SimplexId, ExpandedNode*> nodeMaps;
            for(auto iter = nodePtr->externalTriangleMap_->begin(); iter != nodePtr->externalTriangleMap_->end(); iter++){
              SimplexId nodeId = vertexIndices_[iter->first[0]];
              if(nodeMaps.find(nodeId) == nodeMaps.end()){
                nodeMaps[nodeId] = new ExpandedNode(nodeId);
                getBoundaryCells(nodeMaps[nodeId]);
              }
              if(nodeMaps[nodeId]->boundaryTriangles_->at(iter->second-triangleIntervals_[nodeId-1]-1)){
                if(iter->first[1] > vertexIntervals_[nodePtr->nid-1] && iter->first[1] <= vertexIntervals_[nodePtr->nid]){
                  pair_int edgePair(iter->first[1], iter->first[2]);
                  (*(nodePtr->boundaryEdges_))[internalEdgeMaps_[nodePtr->nid].at(edgePair)-1] = true;
                }
              }
            }

            // release the memory
            for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++){
              delete iter->second;
            }
          }

          // if the boundary vertices are requested
          else if(dim == 0 && nodePtr->boundaryVertices_ == nullptr){
            nodePtr->boundaryVertices_ = new vector<bool>(vertexIntervals_[nodePtr->nid]-vertexIntervals_[nodePtr->nid-1], false);
            if(nodePtr->externalTriangleMap_ == nullptr){
              nodePtr->externalTriangleMap_ = new boost::unordered_map<vector<SimplexId>, SimplexId>();
              buildExternalTriangleMap(nodePtr->nid, nodePtr->externalTriangleMap_);
            }
            // internal triangles
            for(auto iter = internalTriangleMaps_[nodePtr->nid].begin(); iter != internalTriangleMaps_[nodePtr->nid].end(); iter++){
              if((*(nodePtr->boundaryTriangles_))[iter->second-1]){
                for(int j = 0; j < 3; j++){
                  SimplexId vid = iter->first[j];
                  if(vid <= vertexIntervals_[nodePtr->nid]){
                    (*(nodePtr->boundaryVertices_))[vid-vertexIntervals_[nodePtr->nid-1]-1] = true;
                  }
                }
              }
            }
            // external triangles
            boost::unordered_map<SimplexId, ExpandedNode*> nodeMaps;
            for(auto iter = nodePtr->externalTriangleMap_->begin(); iter != nodePtr->externalTriangleMap_->end(); iter++){
              SimplexId nodeId = vertexIndices_[iter->first[0]];
              if(nodeMaps.find(nodeId) == nodeMaps.end()){
                nodeMaps[nodeId] = new ExpandedNode(nodeId);
                getBoundaryCells(nodeMaps[nodeId]);
              }
              if((*(nodeMaps[nodeId]->boundaryTriangles_))[iter->second-triangleIntervals_[nodeId-1]-1]){
                if(iter->first[1] > vertexIntervals_[nodePtr->nid-1] && iter->first[1] <= vertexIntervals_[nodePtr->nid]){
                  (*(nodePtr->boundaryVertices_))[iter->first[1]-vertexIntervals_[nodePtr->nid-1]-1] = true;
                }
                if(iter->first[2] > vertexIntervals_[nodePtr->nid-1] && iter->first[2] <= vertexIntervals_[nodePtr->nid]){
                  (*(nodePtr->boundaryVertices_))[iter->first[2]-vertexIntervals_[nodePtr->nid-1]-1] = true;
                }
              }
            }

            // release the memory
            for(auto iter = nodeMaps.begin(); iter != nodeMaps.end(); iter++){
              delete iter->second;
            }
          }
        }
        else{
          return -1;
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
      vector<boost::unordered_map<pair_int, SimplexId>> internalEdgeMaps_;
      vector<boost::unordered_map<vector<SimplexId>, SimplexId>> internalTriangleMaps_;

      // Cache system
      size_t cacheSize_;
      mutable list<ExpandedNode*> cache_;
      mutable boost::unordered_map<SimplexId, list<ExpandedNode*>::iterator> cacheMap_;

      friend class TestTopoCluster;
  };
}
