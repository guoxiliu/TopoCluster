#include <list>
#include <iostream>
#include <AbstractTriangulation.h>
#include <boost/unordered_map.hpp>

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

    friend class FIFOCache;
    friend class LRUCache;
    friend class LFUCache;
    friend class ExplicitTopoCluster;
  };


  /**
   * FIFO cache implementation.
   */  
  class FIFOCache
  {
  private:
    size_t cacheSize_;
    int missCount_;
    list<ExpandedNode*> cache_;
    boost::unordered_map<SimplexId, list<ExpandedNode*>::iterator> cacheMap_;
  public:
    FIFOCache(){
      FIFOCache(100);
    }
    FIFOCache(int size){
      cacheSize_ = size;
      missCount_ = 0;
    }
    ~FIFOCache(){
      for(auto iter = cache_.begin(); iter != cache_.end(); iter++){
        delete *iter;
      }
      cout << "[FIFOCache] Miss count: " << missCount_ << endl;
    }
    ExpandedNode* get(const SimplexId &nodeId){
      // cannot find the expanded node in the cache
      if(cacheMap_.find(nodeId) == cacheMap_.end()){
        missCount_++;
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
  };

  class LRUCache
  {
  private:
    size_t cacheSize_;
    int missCount_;
    list<ExpandedNode*> cache_;
    boost::unordered_map<SimplexId, list<ExpandedNode*>::iterator> cacheMap_;
  public:
    LRUCache(){
      LRUCache(100);
    }
    LRUCache(int size){
      cacheSize_ = size;
      missCount_ = 0;
    }
    ~LRUCache(){
      for(auto iter = cache_.begin(); iter != cache_.end(); iter++){
        delete *iter;
      }
      cout << "[LRUCache] Miss count: " << missCount_ << endl;
    }
    ExpandedNode* get(const SimplexId &nodeId){
      if(cacheMap_.find(nodeId) == cacheMap_.end()){
        missCount_++;
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
  };

  class LFUCache
  {
    struct CacheNode{
      ExpandedNode *value;
      int frequency;
      list<int>::const_iterator iter;
      ~CacheNode() {delete value;}
    };

  private:
    size_t cacheSize_;
    int missCount_;
    int minFrequency_;
    boost::unordered_map<SimplexId, CacheNode> cacheMap_;
    boost::unordered_map<SimplexId, list<int>> freqMap_;

    void touch(CacheNode &node){
      // step 1: update the frequency
      const int prevFreq = node.frequency;
      const int freq = ++(node.frequency);
      // step 2: remove the entry from old frequency list
      freqMap_[prevFreq].erase(node.iter);
      // step 3: remove the frequency list if it is empty
      if(freqMap_[prevFreq].empty() && prevFreq == minFrequency_){
        freqMap_.erase(prevFreq);
        ++minFrequency_;
      }
      // step 4: insert the key into the front of the new frequency list
      freqMap_[freq].push_front(node.value->nid);
      // step 5: update the pointer
      node.iter = freqMap_[freq].cbegin();
    }
  public:
    LFUCache(){
      LFUCache(100);
    }
    LFUCache(int size){
      cacheSize_ = size;
      missCount_ = 0;
      minFrequency_ = 0;
    }
    ~LFUCache(){
      cout << "[LFUCache] Miss count: " << missCount_ << endl;
    }
    ExpandedNode* get(const SimplexId &nodeId){
      auto it = cacheMap_.find(nodeId);
      if(it != cacheMap_.cend()){
        touch(it->second);
        return it->second.value;
      }
      missCount_++;
      if(cacheMap_.size() >= cacheSize_){
        const int keyToRemove = freqMap_[minFrequency_].back();
        freqMap_[minFrequency_].pop_back();
        cacheMap_.erase(keyToRemove);
      }
      const int freq = 1;
      minFrequency_ = freq;
      freqMap_[freq].push_front(nodeId);
      cacheMap_[nodeId] = {new ExpandedNode(nodeId), freq, freqMap_[freq].cbegin()};
      return cacheMap_[nodeId].value;
    }
  };
}

