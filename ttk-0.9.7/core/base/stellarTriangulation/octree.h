#include <iostream>
#include <assert.h>
#include <unordered_map>

using namespace std;

class OctreeNode
{
    public:
        OctreeNode(){OctreeNode(0);}
        OctreeNode(uint32_t lcode){OctreeNode(lcode, false);}
        OctreeNode(uint32_t lcode, bool hasChild){
            locCode = lcode; 
            childExists = hasChild; 
        }

    protected:
        uint32_t locCode;
        bool childExists;
        // TODO: add necessary fileds, like the vector of vertex ids?

    friend class Octree;
};


class Octree
{
    public:
        Octree(){
            root = OctreeNode(1);
            allNodes[1] = root;
            capacity = 1000; 
        }

        ~Octree(){}

        bool empty() const{
            return root.childExists;
        }

        size_t getNodeTreeDepth(const OctreeNode *node){
            assert(node->locCode);
            for(uint32_t lc=node->locCode, depth=0; lc!=1; lc>>3, depth++)
                return depth;
        }

        OctreeNode* getParentNode(OctreeNode *node){
            const uint32_t locCodeParent = node->locCode >> 3;
            return lookupNode(locCodeParent);
        }

        void visitAll(OctreeNode *node){
            for(int i = 0; i < 8; i++){
                if(node->childExists & (1 << i)){
                    const uint32_t locCodeChild = (node->locCode << 3) | i;
                    const auto child = lookupNode(locCodeChild);
                    visitAll(child);
                }
            }
        }

        // TODO: subdivide function for the octree construction.
        void subdivde(){

        }


    private:
        unordered_map<uint32_t, OctreeNode> allNodes;
        OctreeNode root;
        uint32_t capacity;

        OctreeNode* lookupNode(uint32_t locCode){
            const auto iter = allNodes.find(locCode);
            return (iter == allNodes.end() ? nullptr : &(iter->second));
        }
};
