

#include <Triangulation.h>
#include <Wrapper.h>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <queue>


using namespace std;
using namespace ttk;

typedef std::pair<int,int> ArcLenght;
typedef std::unordered_set<ArcLenght, boost::hash<ArcLenght > > ArcSet;

typedef std::pair<short int,SimplexId> Simplex;


class Cluster{

    vector<Simplex > cells;
    unordered_map<int, int> outgoing_arcs;
    unordered_map<int, int> ingoing_arcs;
    
    unordered_set<int> merged_clusters;
    
    int best_draw;
    char type;

public:
    // inline void addArcs(ArcSet& arcs_in){arcs = vector<ArcLenght>(arcs_in.begin(), arcs_in.end());}
    // inline void addArcs(vector<ArcLenght>& arcs_in){arcs.insert(arcs.end(), arcs_in.begin(), arcs_in.end());}
    // inline vector<ArcLenght>& getArcs(){ return arcs;}

    Cluster(){
        merged_clusters = unordered_set<int>();
        outgoing_arcs= unordered_map<int, int>();
        ingoing_arcs= unordered_map<int, int>();
    }

    inline bool addArc(int label, int lenght, bool outgoing){

        if(outgoing){
            unordered_map<int, int>::iterator it = outgoing_arcs.find(label);
            if(it == outgoing_arcs.end()){
                outgoing_arcs[label] = lenght;
                return false;
            }
            // else
            // {
            //     it->second = lenght;
            // }
        }
        else{
            unordered_map<int, int>::iterator it = ingoing_arcs.find(label);
            if(it == ingoing_arcs.end()){
                ingoing_arcs[label] = lenght;
                return false;
            }
            // else
            // {
            //     if(it->second > lenght)
            //         it->second = lenght;
            // }            
        }

        return true;
    }

    inline void getArcs(unordered_map<int, int>& arcs_out, bool outgoing){
        if(outgoing)
            arcs_out = outgoing_arcs;
        else
            arcs_out = ingoing_arcs;
        
    }

    inline void removeArc(int arc){
        outgoing_arcs.erase(arc);
        ingoing_arcs.erase(arc);
    }

    inline void addMergedCluster(int cluster){
        merged_clusters.insert(cluster);
    }

    inline void getMergedClusters(unordered_set<int>& clusters){
        clusters = merged_clusters;
    }
        

    inline void addBest(int best_draw_in){best_draw=best_draw_in;}
    inline Simplex* getBestDraw(){ 
        return &cells[best_draw];}

    inline void addCells(list<Simplex>* inCells){cells = vector<Simplex>(inCells->begin(),inCells->end());}
    inline void addCells(vector<Simplex>& inCells){cells.insert(cells.end(), inCells.begin(), inCells.end());}
    inline void getSimplices(vector<Simplex>& simpl){simpl = cells;}
    inline int getClusterSize(){return cells.size();}

    inline char getClusterType(){
        if(ingoing_arcs.size() > 0 && outgoing_arcs.size() > 0)
            return 'a';
        else if(ingoing_arcs.size() > 0)
            return 'b';
        else
            return 'c';
    }
};


class CriticalClusterGraph{

    vector<Cluster> clusters;
    vector<bool> isAlive;

    vector<int> clusterIndexBySize;

public:
    inline CriticalClusterGraph(){};

    inline int size(){ return clusters.size();}

    inline int getClusterSize(int i){return clusters[i].getClusterSize();}

    inline void getClusterBestDraw(int i, Simplex& simplex){
        //compute a representative coordinate for each cluster
        simplex = *clusters[i].getBestDraw();
    }

    inline void getConnections(int i, unordered_map<int, int>& connected, bool outgoing){
        //compute arcs for such cluster
        clusters[i].getArcs(connected, outgoing);
    }

    inline void getMergedClusters(int cluster, unordered_set<int>& ret_clusters){
        clusters[cluster].getMergedClusters(ret_clusters);
    }

    inline char getType(int i){
        //compute cluster's type
        return clusters[i].getClusterType();
    }

    inline void addClusters(unordered_map<int, list<Simplex>* > label_to_simpl, unordered_map<int,int>& modified_label, vector<ArcSet >& arcs, vector<int>& index_to_best){

        clusters = vector<Cluster>(label_to_simpl.size());
        isAlive = vector<bool>(label_to_simpl.size(),true);

        int iter=0;
        for(auto region : label_to_simpl){
            iter = modified_label[region.first];
            clusters[iter].addBest(index_to_best[iter]);
            clusters[iter].addCells(region.second);

            for(auto arc : arcs[iter]){
                clusters[iter].addArc(arc.first,arc.second,true);
                clusters[arc.first].addArc(iter, arc.second,false);
            }

            delete region.second;
        }
    }

    inline void addClusters(unordered_map<int, list<Simplex>* > label_to_simpl, vector<ArcSet >& arcs, vector<int>& index_to_best){

        clusters = vector<Cluster>(label_to_simpl.size());
        isAlive = vector<bool>(label_to_simpl.size(),true);

        int iter=0;
        for(auto region : label_to_simpl){
            iter = region.first;
            clusters[iter].addBest(index_to_best[iter]);
            clusters[iter].addCells(region.second);

            for(auto arc : arcs[iter]){
                clusters[iter].addArc(arc.first,arc.second,true);
                clusters[arc.first].addArc(iter, arc.second,false);
            }

            delete region.second;
        }
    }

    inline void getSimplices(int i, vector<Simplex>& simplices){
        clusters[i].getSimplices(simplices);
    }

    inline bool contractArc(int i_cluster1, int i_cluster2, int length, int minLength){
        vector<Simplex> simplices_c1;
        vector<Simplex> simplices_c2;
        
        clusters[i_cluster1].getSimplices(simplices_c1);
        clusters[i_cluster2].getSimplices(simplices_c2);

        vector<vector<int> > new_arcs; 
        bool newArc = false;

        if(simplices_c1.size() > simplices_c2.size()){
            // cout << "REMOVE " << i_cluster2 << " in " << i_cluster1 << endl;
            clusters[i_cluster1].addCells(simplices_c2);

            unordered_map<int, int> arcs;
            clusters[i_cluster2].getArcs(arcs,false);

            for(auto arc : arcs){
                removeArc(arc.first, i_cluster2);
                if(arc.first == i_cluster1)
                    continue;
                if(!clusters[i_cluster1].addArc(arc.first,arc.second+length,false)){
                    clusters[arc.first].addArc(i_cluster1,arc.second+length,true);
                    newArc = true;
                    if(arc.second+length <= minLength){
                        vector<int> newarc = {i_cluster1,arc.first,arc.second+length};
                        new_arcs.push_back(newarc);
                    }
                }

            }

            clusters[i_cluster2].getArcs(arcs, true);
            for(auto arc : arcs){
                removeArc(arc.first, i_cluster2);
                if(arc.first == i_cluster1)
                    continue;
                if(!clusters[i_cluster1].addArc(arc.first,arc.second+length,true)){
                    clusters[arc.first].addArc(i_cluster1,arc.second+length,false);
                    newArc = true;
                    if(arc.second+length <= minLength){
                        vector<int> newarc = {i_cluster1,arc.first,arc.second+length};
                        new_arcs.push_back(newarc);
                    }
                }
            }

            clusters[i_cluster2] = Cluster();
            isAlive[i_cluster2] = false;
        }
        else{
            // cout << "REMOVE " << i_cluster1 << " in " << i_cluster2 << endl;
            clusters[i_cluster2].addCells(simplices_c1);

            unordered_map<int, int> arcs;
            clusters[i_cluster1].getArcs(arcs, false);

            for(auto arc : arcs){
                removeArc(arc.first, i_cluster1);
                if(arc.first == i_cluster2)
                    continue;
                if (!clusters[i_cluster2].addArc(arc.first,arc.second+length,false)){
                    clusters[arc.first].addArc(i_cluster2,arc.second+length,true);
                    newArc = true;
                    if(arc.second+length <= minLength){
                        vector<int> newarc = {i_cluster2,arc.first,arc.second+length};
                        new_arcs.push_back(newarc);
                    }
                }
                

            }

            clusters[i_cluster1].getArcs(arcs, true);
            for(auto arc : arcs){
                removeArc(arc.first, i_cluster1);
                if(arc.first == i_cluster2)
                    continue;
                if(!clusters[i_cluster2].addArc(arc.first,arc.second+length,true)){
                    clusters[arc.first].addArc(i_cluster2,arc.second+length,false);
                    newArc = true;
                    if(arc.second+length <= minLength){
                        vector<int> newarc = {i_cluster2,arc.first,arc.second+length};
                        new_arcs.push_back(newarc);
                    }
                }
            }


            clusters[i_cluster1] = Cluster();
            isAlive[i_cluster1] = false;
        }
    
        for(auto arc : new_arcs){
            if(!isDead(arc[0]) && !isDead(arc[1]))
                contractArc(arc[0],arc[1],arc[2],minLength);
        }

        return newArc;
    }

    inline void removeCluster(int i){
        unordered_map<int, int> outgoing_arcs;
        unordered_map<int, int> incoming_arcs;
        clusters[i].getArcs(outgoing_arcs,true);
        clusters[i].getArcs(incoming_arcs,false);

        for(auto arc : incoming_arcs){
            for(auto old_arc : outgoing_arcs){
                clusters[arc.first].addArc(old_arc.first, arc.second+old_arc.second, true);
                clusters[old_arc.first].addArc(arc.first, arc.second+old_arc.second, false);
            }

            clusters[arc.first].removeArc(i);
            clusters[arc.first].addMergedCluster(i);
        }

        for(auto old_arc : outgoing_arcs)
            clusters[old_arc.first].removeArc(i);

        outgoing_arcs.clear();
        incoming_arcs.clear();
 
        isAlive[i] = false;

    }

    inline void removeArc(int from, int to){
        clusters[from].removeArc(to);
        clusters[to].removeArc(from);
    }

    inline bool isDead(int i){ return !isAlive[i]; }

    inline void printSimplices(){
        for(int i=0; i<clusters.size(); i++){
            vector<Simplex> arcs;
            clusters[i].getSimplices(arcs);
            
            cout << "Node " << i << " ";
            for(auto simpl : arcs){
                cout << "(" << simpl.first << "," << simpl.second << ") "; 
            }
            cout << endl;        
        }
        cout << endl;
    }

    inline void sort(){

        auto cmp = [](pair<int,int> left, pair<int,int> right) { return left.first < right.first;};
        std::priority_queue<pair<int,int>, vector<pair<int,int>>, decltype(cmp)> q(cmp);

        for(int i=0; i<clusters.size(); i++){
            if( isAlive[i])
                q.push(pair<int,int>(this->getClusterSize(i),i));
        }
        
        int step=0;
        clusterIndexBySize = vector<int>(101,-1);
        while(step < 100 && !q.empty()){
            clusterIndexBySize[step++] = q.top().second;
            q.pop();
        }
    }

    inline int getSortedClusterIndex(int localIndex){ return clusterIndexBySize[localIndex];}
};
