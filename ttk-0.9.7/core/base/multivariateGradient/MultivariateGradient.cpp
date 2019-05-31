#include                  <MultivariateGradient.h>

using namespace std;
using namespace ttk;

MultivariateGradient::MultivariateGradient(){

  triangulation_ = NULL;
  inputData_ = NULL;
  vertexIndexing = NULL;
  gradient = NULL;
  numberCriticalSimplices = 0;
}

MultivariateGradient::~MultivariateGradient(){

}

void MultivariateGradient::simplexToVertices(Simplex simplex, vector<SimplexId>& vertices){

    vertices = vector<SimplexId>(simplex.first+1);
    for(int i=0; i<simplex.first+1;i++){
        SimplexId vertex;

        switch (simplex.first) {
            case 0:
                vertex = simplex.second;
                break;
            case 1:
                triangulation_->getEdgeVertex(simplex.second,i,vertex);
                break;
            case 2:
                triangulation_->getTriangleVertex(simplex.second,i,vertex);
                break;
            case 3:
                triangulation_->getCellVertex(simplex.second,i,vertex);
                break;
        }

        vertices[i]=vertex;
    }
}


SimplexId MultivariateGradient::getIndexing(Simplex simplex){

    vector<SimplexId> vertices;
    simplexToVertices(simplex,vertices);

    //from the vertices select the field value of the simplex
    SimplexId index = (*vertexIndexing)[vertices[0]];

    for(SimplexId v=1; v<vertices.size(); v++){
        index = index > (*vertexIndexing)[vertices[v]] ? index : (*vertexIndexing)[vertices[v]];
    }

    return index;
}

void MultivariateGradient::getVertexIndexingOfVertices(Simplex simplex, vector<SimplexId>& vec){

    simplexToVertices(simplex,vec);

    for(SimplexId v=0; v<vec.size(); v++){
        vec[v] = (*vertexIndexing)[vec[v]];
    }

}

void MultivariateGradient::computeLowerStar(SimplexId simpl, SimplexesVec& vec){

    vec = SimplexesVec(dimensionality_+1, vector<SimplexId>());

    SimplexId filtrV = getIndexing(Simplex(0,simpl));
    vec[0].push_back(simpl);


    for(int d=1; d<=dimensionality_; d++){

        int number=0;
        switch(d){
            case 1:
                number = triangulation_->getVertexEdgeNumber(simpl);
                break;
            case 2:
                number = triangulation_->getVertexTriangleNumber(simpl);
                break;
            case 3:
                number = triangulation_->getVertexStarNumber(simpl);
                break;
        }


        for(int i=0;i<number;i++){
            SimplexId id=-1;

            switch(d){
            case 1:
                triangulation_->getVertexEdge(simpl,i,id);
                break;
            case 2:
                triangulation_->getVertexTriangle(simpl,i,id);
                break;
            case 3:
                triangulation_->getVertexStar(simpl,i,id);
                break;
            }

            if(getIndexing(Simplex(d,id)) == filtrV)
                vec[d].push_back(id);
        }

    }

}

bool MultivariateGradient::cmpSimplexesLevelSets(const Simplex& s1, const Simplex& s2){

    if(s1.first == s2.first){

        vector<SimplexId> indexing1;
        vector<SimplexId> indexing2;

        getVertexIndexingOfVertices(s1,indexing1);
        getVertexIndexingOfVertices(s2,indexing2);

        sort(indexing1.begin(),indexing1.end(), std::greater<int>());
        sort(indexing2.begin(),indexing2.end(), std::greater<int>());

        for(int i=0; i<indexing1.size(); i++){
            if(indexing1[i] == indexing2[i]){
                continue;
            }

            return indexing1[i] < indexing2[i];
        }
    }


    return s1.first < s2.first;
}


void MultivariateGradient::setPair(const Simplex& tail, const Simplex& head){


    if(tail.first == 0){
        (*gradient)[tail.first][tail.second][0] = head.second;
    }
    else{
       (*gradient)[tail.first][tail.second][1] = head.second;
    }

    (*gradient)[head.first][head.second][0] = tail.second;
}

bool MultivariateGradient::getPair(const Simplex& simpl, Simplex& paired){

    if(simpl.first == 0){
        if((*gradient)[simpl.first][simpl.second][0] != -1){
            paired = Simplex(1,(*gradient)[simpl.first][simpl.second][0]);
            return true;
        }
    }
    else if(simpl.first == dimensionality_){
        if((*gradient)[simpl.first][simpl.second][0] != -1){
            paired = Simplex(dimensionality_-1,(*gradient)[simpl.first][simpl.second][0]);
            return true;
        }
    }
    else{
        if((*gradient)[simpl.first][simpl.second][0] != -1){
            paired = Simplex(simpl.first-1,(*gradient)[simpl.first][simpl.second][0]);
            return true;
        }
        else if((*gradient)[simpl.first][simpl.second][1] != -1){
            paired = Simplex(simpl.first+1,(*gradient)[simpl.first][simpl.second][1]);
            return true;
        }
    }

    return false;
}


bool MultivariateGradient::isCritical(const Simplex& simpl){

    for(auto i : (*gradient)[simpl.first][simpl.second]){
        if(i != -1)
            return false;
    }

    return true;
}

bool MultivariateGradient::isTopCritical(const Simplex& simplex){

    SimplexId simpl = simplex.second;

    for(int d=simplex.first+1; d<=dimensionality_; d++){

        int number=0;
        switch(d){
            case 1:
                number = triangulation_->getVertexEdgeNumber(simpl);
                break;
            case 2:
                if(simplex.first == 0)
                    number = triangulation_->getVertexTriangleNumber(simpl);
                else
                    number = triangulation_->getEdgeTriangleNumber(simpl);
                break;
            case 3:
                if(simplex.first == 0)
                    number = triangulation_->getVertexStarNumber(simpl);
                else if(simplex.first == 1)
                    number = triangulation_->getEdgeStarNumber(simpl);
                else
                    number = triangulation_->getTriangleStarNumber(simpl);
                break;
        }


        for(int i=0;i<number;i++){
            SimplexId id=-1;

            switch(d){
                case 1:
                triangulation_->getVertexEdge(simpl,i,id);
                break;
            case 2:
                if(simplex.first == 0)
                    triangulation_->getVertexTriangle(simpl,i,id);
                else
                    triangulation_->getEdgeTriangle(simpl,i,id);
                break;
            case 3:
                if(simplex.first == 0)
                    triangulation_->getVertexStar(simpl,i,id);
                else if(simplex.first == 1)
                    triangulation_->getEdgeStar(simpl,i,id);
                else
                    triangulation_->getTriangleStar(simpl,i,id);
                break;
            }

            if(isCritical(Simplex(d,id)))
                return false;
        }

    }


    return true;
}

void MultivariateGradient::computeBarycenter(const Simplex& simplex, vector<float>& coords){

    std::vector<SimplexId> vertices;
    simplexToVertices(simplex, vertices);

    coords = vector<float>(3,0);

    for(auto vertex : vertices){
        vector<float> vCoords(3);
        triangulation_->getVertexPoint(vertex, vCoords[0], vCoords[1], vCoords[2]);

        for(int i=0;i<3;i++)
            coords[i] += vCoords[i];
    }

    for(int i=0;i<3;i++)
        coords[i] = coords[i]/static_cast<float>(vertices.size());

}

float MultivariateGradient::computeDistance(const vector<float>& point1, const vector<float>& point2){

    float val=0;
    for(int i=0; i<point1.size(); i++)
        val += pow(point1[i]-point2[i],2);

    return sqrt(val);
}

int MultivariateGradient::computeCriticalClusters(){

    Timer t;

    unordered_map<Simplex, int, boost::hash<Simplex>> simpl_to_label;
    unordered_map<int, list<Simplex>* > label_to_simpl;
    int labels=0;

    std::cout << "Preparing for the critical clusters" << std::endl;

    //Label all critical simplices as belonging to a different cluster
    for(int i=0;i<=dimensionality_;i++){
        for(int j=0; j<(*gradient)[i].size(); j++){
            Simplex simplex = Simplex(i,j);
            if(isCritical(simplex)){
                if(simpl_to_label.find(simplex) == simpl_to_label.end()){

                    //start a new search
                    list<Simplex>* cluster_simplices = new list<Simplex>();
                    cluster_simplices->push_back(simplex);

                    label_to_simpl[labels] = cluster_simplices;
                    simpl_to_label[simplex] = labels++;
                }
            }
        }
    }

    std::cout << "Found " << simpl_to_label.size() << " critical simplices in -" << t.getElapsedTime() << std::endl;





    // computeClusterOfSimplicesDim3(label_to_simpl, simpl_to_label);

    for(int i=0; i<triangulation_->getNumberOfVertices(); i++){
        Simplex simplex(0,i);
        bool first_found=false;
        Simplex lead;

        if(isCritical(simplex)){
            first_found=true;
            lead=simplex;
        }

        for(int j=1; j<=dimensionality_; j++){
            for(int k=0; k<sizeCofacets(simplex,j); k++){
                SimplexId cofacetId = extractCoboundary(simplex, j, k);
                Simplex cofacet(j, cofacetId);

                if(!isCritical(cofacet))
                    continue;

                if(!first_found){
                    lead = cofacet;
                    first_found=true;
                }
                else{

                    int region_s = simpl_to_label[lead];
                    int region_b = simpl_to_label[cofacet];

                    //if the two regions they belong to arre different
                    if(region_s != region_b){
                        //take all simplices with simplex
                        //bring them to be part of the other cluster
                        list<Simplex>* cluster_s = label_to_simpl[region_s];
                        list<Simplex>* cluster_b = label_to_simpl[region_b];

                        if(cluster_s->size() > cluster_b->size()){
                            //bring cluster of boundary to be part of cluster of simplex

                            for(auto iter = cluster_b->begin(); iter != cluster_b->end(); iter++){
                                cluster_s->push_back(*iter);
                                simpl_to_label[*iter]=region_s;
                            }

                            label_to_simpl.erase(region_b);
                            delete cluster_b;

                        }
                        else{
                            //bring cluster of simplex to be part of cluster of boundary
                            for(auto iter = cluster_s->begin(); iter != cluster_s->end(); iter++){
                                cluster_b->push_back(*iter);
                                simpl_to_label[*iter]=region_b;
                            }

                            label_to_simpl.erase(region_s);
                            delete cluster_s;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Union find for the critical clusters - " << t.getElapsedTime() << std::endl;

    std::cout << "Rename labels " << label_to_simpl.size() << std::endl;
    int new_label=0;
    unordered_map<int,int> modified_label;
    for(auto region : label_to_simpl){
        modified_label[region.first] = new_label++;
    }


    //create arcs here only based on descending paths
    vector<ArcSet > list_of_arcs = vector<ArcSet>(label_to_simpl.size(), ArcSet());

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber)
    #endif

    for(int i=0; i<labels; i++){
    // for(auto it = label_to_simpl.begin(); it != label_to_simpl.end(); it++){
        if(label_to_simpl.find(i) == label_to_simpl.end())
            continue;

        int label = modified_label[i];

        ArcSet arcs;
        list<Simplex>* cSimplices = label_to_simpl[i];
        for(list<Simplex>::iterator it = cSimplices->begin(); it != cSimplices->end(); it++){
            if(it->first > 0){
                CMap critical_reached;
                unordered_set<int> geometry;

                extractOutgoingPath(*it, critical_reached, false, geometry);

                for(auto reached : critical_reached){
                    Simplex reachedSimplex(it->first-1,reached.first);
                    //get simplex real index
                    // assert(simpl_to_label.find(reachedSimplex) != simpl_to_label.end());
                    int other_label = modified_label[simpl_to_label[reachedSimplex]];
                    if(other_label != label)
                        arcs.insert(ArcLenght(other_label, reached.second));
                }
            }
        }
        list_of_arcs[modified_label[i]]=arcs;
    }
    std::cout << "Created arcs in - " << t.getElapsedTime() << std::endl;

    std::cout << "Selecting bestrepresentative" << std::endl;
    vector<int > index_to_best = vector<int >(label_to_simpl.size(), 0);
    for(auto region : label_to_simpl){
        // list<Simplex>* cSimplices = region.second;

        // vector<float> distances = vector<float>(cSimplices->size(),0);
        // list<Simplex>::iterator it = cSimplices->begin();
        // for(int i=0;i<cSimplices->size(); i++){
        //     list<Simplex>::iterator it2 = it;
        //     it2++;

        //     for(int j=i+1;j<cSimplices->size(); j++){

        //         vector<float> coords1;
        //         computeBarycenter(*it, coords1);

        //         vector<float> coords2;
        //         computeBarycenter(*it2, coords2);

        //         float distance = computeDistance(coords1,coords2);

        //         distances[i] += distance;
        //         distances[j] += distance;

        //         it2++;
        //     }

        //     it++;
        // }
        // index_to_best[modified_label[region.first]] = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));

        index_to_best[modified_label[region.first]] = 0;
    }

    cout << "Creating cluster graph" << endl;
    cclusters = new CriticalClusterGraph();
    cclusters->addClusters(label_to_simpl, modified_label, list_of_arcs, index_to_best);

    cout << "Critical clusters identified " << cclusters->size() << t.getElapsedTime() << endl;
}


void MultivariateGradient::extractOutgoingPath(Simplex simpl, CMap& critical_reached,
                                                     bool saveGeometry,
                                                     unordered_set<int>& geometry){

    //like descending

    stack<pair<SimplexId,uint> > q;

    unordered_set<int> visited;

    SimplexId dim_simpl = simpl.first;
    q.push(pair<SimplexId,uint>(simpl.second,0));

    while(!q.empty()){
        pair<SimplexId,uint> next = q.top();
        q.pop();

        for(int i=0; i<dim_simpl+1; i++){

            int facet = extractBoundary(Simplex(dim_simpl,next.first), dim_simpl-1, i);

            Simplex tail(dim_simpl-1,facet);

            Simplex arrow;
            if( getPair(tail, arrow)){

                if(arrow.first > tail.first && next.first != arrow.second){

                    if(visited.find(arrow.second) == visited.end()){
                        visited.insert(arrow.second);
                        q.push(pair<SimplexId,uint>(arrow.second,next.second+1));
                        if(saveGeometry)
                            geometry.insert(arrow.second);
                    }
                }
                else{
                    //here frontier if we want
                }
            }
            else{
                auto inside = critical_reached.find(facet);
                if(inside == critical_reached.end())
                    critical_reached[facet] = next.second;
                else{
                    if(inside->second > next.second)
                        inside->second = next.second;
                }
            }
        }
    }
}

void MultivariateGradient::extractIngoingPath(Simplex simpl, unordered_set<int>& critical_reached,
                                                    bool saveGeomettry,
                                                    unordered_set<int>& geometry){
    assert(simpl.first < dimensionality_);
    //like ascending

    queue<SimplexId> q;

    SimplexId dim_simpl = simpl.first;
    q.push(simpl.second);

    while(!q.empty()){
        SimplexId next = q.front();
        q.pop();

        for(int i=0; i<sizeCofacets(Simplex(dim_simpl,next),dim_simpl+1); i++){

            int cofacet = extractCoboundary(Simplex(dim_simpl,next), dim_simpl+1, i);
            Simplex arrow(dim_simpl+1,cofacet);
            Simplex tail;
            if( getPair(arrow,tail)){

                if(arrow.first > tail.first && next != tail.second){
                    if(geometry.find(arrow.second) == geometry.end() ){
                      q.push(tail.second);
                      geometry.insert(arrow.second);
                    }

                }
                else{
                    //add something in frontier here when you are saving the geometry
                }
            }
            else{
                //is critical
                critical_reached.insert(cofacet);
            }
        }
    }

}


int MultivariateGradient::extractBoundary(Simplex simpl,int dimension_b, int index_b){

    SimplexId id=-1;
    switch(simpl.first){
        case 0:
            return id;
        case 1:
            triangulation_->getEdgeVertex(simpl.second,index_b,id);
            break;
        case 2:
            if(dimension_b==0)
                triangulation_->getTriangleVertex(simpl.second, index_b, id);
            else{
                triangulation_->getTriangleEdge(simpl.second, index_b, id);
            }
            break;
        case 3:
            if(dimension_b==0)
                triangulation_->getCellVertex(simpl.second,index_b,id);
            else if(dimension_b==0)
                triangulation_->getCellEdge(simpl.second,index_b,id);
            else
                triangulation_->getCellTriangle(simpl.second,index_b,id);
            break;
    }

    return id;
}


int MultivariateGradient::sizeCofacets(Simplex simpl, int dimension_cb){

    switch(simpl.first){
        case 0:{
            if(dimension_cb==1)
                return triangulation_->getVertexEdgeNumber(simpl.second);
            else if(dimension_cb==2)
                return triangulation_->getVertexTriangleNumber(simpl.second);
            else
                return triangulation_->getVertexStarNumber(simpl.second);
        }

        case 1:{
            if(dimension_cb==2)
                return triangulation_->getEdgeTriangleNumber(simpl.second);
            else
                return triangulation_->getEdgeStarNumber(simpl.second);
        }


        case 2:
            return triangulation_->getTriangleStarNumber(simpl.second);
        case 3:
            return 0;
    }

}

int MultivariateGradient::extractCoboundary(Simplex simpl, int dimension_cb, int index_cb){

    SimplexId id=-1;

    switch(simpl.first){
        case 0:
            if(dimension_cb==1)
                triangulation_->getVertexEdge(simpl.second,index_cb,id);
            else if(dimension_cb==2)
                triangulation_->getVertexTriangle(simpl.second,index_cb,id);
            else
                triangulation_->getVertexStar(simpl.second,index_cb,id);
            break;

        case 1:
            if(dimension_cb==2)
                triangulation_->getEdgeTriangle(simpl.second,index_cb,id);
            else
                triangulation_->getEdgeStar(simpl.second,index_cb,id);
            break;

        case 2:
            triangulation_->getTriangleStar(simpl.second,index_cb,id);
            break;
        case 3:
            break;
    }

    return id;
}


void MultivariateGradient::simplifyGraph(int minLenght){
    // cclusters->printSimplices();
    cout << "Simplify path smaller than " << minLenght << endl;
    for( int i=0; i<cclusters->size(); i++){
        if(!cclusters->isDead(i)){
            unordered_map<int,int> arcs;
            cclusters->getConnections(i,arcs,true);
            for(auto arc : arcs){
                if(arc.second <= minLenght){
                    bool newArc = cclusters->contractArc(i,arc.first,arc.second,minLenght);
                    if(cclusters->isDead(i))
                        break;
                }
            }
        }
    }

    // cclusters->printArcs();
    // cclusters->printSimplices();
}


float MultivariateGradient::simplifyClusters(int maxSize){

    int max = cclusters->getClusterSize(0);
    int min = cclusters->getClusterSize(0);
    for( int i=0; i<cclusters->size(); i++){
        int val = cclusters->getClusterSize(i);
        if(val > max)
            max = val;
        if(val < min)
            min = val;
    }

    float val = float(float(maxSize)/100.0*float(max-min)) + float(min);

    for( int i=0; i<cclusters->size(); i++){
        if(cclusters->getClusterSize(i) < val){
            cclusters->removeCluster(i);
        }
    }

    return val;
}


void MultivariateGradient::mergeRegions(int region_s, int region_b, unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    if(region_s != region_b){
        //take all simplices with simplex
        //bring them to be part of the other cluster
        list<Simplex>* cluster_s = label_to_simpl[region_s];
        list<Simplex>* cluster_b = label_to_simpl[region_b];

        if(cluster_s->size() > cluster_b->size()){
            //bring cluster of boundary to be part of cluster of simplex

            for(auto iter = cluster_b->begin(); iter != cluster_b->end(); iter++){
                cluster_s->push_back(*iter);
                simpl_to_label[*iter]=region_s;
            }

            label_to_simpl.erase(region_b);
            delete cluster_b;

        }
        else{
            //bring cluster of simplex to be part of cluster of boundary
            for(auto iter = cluster_s->begin(); iter != cluster_s->end(); iter++){
                cluster_b->push_back(*iter);
                simpl_to_label[*iter]=region_b;
            }

            label_to_simpl.erase(region_s);
            delete cluster_s;
        }
    }

}

void MultivariateGradient::unionFindEdgeVertices(unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    for(int i=0; i<triangulation_->getNumberOfEdges(); i++){
        Simplex simplex(1,i);

        if(!isCritical(simplex)){
            continue;
        }

        for(int k=0; k<2; k++){
            SimplexId facetId = extractBoundary(simplex, 0, k);
            Simplex facet(0, facetId);

            if(!isCritical(facet))
                continue;

            int region_s = simpl_to_label[simplex];
            int region_b = simpl_to_label[facet];
            mergeRegions(region_s,region_b,simpl_to_label,label_to_simpl);
        }
    }
}

void MultivariateGradient::unionFindTriangleTetras(unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    for(int i=0; i<triangulation_->getNumberOfCells(); i++){
        Simplex simplex(3,i);

        if(!isCritical(simplex)){
            continue;
        }

        for(int k=0; k<4; k++){
            cout << "Got here" << endl;
            SimplexId facetId = extractBoundary(simplex, 2, k);
            cout << "FACE" << endl;
            Simplex facet(2, facetId);

            if(!isCritical(facet))
                continue;

            cout << "Got here" << endl;
            int region_s = simpl_to_label[simplex];
            int region_b = simpl_to_label[facet];
            cout << "Got here?" << endl;
            //if the two regions they belong to arre different
            mergeRegions(region_s,region_b,simpl_to_label,label_to_simpl);
            cout << "Got here" << endl;
        }
    }
}

void MultivariateGradient::unionFindTetraEdges(unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    for(int i=0; i<triangulation_->getNumberOfCells(); i++){
        Simplex simplex(3,i);

        if(!isCritical(simplex)){
            continue;
        }

        for(int k=0; k<6; k++){
            SimplexId facetId = extractBoundary(simplex, 1, k);
            Simplex facet(1, facetId);

            if(!isCritical(facet))
                continue;

            int region_s = simpl_to_label[simplex];
            int region_b = simpl_to_label[facet];

            //if the two regions they belong to arre different
            mergeRegions(region_s,region_b,simpl_to_label,label_to_simpl);
        }
    }
}

void MultivariateGradient::unionFindTriangleVertices(unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    for(int i=0; i<triangulation_->getNumberOfTriangles(); i++){
        Simplex simplex(2,i);

        if(!isCritical(simplex)){
            continue;
        }

        for(int k=0; k<3; k++){
            SimplexId facetId = extractBoundary(simplex, 0, k);
            Simplex facet(0, facetId);

            if(!isCritical(facet))
                continue;

            int region_s = simpl_to_label[simplex];
            int region_b = simpl_to_label[facet];

            //if the two regions they belong to arre different
            mergeRegions(region_s,region_b,simpl_to_label,label_to_simpl);
        }
    }
}

void MultivariateGradient::unionFindTetraVertices(unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    for(int i=0; i<triangulation_->getNumberOfCells(); i++){
        Simplex simplex(3,i);

        if(!isCritical(simplex)){
            continue;
        }

        for(int k=0; k<4; k++){
            SimplexId facetId = extractBoundary(simplex, 0, k);
            Simplex facet(0, facetId);

            if(!isCritical(facet))
                continue;

            int region_s = simpl_to_label[simplex];
            int region_b = simpl_to_label[facet];

            //if the two regions they belong to arre different
            mergeRegions(region_s,region_b,simpl_to_label,label_to_simpl);
        }
    }
}

void MultivariateGradient::unionFindTriangleEdges(unordered_map<Simplex, int, boost::hash<Simplex> >& simpl_to_label, unordered_map<int, list<Simplex>* >& label_to_simpl){

    for(int i=0; i<triangulation_->getNumberOfTriangles(); i++){
        Simplex simplex(2,i);

        if(!isCritical(simplex)){
            continue;
        }

        for(int k=0; k<3; k++){
            SimplexId facetId = extractBoundary(simplex, 1, k);
            Simplex facet(1, facetId);

            if(!isCritical(facet))
                continue;

            int region_s = simpl_to_label[simplex];
            int region_b = simpl_to_label[facet];

            //if the two regions they belong to arre different
            mergeRegions(region_s,region_b,simpl_to_label,label_to_simpl);
        }
    }
}


int MultivariateGradient::computeCriticalClustersAlt(){

    Timer t;

    vector< list<Simplex>* > vertex_to_simpl(triangulation_->getNumberOfVertices(),NULL);

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber)
    #endif
    for(int i=0; i<triangulation_->getNumberOfVertices(); i++){
        Simplex simplex(0,i);

        vertex_to_simpl[i] = new list<Simplex>();
        if(isCritical(simplex))
            vertex_to_simpl[i]->push_back(simplex);

        for(int j=1; j<=dimensionality_;j++){
            for(int k=0; k<sizeCofacets(simplex,j); k++){
                SimplexId cofacetId = extractCoboundary(simplex, j, k);
                Simplex cofacet(j, cofacetId);

                if(isCritical(cofacet))
                    vertex_to_simpl[i]->push_back(cofacet);
            }
        }
    }

    unordered_map<Simplex, int, boost::hash<Simplex>> simpl_to_label;
    unordered_map<int, list<Simplex>* > label_to_simpl;

    int label=0;
    vector<bool > visited(triangulation_->getNumberOfVertices(), false);
    for(int i=0; i<triangulation_->getNumberOfVertices(); i++){
        if(!visited[i] && vertex_to_simpl[i]->size() > 0){

            set<int> touched;
            touched.insert(i);

            visited[i]=true;

            queue<int> q;
            q.push(i);

            while(!q.empty()){
                int v = q.front();
                q.pop();

                //collect other vertices;
                unordered_set<int> adjacent_v;
                collectVertices(vertex_to_simpl[v],adjacent_v);

                //if they are not visited continue;
                for(auto v : adjacent_v){
                    if(!visited[v]){
                        q.push(v);
                        touched.insert(v);
                        visited[v]=true;
                    }
                }
            }

            //collect all the vertices touched this time and return those items as a cluster
            unordered_set<Simplex, boost::hash<Simplex> > simplices_in_cluster;
            for(auto v : touched){
                simplices_in_cluster.insert(vertex_to_simpl[v]->begin(), vertex_to_simpl[v]->end());
            }


            for(auto s : simplices_in_cluster){
                simpl_to_label[s] = label;
            }

            label_to_simpl[label++] = new list<Simplex>(simplices_in_cluster.begin(), simplices_in_cluster.end());
        }
    }

    std::cout << "Union find for the critical clusters - " << t.getElapsedTime() << std::endl;


    //create arcs here only based on descending paths
    vector<ArcSet > list_of_arcs = vector<ArcSet>(label_to_simpl.size(), ArcSet());

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber)
    #endif
    for(int i=0; i<label; i++){
    // for(auto it = label_to_simpl.begin(); it != label_to_simpl.end(); it++){
        if(label_to_simpl.find(i) == label_to_simpl.end())
            continue;

        int label = i;

        ArcSet arcs;
        list<Simplex>* cSimplices = label_to_simpl[i];
        for(list<Simplex>::iterator it = cSimplices->begin(); it != cSimplices->end(); it++){
            if(it->first > 0){
                CMap critical_reached;
                unordered_set<int> geometry;

                extractOutgoingPath(*it, critical_reached, false, geometry);

                for(auto reached : critical_reached){
                    Simplex reachedSimplex(it->first-1,reached.first);
                    //get simplex real index
                    // assert(simpl_to_label.find(reachedSimplex) != simpl_to_label.end());
                    int other_label = simpl_to_label[reachedSimplex];
                    if(other_label != label)
                        arcs.insert(ArcLenght(other_label, reached.second));
                }
            }
        }
        list_of_arcs[i]=arcs;
    }
    std::cout << "Created arcs in - " << t.getElapsedTime() << std::endl;

    std::cout << "Selecting bestrepresentative" << std::endl;

    vector<int > index_to_best = vector<int >(label_to_simpl.size(), 0);
    for(auto region : label_to_simpl){
        // list<Simplex>* cSimplices = region.second;

        // vector<float> distances = vector<float>(cSimplices->size(),0);
        // list<Simplex>::iterator it = cSimplices->begin();
        // for(int i=0;i<cSimplices->size(); i++){
        //     list<Simplex>::iterator it2 = it;
        //     it2++;

        //     for(int j=i+1;j<cSimplices->size(); j++){

        //         vector<float> coords1;
        //         computeBarycenter(*it, coords1);

        //         vector<float> coords2;
        //         computeBarycenter(*it2, coords2);

        //         float distance = computeDistance(coords1,coords2);

        //         distances[i] += distance;
        //         distances[j] += distance;

        //         it2++;
        //     }

        //     it++;
        // }
        // index_to_best[modified_label[region.first]] = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));

        index_to_best[region.first] = 0;
    }

    cout << "Creating cluster graph" << endl;
    cclusters = new CriticalClusterGraph();
    cclusters->addClusters(label_to_simpl, list_of_arcs, index_to_best);
    backup = new CriticalClusterGraph();
    *backup = *cclusters;
    cout << "Critical clusters identified " << cclusters->size() << " - " << t.getElapsedTime() << endl;
}

void MultivariateGradient::collectVertices(list<Simplex>* star, unordered_set<int>& adjacent_v){

    for(auto it = star->begin(); it != star->end(); it++){
        vector<SimplexId> vertices;
        simplexToVertices(*it, vertices);
        adjacent_v.insert(vertices.begin(), vertices.end());
    }
}
