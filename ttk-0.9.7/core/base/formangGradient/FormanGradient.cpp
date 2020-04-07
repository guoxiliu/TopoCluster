#include                  <FormanGradient.h>

using namespace std;
using namespace ttk;

FormanGradient::FormanGradient(){

  indexing_ = NULL;
  triangulation_ = NULL;
}

FormanGradient::~FormanGradient(){
  
}


void FormanGradient::computeGradient(vector<int>* indexing){

    indexing_ = indexing;

    Memory m;
    //prepare gradient
    gradient_ = new vector<boost::dynamic_bitset<>* >(dimensionality_);

    (*gradient_)[0] = new boost::dynamic_bitset<>(triangulation_->getNumberOfEdges()*2);

    if(dimensionality_ >= 2)
        (*gradient_)[1] = new boost::dynamic_bitset<>(triangulation_->getNumberOfTriangles()*3);

    if(dimensionality_ == 3)
        (*gradient_)[2] = new boost::dynamic_bitset<>(triangulation_->getNumberOfCells()*4);

    criticalSimplices = vector<set<SimplexId> >(dimensionality_+1, set<SimplexId>());

    //compute indexing on the data
    SimplexId vertexNumber = triangulation_->getNumberOfVertices();

    cout << "Current memory usage with gradient" << m.getInstantUsage() << endl;
    cout << "Start computing homotopy expansion " << endl;

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for default(shared)
    #endif
    for(SimplexId i=0; i<vertexNumber; i++){
        //compute homotopy expansion
        int local = homotopyExpansion(i);
    }
}


int FormanGradient::getIndex(Simplex simplex){

    //from the simplex extract its boundary vertices
    vector<SimplexId> vertices;
    simplexToVertices(simplex,vertices);

    int index_val = (*indexing_)[vertices[0]];

    for(SimplexId v=1; v<vertices.size(); v++){
      index_val = index_val > (*indexing_)[vertices[v]] ? index_val : (*indexing_)[vertices[v]];
    }

    return index_val;
}

void FormanGradient::simplexToVertices(Simplex simplex, vector<SimplexId>& vertices){

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

void FormanGradient::computeBarycenter(Simplex& simplex, vector<float>& coords){

    coords = vector<float>(3,0);
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

        float x,y,z;
        triangulation_->getVertexPoint(vertex, x,y,z);
        coords[0] += x;
        coords[1] += y;
        coords[2] += z;

    }

    for(int i=0;i<3;i++)
        coords[i] /= (float)simplex.first+1;
}



int FormanGradient::homotopyExpansion(int v){

    auto lexycographic = bind(&FormanGradient::cmpSimplexesLevelSets, this, _1, _2);
    //organize the lower star
    Simplex vertex = Simplex(0,v);

    SimplexesSet order = SimplexesSet(lexycographic);
    for(int j=0; j<sizeCofacets(vertex,1); j++){
        int coface = extractCoboundary(vertex,1,j);
        int index = getIndex(Simplex(1,coface));
        if(index == (*indexing_)[v]){
            order.insert(Simplex(1,coface));
        }
    }

    if(order.size() == 0) {
        criticalSimplices[0].insert(v);
        return 0;
    }

    else{

        Simplex top = *order.begin();
        order.erase(top);
        setPair(vertex, top);

        map<SimplexId, int> facesName;
        vector<SimplexId> facesId(order.size(),0);

        int index=0;
        for(auto simpl : order){
            facesName[simpl.second] = index;
            facesId[index++] = simpl.second;
        }


        order.clear();
        order = SimplexesSet(lexycographic);
        for(int j=0; j<sizeCofacets(vertex,2); j++){
            int coface = extractCoboundary(vertex,2,j);
            int index = getIndex(Simplex(2,coface));
            if(index == (*indexing_)[v]){
                order.insert(Simplex(2,coface));
            }
        }

        map<SimplexId, int> cofacesName;
        vector<SimplexId> cofacesId(order.size(),0);

        index=0;
        for(auto simpl : order){
            cofacesName[simpl.second] = index;
            cofacesId[index++] = simpl.second;
        }

        //now check for pairings starting from the triangles
        while(!facesName.empty()){

            bool someonePaired = false;
            for(int i=0; i<cofacesId.size(); i++){

                if(cofacesId[i] >= 0){
                    SimplexId pairable;
                    int num = numPairableFaces( Simplex(2,cofacesId[i]), facesName, pairable);

                    if(num == 1){
                        setPair( Simplex(1,pairable), Simplex(2,cofacesId[i]));

                        someonePaired = true;

                        cofacesName.erase(cofacesId[i]);
                        cofacesId[i] = -2;

                        facesId[facesName[pairable]] = -2;
                        facesName.erase(pairable);
                    }
                }
            }

            if(!someonePaired){
                for(int i=0; i<facesId.size(); i++){
                    if(facesId[i] >= 0){
                        criticalSimplices[1].insert(facesId[i]);
                        facesName.erase(facesId[i]);
                        facesId[i] = -2;
                     }
                }
            }
        }
        if(dimensionality_ == 2){
            for(int i=0; i<cofacesId.size(); i++){
                if(cofacesId[i] >= 0){
                    criticalSimplices[2].insert(cofacesId[i]);
                }
            }
        }
        else if(dimensionality_ == 3){

//            cout << "Higher dimension" << endl;

            facesId = cofacesId;
            facesName = cofacesName;

            order.clear();
            order = SimplexesSet(lexycographic);
            for(int j=0; j<sizeCofacets(vertex,3); j++){
                int coface = extractCoboundary(vertex,3,j);
                int index = getIndex(Simplex(3,coface));
                if(index == (*indexing_)[v]){
                    order.insert(Simplex(3,coface));
                }
            }

            cofacesName.clear();
            vector<SimplexId> cofacesId(order.size(),0);

            index=0;
            for(auto simpl : order){
                cofacesName[simpl.second] = index;
                cofacesId[index++] = simpl.second;
            }

//            cout << "Populated ready to start" << endl;

            while(!facesName.empty()){

                bool someonePaired = false;
                for(int i=0; i<cofacesId.size(); i++){

                    if(cofacesId[i] >= 0){
                        SimplexId pairable;
                        int num = numPairableFaces( Simplex(3,cofacesId[i]), facesName, pairable);

                        if(num == 1){
                            setPair( Simplex(2,pairable), Simplex(3,cofacesId[i]));
                            someonePaired = true;

                            cofacesName.erase(cofacesId[i]);
                            cofacesId[i] = -2;

                            facesId[facesName[pairable]] = -2;
                            facesName.erase(pairable);
                        }
                    }
                }

                if(!someonePaired){
                    for(int i=0; i<facesId.size(); i++){
                        if(facesId[i] >= 0){
                            criticalSimplices[2].insert(facesId[i]);
                            facesName.erase(facesId[i]);
                            facesId[i] = -2;
                        }
                    }
                }
            }

            for(int i=0; i<cofacesId.size(); i++){
                if(cofacesId[i] >= 0){
                    criticalSimplices[3].insert(cofacesId[i]);
                }
            }

//            cout << "Done second while" << endl;

        }
    }

    order.clear();

}


int FormanGradient::numPairableFaces(const Simplex& simplex, const map<SimplexId, int>& facesName, SimplexId& chosen){

    int tot=0;
    for(int i=0; i<simplex.first+1; i++){
        int face = extractBoundary(simplex,simplex.first-1,i);
        if(facesName.find(face) != facesName.end()){
            chosen = face;
            tot++;
        }
    }

  return tot;
}

int FormanGradient::extractBoundary(const Simplex& simpl,int dimension_b, int index_b){

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


int FormanGradient::sizeCofacets(const Simplex& simpl, int dimension_cb){

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
            if(dimensionality_ > 2)
                return triangulation_->getTriangleStarNumber(simpl.second);
            else
                return 0;
        case 3:
            return 0;
    }

}

int FormanGradient::extractCoboundary(const Simplex& simpl, int dimension_cb, int index_cb){

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

int FormanGradient::indexInsideCoface(const Simplex &simpl, const Simplex &coface) {

    for(int i=0; i<=coface.first; i++){
        if(extractBoundary(coface,simpl.first,i) == simpl.second)
            return i;
    }

    return -1;
}

void FormanGradient::getVertexIndexingOfVertices(Simplex simplex, vector<SimplexId>& vec){

    simplexToVertices(simplex,vec);

    for(SimplexId v=0; v<vec.size(); v++){
        vec[v] = (*indexing_)[vec[v]];
    }

}

bool FormanGradient::cmpSimplexesLevelSets(const Simplex& s1, const Simplex& s2){

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

void FormanGradient::setPair(const Simplex& tail, const Simplex& head){

    int index = indexInsideCoface(tail,head);
    (*(*gradient_)[tail.first])[(head.second*(head.first+1)) + index ] = true;
}

bool FormanGradient::getPairLower(const Simplex &simpl, Simplex &paired) {

    if(simpl.first != 0) {
        for (int i = 0; i <= simpl.first; i++) {
            if ((*(*gradient_)[simpl.first - 1])[(simpl.second * (simpl.first + 1)) + i]) {
                paired = Simplex(simpl.first - 1, extractBoundary(simpl, simpl.first - 1, i));
                return true;
            }
        }
    }

    return false;
}

bool FormanGradient::getPairHigher(const Simplex &simpl, Simplex &paired) {

    Simplex check;
    if(simpl.first<dimensionality_){
        for(int i=0; i<sizeCofacets(simpl,simpl.first+1); i++){
            Simplex coface = Simplex(simpl.first+1, extractCoboundary(simpl,simpl.first+1,i));
            if(getPairLower(coface,check)){
                if(check.second == simpl.second) {
                    paired = coface;
                    return true;
                }
            }
        }
    }

    return false;
}

bool FormanGradient::getPair(const Simplex& simpl, Simplex& paired){


    Simplex check;
    if(getPairLower(simpl, check)) {
        paired = check;
        return true;
    }

    if(getPairHigher(simpl, check)) {
        paired = check;
        return true;
    }

    return false;
}


bool FormanGradient::isCritical(const Simplex& simpl){

    Simplex check;
    return !getPairHigher(simpl,check) && !getPairLower(simpl,check);

}
