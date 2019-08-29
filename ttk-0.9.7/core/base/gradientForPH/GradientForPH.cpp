#include                  <GradientForPH.h>

using namespace std;
using namespace ttk;

GradientForPH::GradientForPH(){

  triangulation_ = NULL;
  inputData_ = NULL;
  vertexIndexing = NULL;
  gradient = NULL;
  numberCriticalSimplices = 0;
}

GradientForPH::~GradientForPH(){

}

void GradientForPH::simplexToVertices(Simplex simplex, vector<SimplexId>& vertices){

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


SimplexId GradientForPH::getIndexing(Simplex simplex){

    vector<SimplexId> vertices;
    simplexToVertices(simplex,vertices);

    //from the vertices select the field value of the simplex
    SimplexId index = (*vertexIndexing)[vertices[0]];

    for(SimplexId v=1; v<vertices.size(); v++){
        index = index > (*vertexIndexing)[vertices[v]] ? index : (*vertexIndexing)[vertices[v]];
    }

    return index;
}

void GradientForPH::getVertexIndexingOfVertices(Simplex simplex, vector<SimplexId>& vec){

    simplexToVertices(simplex,vec);

    for(SimplexId v=0; v<vec.size(); v++){
        vec[v] = (*vertexIndexing)[vec[v]];
    }

}

void GradientForPH::computeLowerStar(SimplexId simpl, SimplexesVec& vec){

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

bool GradientForPH::cmpSimplexesLevelSets(const Simplex& s1, const Simplex& s2){

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


void GradientForPH::setPair(const Simplex& tail, const Simplex& head){


    if(tail.first == 0){
        (*gradient)[tail.first][tail.second][0] = head.second;
    }
    else{
       (*gradient)[tail.first][tail.second][1] = head.second;
    }

    (*gradient)[head.first][head.second][0] = tail.second;
}

bool GradientForPH::getPair(const Simplex& simpl, Simplex& paired){

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


bool GradientForPH::isCritical(const Simplex& simpl){

    for(auto i : (*gradient)[simpl.first][simpl.second]){
        if(i != -1)
            return false;
    }

    return true;
}

bool GradientForPH::isTopCritical(const Simplex& simplex){

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

void GradientForPH::computeBarycenter(const Simplex& simplex, vector<float>& coords){

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

float GradientForPH::computeDistance(const vector<float>& point1, const vector<float>& point2){

    float val=0;
    for(int i=0; i<point1.size(); i++)
        val += pow(point1[i]-point2[i],2);

    return sqrt(val);
}


int GradientForPH::extractBoundary(Simplex simpl,int dimension_b, int index_b){

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


int GradientForPH::sizeCofacets(Simplex simpl, int dimension_cb){

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

int GradientForPH::extractCoboundary(Simplex simpl, int dimension_cb, int index_cb){

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


