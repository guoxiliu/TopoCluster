#include                  <FG_Segmentation.h>
#include                <FG_Segmentation_template.h>

using namespace std;
using namespace ttk;

FG_Segmentation::FG_Segmentation(){

}

FG_Segmentation::~FG_Segmentation(){
  
}

void FG_Segmentation::readCriticalPoints(vector<float>& points, 
                   vector<char>& dimension, 
                   vector<vector<int> >& indexes){

    for(int i=0;i<=dimensionality_; i++)
        indexes.push_back(vector<int>());

    for(int i=0; i<=dimensionality_; i++){

        int numSimplices;
        switch (i)
        {
        case 0:
            numSimplices = triangulation_->getNumberOfVertices();
            break;
        case 1:
            numSimplices = triangulation_->getNumberOfEdges();
            break;
        case 2:
            numSimplices = triangulation_->getNumberOfTriangles();
            break;
        case 3:
            numSimplices = triangulation_->getNumberOfCells();
            break;
        
        default:
            break;
        }
    
        for(int j=0; j<numSimplices; j++){
            Simplex simpl = Simplex(i,j);
            if(isCritical(simpl)){
                vector<float> barycenter;
                computeBarycenter(simpl,barycenter);
                points.push_back(barycenter[0]);
                points.push_back(barycenter[1]);
                points.push_back(barycenter[2]);
                dimension.push_back(i);
                indexes[i].push_back(j);  
            }
        }
    }
}

void FG_Segmentation::extractDescendingCell(Simplex simpl, list<Simplex>& descendingCell, set<SimplexId>& vertices){
    
    vector<SimplexId> verts;
    int dim = simpl.first;
    queue<SimplexId> visitCell;
    visitCell.push(simpl.second);

    unordered_set<SimplexId> visited;
    visited.insert(simpl.second);

    descendingCell.push_back(simpl);
    simplexToVertices(simpl,verts);
    vertices.insert(verts.begin(), verts.end());

    while(!visitCell.empty()){

        SimplexId simpl = visitCell.front();
        visitCell.pop();

        for(int i=0; i<dim+1; i++){
            SimplexId face = extractBoundary(Simplex(dim,simpl),dim-1,i);
            
            Simplex pair;
            bool isPaired = getPair(Simplex(dim-1,face), pair);

            if(isPaired && 
               pair.first == dim && 
               pair.second != simpl && 
               visited.find(pair.second) == visited.end()){

                visitCell.push(pair.second);
                descendingCell.push_back(pair);
                simplexToVertices(pair,verts);
                vertices.insert(verts.begin(), verts.end());
                visited.insert(pair.second);
            }
        }
    }
}


void FG_Segmentation::extractAscendingCell(Simplex simpl, list<Simplex>& ascendingCell, set<SimplexId>& vertices){

    vector<SimplexId> verts;
    int dim = simpl.first;

    queue<SimplexId> visitCell;
    visitCell.push(simpl.second);

    unordered_set<SimplexId> visited;
    visited.insert(simpl.second);

    while(!visitCell.empty()){

        SimplexId simplex = visitCell.front();
        visitCell.pop();

        for(int i=0; i<sizeCofacets(Simplex(dim,simplex), dim+1); i++){
            SimplexId coface = extractCoboundary(Simplex(dim,simplex),dim+1,i);
            
            Simplex pair;
            Simplex head = Simplex(dim+1,coface);
            bool isPaired = getPair(head, pair);

            if(isPaired && pair.first == dim && pair.second != simplex && visited.find(pair.second) == visited.end()){
                visitCell.push(pair.second);
                visited.insert(head.second);

                ascendingCell.push_back(head);
                simplexToVertices(head,verts);
                vertices.insert(verts.begin(), verts.end());
            }
            else if(!isPaired){
                simplexToVertices(head,verts);
                vertices.insert(verts.begin(), verts.end());
                ascendingCell.push_back(head);
            }
        }
    }
}

void FG_Segmentation::computeMorseSmale() {

    boost::dynamic_bitset<> ascending1Cells(triangulation_->getNumberOfTriangles());
    boost::dynamic_bitset<> ascending2Cells(triangulation_->getNumberOfEdges());
    boost::dynamic_bitset<> ascending3Cells(triangulation_->getNumberOfVertices());
    boost::dynamic_bitset<> descending1Cells(triangulation_->getNumberOfEdges());
    boost::dynamic_bitset<> descending2Cells(triangulation_->getNumberOfTriangles());
    boost::dynamic_bitset<> descending3Cells(triangulation_->getNumberOfCells());

    for(int i=0; i<criticalSimplices.size(); i++){
        cout << "Dimension " << i << endl;
        for(auto simpl_id : criticalSimplices[i]){

            Simplex simpl(i,simpl_id);
            switch(i){
                case 0:
                    markAscendingCell(simpl,ascending3Cells);
                    break;
                case 1:
                    markDescendingCell(simpl,descending1Cells);
                    markAscendingCell(simpl,ascending2Cells);
                    break;
                case 2:
                    markDescendingCell(simpl,descending2Cells);
                    markAscendingCell(simpl,ascending1Cells);
                    break;
                case 3:
                    markDescendingCell(simpl,descending3Cells);
                    break;
                default:
                    cout << "Wrong dimension " << endl;
                    break;
            }
        }
    }

    cout << "For ascending cells " << endl;
    cout << "   vertices touched " << ascending3Cells.count() << " " << ascending3Cells.size() << endl;
    cout << "   edges touched " << ascending2Cells.count() << " " << ascending2Cells.size() << endl;
    cout << "   triangles touched " << ascending1Cells.count() << " " << ascending1Cells.size() << endl;

    cout << "For descending cells " << endl;
    cout << "   tetra touched " << descending3Cells.count() << " " << descending3Cells.size() << endl;
    cout << "   triangles touched " << descending2Cells.count() << " " << descending2Cells.size() << endl;
    cout << "   edges touched " << descending1Cells.count() << " " << descending1Cells.size() << endl;

}

void FG_Segmentation::markDescendingCell(Simplex simpl, boost::dynamic_bitset<>& descendingCells){

    int dim = simpl.first;
    queue<SimplexId> visitCell;
    visitCell.push(simpl.second);

    descendingCells[simpl.second]=true;

    while(!visitCell.empty()) {

        SimplexId simpl = visitCell.front();
        visitCell.pop();
        descendingCells[simpl] = true;

        for (int i = 0; i < dim + 1; i++) {
            SimplexId face = extractBoundary(Simplex(dim, simpl), dim - 1, i);

            Simplex pair;
            bool isPaired = getPair(Simplex(dim - 1, face), pair);

            if (isPaired && pair.first == dim && pair.second != simpl) {
                visitCell.push(pair.second);
            } else if (!isPaired) {
//                descendingCells[pair.second] = true;
            }
        }
    }
}

void FG_Segmentation::markAscendingCell(Simplex simpl, boost::dynamic_bitset<>& ascendingCells){

    int dim = simpl.first;
    queue<SimplexId> visitCell;
    visitCell.push(simpl.second);

    ascendingCells[simpl.second]=true;

    while(!visitCell.empty()) {

        SimplexId simpl_id = visitCell.front();
        visitCell.pop();
        ascendingCells[simpl_id] = true;

        Simplex simpl2(dim,simpl_id);
        for (int i = 0; i < sizeCofacets(simpl2, dim+1); i++) {
            SimplexId coface = extractCoboundary(simpl2, dim + 1, i);

            Simplex pair;
            bool isPaired = getPair(Simplex(dim + 1, coface), pair);

            if (isPaired && pair.first == dim && pair.second != simpl_id) {
                visitCell.push(pair.second);
            } else if (!isPaired) {
//                ascendingCells[pair.second] = true;
            }
        }
    }
}