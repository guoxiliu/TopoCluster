#include <vector>
#include <utility>
#include <algorithm>
#include <map>

using namespace std;


class BoundaryMatrix{

    vector<vector<int> > matrix;

    vector<float> values;

    inline void addValue(int simplex, int boundary, float filtr){
        if(matrix.size() < simplex){
            matrix.push_back(vector<int>());
            values.push_back(filtr);
        }

        auto lb = std::lower_bound(matrix[simplex].begin(), matrix[simplex].end(), boundary);
        if (lb != matrix[simplex].end() && *lb == boundary) {
            auto ub = std::upper_bound(lb, matrix[simplex].end(), boundary);
            matrix[simplex].erase(lb, ub);
        }
        else{
            matrix[simplex].insert(lb,boundary);
        }
    }


    inline void sort(){

        //organize the simplices according to the filtration value
        map<float, vector<int> > organize_values;
        for(int i=0; i<values.size(); i++){
            auto it = organize_values.find(values[i]);
            if(it == organize_values.end())
                organize_values.insert(it, std::pair<float,vector<int> >(values[i],vector<int>()));

            it->second.push_back(i);
        }

        //compute the filtration value they should have in the matrix
        vector<int> filtration = vector<int>(values.size(), -1);
        int fI = 0;
        for(auto levelset : organize_values){
            vector<int> simplices = levelset.second;
            //here i should sort these guys

            for(int s=0; s < simplices.size(); s++){
                if(matrix[s].size() == 0){
                    filtration[s]=fI++;
                }
                else{
                    if( boundary_done(s, filtration)){
                        filtration[s]=fI++;
                    }
                    else{
                        int temp = simplices[s];
                        simplices[s] = simplices[s+1];
                        simplices[s+1] = temp;
                        s--;
                    }
                }
            }
        }
        organize_values.clear();

        //debug purposes
        for(auto s : filtration){
            assert(s != -1);
        }

        //sort the structure accordingly
        for(int i=0; i<matrix.size(); i++){
            if(i == filtration[i]){
                update_indices(matrix[i], filtration);
            }
            else{
                vector<int> temp1 = matrix[i];
                matrix[i] = matrix[filtration[i]];
                matrix[filtration[i]] = temp1;
                
                float val = values[i];
                values[i] = values[filtration[i]];
                values[filtration[i]] = val;
                i--;
            }
        }
    }

    inline bool boundary_done(int simpl, vector<int>& filtration){
        for(auto bd : matrix[simpl]){
            if(filtration[bd] == -1)
                return false;
        }
        return true;
    }

    inline void update_indices(vector<int>& vec, vector<int>& filtration){
        for(int i=0; i<vec.size(); i++)
            vec[i] = filtration[vec[i]];
        std::sort(vec.begin(), vec.end());
    }
};