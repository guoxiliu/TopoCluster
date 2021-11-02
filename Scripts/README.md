## Testing Scripts

### In This Folder
The folder contains the convenient scripts to run the experiments of scalar field critical points plugin with TopoCluster and TTK Triangulation data structures. 

- `runCriticalPoints.sh` is a bash script to run the critical points plugin on different data sets and saves all the output to a text file. The user has to specify which data strcuture to test when running the script. 
- `parseCriticalPoints.py` is a Python script to parse the text file saved by `runCriticalPoints.sh`, which extracts the useful information from the text file and saves it to a csv file. It also provides a function to combine the results of all three data structures: Explicit TopoCluster, Implicit TopoCluster, and TTK Triangulation. 
- `drawFigure.py` is a Python script to generate the comparison figure from the combined csv file. 

### How to Use
Before using the provided scripts to test TopoCluster, please make sure you follow the [installation guide](../README.md) to install ParaView and TTK. 

#### Download the experimental data sets
Please check out the information of the data sets [here](../Datasets/README.md), and download them from the provided link. 

### Add execution rights to scripts (optional)
Use command `chmod +x <script_name>` to add the execution permission to the script. After this, you can directly use `./ <script_name>` in the terminal to run the script. Otherwise you might need to use command `bash <bash_script>` or `python3 <python3_script>` to run the script.

#### Get the text files of results

- Check the variable `DATASET_FOLDER` in file `runCriticalPoints.sh` to make sure the path of the downloaded data sets is correct. 
- Set `ENABLE_IMPLICIT_TOPOCLUSTER` to `OFF` to install Explicit TopoCluster and run `./runCritialPoints.sh 1` to get the results for Explicit TopoCluster. The results are saved as `results_cp_explicit.txt`.
- Set `ENABLE_IMPLICIT_TOPOCLUSTER` to `ON` to install Implicit TopoCluster and run `./runCritialPoints.sh 2` to get the results for Implicit TopoCluster. The results are saved as `results_cp_implicit.txt`. 
- run `./runCriticalPoints.sh 3` to get the results for TTK Triangulation. The results are saved as `results_cp_ttk.txt`. 

#### Generate the csv files

- Run command `./parseCriticalPoints.py 1 <text_file>` to get the corresponding csv file for each text file obtained in the last step. For example, `./parseCriticalPoints.py results_cp_explicit.txt` will generate the csv file for Explicit TopoCluster. 
- After all three csv files are generated, you can use command `./parseCriticalPoints.py 2` to combine all three csv files to a single csv file named `results_cmp_cp.csv`.

#### Draw figures! 
Just run command `./drawFigure.py` and the scripts will generate the figures similar to Fig. 7 in the paper. 

