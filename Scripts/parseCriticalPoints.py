#!/usr/bin/python3
from os import error
import sys
import csv

def parseCPResults(inputFile, outputFile):
  """
  Parse the result text file from running the TTK critical points plugin.
  The output file is in csv format.
  """
  allData = []

  # read the file
  with open(inputFile, "r") as fp:
    info = []
    row = []
    cache = 0
    preprocess = 0
    for line in fp:
      # here we get ttk result
      if line.startswith("Dataset"):
        dataset = line[line.find(' ')+1 : -1]
        if len(info) > 0:
          if len(row) > 0:
            row.insert(0, preprocess)
            info.extend(row)
          allData.append(info)
        # reset variables
        cache = 0
        preprocess = 0
        row = []
        info = [dataset]
      # here we have the cache size
      # elif cache == 0 and "Cache size" in line:
      #   cache = int(line[line.find(':')+1 : -1])
      #   info.append(cache)
      # here we get the time for preprocessing
      elif "Triangles processed in " in line:
        preprocess += float(line[line.find(' in ')+3 : line.find(' s')])
      elif "Time usage for preprocessing" in line:
        preprocess += (float(line[line.find(':')+1 : -3]))
      # here we get the memory usage
      elif line.startswith("[ScalarFieldCriticalPoints]") and "processed in " in line:
        row.append(float(line[line.find(' in ')+3 : line.find(' s')]))
      elif line.startswith("[ttkScalarFieldCriticalPoints] Memory usage"):
        row.append(float(line[line.find(':')+1 : -4]))

    # take care of the last one
    if len(info) > 0:
      if len(row) > 0:
        row.insert(0, preprocess)
        info.extend(row)
      allData.append(info)

    with open(outputFile, 'w', newline='') as fp:
      # header = ["Dataset", "Cache Size", "Preprocessing", "Total Time", "Memory Usage"]
      # for ttk
      header = ["Dataset", "Preprocessing", "Total Time", "Memory Usage"]
      csvWriter = csv.writer(fp, delimiter=',')
      csvWriter.writerow(header)
      csvWriter.writerows(allData)

def combineResults(prefix): 
  """
  Combine different csv files to a new csv file for comparison. 
  """
  pos = prefix.index('_')
  outputCSV = prefix[:pos] + "_cmp" + prefix[pos:] + ".csv"
  dataStructures = ["Explicit", "Implicit", "TTK"]
  with open(outputCSV, 'w') as wp: 
    csvWriter = csv.writer(wp, delimiter=',')
    # write header
    header = ["Dataset", "Preprocessing", "Total Time", "Memory Usage", "Data Structure"]
    csvWriter.writerow(header)
    # add the corresponding data structure
    for type in dataStructures:
      inputCSV = prefix+'_'+type.lower()+".csv"
      try: 
        with open(inputCSV, 'r') as rp:
          csvReader = csv.reader(rp)
          next(csvReader)   # skip the header
          for row in csvReader:
            row.append(type)
            csvWriter.writerow(row)
      except FileNotFoundError:
        print("Cannot find file: ", inputCSV)


if __name__ == "__main__":
  if len(sys.argv) < 2: 
    print("Usage: python3", sys.argv[0], "<type> [input_file (prefix)]  [output_file]")
    print("Type: ")
    print("1. Critical points results")
    print("2. Combine results from csv file")
    sys.exit()

  type = int(sys.argv[1])
  if type == 1: 
    if len(sys.argv) < 3: 
      print("Usage: python3", sys.argv[0], "1  <input_file>  [output_file]")
      sys.exit()
    inFile = sys.argv[2]
    if len(sys.argv) > 3:
      outFile = sys.argv[3]
    else:
      outFile = inFile.rsplit('.', 1)[0] + '.csv'
    parseCPResults(inFile, outFile)
  
  elif type == 2:
    filePrefix = "results_cp"
    if len(sys.argv) > 2:
      filePrefix = sys.argv[2]
    combineResults(filePrefix)
  
  else: 
    print("Wrong type argument, please double check it!")
    sys.exit()
  
  print("Finished!\n")
