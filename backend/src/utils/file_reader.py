import csv
import numpy as np

def read_csv_file(filename):
  with open(filename, 'r') as f:
    reader = csv.reader(f)
    array = []
    for row in reader:
      array.append(list(map(float,row)))
    return np.array(array)