import csv
import numpy as np
import json


def read_csv_file(filename):
  with open(filename, 'r') as f:
    reader = csv.reader(f)
    array = []
    for row in reader:
      #ignore row if header or footer
      if len(row[0]) > 0 and row[0] == "#":
        continue
      array.append(list(map(float,row)))
    return np.array(array)

def write_csv_file(ndarray, filename, header=""):
  with open(filename, 'w') as f:
    if len(header) > 0:
      f.write("#" + header + "\n")
    writer = csv.writer(f)
    if len(ndarray.shape) == 1:
      writer.writerow(ndarray)
    else:
      writer.writerows(ndarray)

def dict_to_json(data, filename):
  with open(filename, 'w') as f:
    f.write(json.dumps(data))