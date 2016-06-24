import os
import sys

filepath = sys.argv[0]
out = sys.argv[1]

with open(filepath, 'r') as input_file, open(out, 'w') as output_file:
    for line in input_file:
    	line = line.replace("\"", "").replace(":","").lower()
    	output_file.write(line)
os.remove(filepath)