import os
import sys

filepath = sys.argv[1]
out = sys.argv[2]

with open(filepath, 'r') as input_file, open(out, 'w') as output_file:
    for line in input_file:
    	line = line.replace("\"", "").replace(":","").lower()
    	output_file.write(line)
input_file.close()
output_file.close()
