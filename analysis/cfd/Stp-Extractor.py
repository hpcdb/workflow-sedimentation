import sys
import os

singleSpace = " "
doubleSpace = "  "

def removeDoubleSpace(line):
	global singleSpace
	global doubleSpace
	if(doubleSpace in line):
		line = line.replace(doubleSpace,singleSpace)
		if(doubleSpace in line):
			line = removeDoubleSpace(line)
	return line

def readDataFromStpFile(inputFilePath,outputFilePath):
	readFile = open(inputFilePath,'r')
	writeFile = open(outputFilePath,'a')

	if os.path.isfile(outputFilePath):
  		writeFile = open(outputFilePath,'w+')
	
	writeFile.write("TIME;NON_LINEAR_ITERS;BACKTRACKING_STEPS;REJECTED_BACKTRACKING_STEPS;LINEAR_ITERS;AVG_TOLER;R;RRO\n")

	for line in readFile:
		line = line.strip()
		line = removeDoubleSpace(line)
		attributes = line.split(singleSpace)
		attributes[0] = str(round(float(attributes[0]),1))
		if(attributes[5] == "*********"):
			attributes[5] = ""
		else:
			attributes[5] = "%.7f" % float(attributes[5])
		attributes[6] = "%.7f" % float(attributes[6])
		attributes[7] = "%.7f" % float(attributes[7])
		writeFile.write(attributes[0] + ";"
			+ attributes[1] + ";"
			+ attributes[2] + ";"
			+ attributes[3] + ";"
			+ attributes[4] + ";"
			+ attributes[5] + ";"
			+ attributes[6] + ";"
			+ attributes[7] + "\n")

	readFile.close()
	writeFile.close()

if __name__ == "__main__":
	if(len(sys.argv) > 1):
		# to debug
		inputFilePath = sys.argv[1]
		outputFilePath = sys.argv[2]

		print "###################"
		print "Stp Extractor"
		print "INPUT_FILE_PATH = " + inputFilePath
		print "OUTPUT_FILE_PATH = " + outputFilePath

		readDataFromStpFile(inputFilePath,outputFilePath)


