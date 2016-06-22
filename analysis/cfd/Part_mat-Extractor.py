import sys

singleSpace = " "
doubleSpace = "  "
intAttributes = ["ID"]
strAttributes = ["NAME"]
separator = ";"

startTerm = 1
# endTerm = 10
endTerm = 1

def extractData(path, inputFileName, outputFileName, attributes, dfa):
	global singleSpace

	inputFilePath = path + "/" + inputFileName
	outputFilePath = path + "/" + outputFileName
		
	readFile = open(inputFilePath,'r')
	writeFile = open(outputFilePath,'w+')
	
	indexes = []
	for line in readFile:
		if(line):
			if(("$" in line) 
				and ("$=" not in line)):
				line = line.replace("$","").strip()
				line = removeDoubleSpace(line)
				atts = ["MAT"]
				atts += line.split(singleSpace)

				first = True
				for att in attributes:
					indexes += [atts.index(att)]
					if(first):
						first = False
					else:
						writeFile.write(";")
					writeFile.write(att)
				writeFile.write(";MMODEL")
				writeFile.write("\n")
			elif(len(indexes) > 0 
				and ("$=" not in line)):
				line = line.strip()
				line = removeDoubleSpace(line)
				if(line):
					viscIndex = atts.index("VISC")
					for visc in range(startTerm,endTerm+1):
						values = [inputFilePath]
						values += line.split(singleSpace)
						first = True
						for i in indexes:
							prefix = ""
							suffix = ""
							if(first):
								first = False
								if(dfa):
									prefix = "'"
									suffix = "'"
							else:
								writeFile.write(";")

							if i==viscIndex:
								writeFile.write(prefix + str(visc*0.0001) + suffix)
							else:
								writeFile.write(prefix + values[i] + suffix)
						writeFile.write(";" + path + "/model.part.mat" + "\n")

	readFile.close()
	writeFile.close()

def removeDoubleSpace(line):
	global singleSpace
	global doubleSpace
	if(doubleSpace in line):
		line = line.replace(doubleSpace,singleSpace)
		if(doubleSpace in line):
			line = removeDoubleSpace(line)
	return line



if __name__ == "__main__":
	if(len(sys.argv) > 3):
		# to debug
		# filePath = "/Users/vitor/Documents/Repository/Thesis/EdgeCFD-trunk/analysis/paraview/hdf5"
		# fileName = "cav.2_1.part.in"
		dfaArg = sys.argv[1]
		dfa = False
		if(dfaArg == "true"):
			dfa = True

		filePath = sys.argv[2]
		inputFileName = sys.argv[3]
		outputFileName = sys.argv[4]
		attributes = ["MAT"]
		attributes += sys.argv[5].replace("[","").replace("]","").split(",")
		print "FILE_PATH = " + filePath
		print "INPUT_FILE_NAME = " + inputFileName
		print "OUTPUT_FILE_NAME = " + outputFileName
		print "ATTRIBUTES = ", attributes

		extractedData = extractData(filePath,inputFileName,outputFileName,attributes,dfa)


