import sys

singleSpace = " "
doubleSpace = "  "
floatAttributes = ["CSMAG","FLSIC","RTOL","DUTOL","DT","TMAX",
"CSI","SIGMA0","SIGMA1","GAMMA","ETAF","ETA",
"CFLMIN","CFLMAX","TOL_V","DTMIN","DTMAX","GRAVX","GRAVY","GRAVZ",
"BETA","TREF"]

startTerm = -4
# endTerm = 4
endTerm = -4

def extractData(filePath, inputFileName):
	global singleSpace

	inputFilePath = filePath + "/" + inputFileName
	readFile = open(inputFilePath,'r')
	
	data = {}
	attributes = []
	stop = False
	
	filePath = "model.part.in"
	writeFile = open(filePath,'w+')
	
	for line in readFile:
		if(len(attributes) == 5 and attributes[3] == 'FORCING'):
			line = line.strip()
			line = removeDoubleSpace(line)
			values = line.split(singleSpace)
			writeFile.write("  " + values[0] + "       " + values[1] + "             " + values[2] + "         %=FORCING%        " + values[4] + "\n")
		else:
			writeFile.write(line)

		if(line and not stop):
			if("$     TIME FUNCTIONS" in line):
				stop = True
			if(("$" in line) 
				and ("TITLE" not in line)):
				line = line.replace("$","").strip()
				line = removeDoubleSpace(line)
				attributes = line.split(singleSpace)
				if(len(attributes)==5 and attributes[3] == 'UNUSED'):
					attributes[3] = "FORCING"
			elif(len(attributes) > 0):
				line = line.strip()
				line = removeDoubleSpace(line)
				values = line.split(singleSpace)

				if(len(values) > 0 and len(attributes) >= len(values)):
					for i in range(len(values)):
						if(attributes[i] in floatAttributes):
							data[attributes[i]] = float(values[i])
						else:
							if(attributes[i] not in data.keys()):
								data[attributes[i]] = int(values[i])
				attributes = []
		else:
			attributes = []

	writeFile.close()

	data["INN"] = inputFilePath

	readFile.close()
	return data

def removeDoubleSpace(line):
	global singleSpace
	global doubleSpace
	if(doubleSpace in line):
		line = line.replace(doubleSpace,singleSpace)
		if(doubleSpace in line):
			line = removeDoubleSpace(line)
	return line

def generateOutputFile(path, fileName, attributes, data, dfa):
	filePath = path + "/" + fileName
	writeFile = open(filePath,'w+')

	header = ""
	values = ""
	
	first = True
	for att in attributes:
		if(first):
			first = False
		else:
			header += ";"
		header += att
	if 'FORCING' in attributes:
		header += ';FMODEL'
	writeFile.write(header + "\n")

	if 'FORCING' in attributes:
		for i in range(startTerm,endTerm+1):
			if i!=0:
				first = True
				values = ""
				for att in attributes:
					if(first):
						first = False
					else:
						values += ";"

					if(dfa):
						if(att == "INN"):
							values += "'" + str(data[att]) + "'"
						elif(att == 'FORCING'):
							values += str(i)
						else:	
							values += str(data[att])
					else:
						if(att == 'FORCING'):
							values += str(i)
						else:
							values += str(data[att])
				values += ";" + path + "/model.part.in"
				writeFile.write(values + "\n")
	else:
		first = True
		for att in attributes:
			if(first):
				first = False
			else:
				values += ";"

			if(dfa):
				if(att == "INN"):
					values += "'" + str(data[att]) + "'"
				else:	
					values += str(data[att])
			else:
				values += str(data[att])
		writeFile.write(values + "\n")
	writeFile.close()



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
		attributes = sys.argv[5].replace("[","").replace("]","").split(",")
		attributes += ["INN"]
		print "FILE_PATH = " + filePath
		print "INPUT_FILE_NAME = " + inputFileName
		print "OUTPUT_FILE_NAME = " + outputFileName

		extractedData = extractData(filePath,inputFileName)
		generateOutputFile(filePath,outputFileName,attributes,extractedData,dfa)


