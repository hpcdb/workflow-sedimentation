import sys

def replace(inputFilePath,outputFilePath,forcing):
	readFile = open(inputFilePath,'r')
	writeFile = open(outputFilePath,'w+')

	for line in readFile:
		if "%=FORCING%" in line:
			line = line.replace("%=FORCING%",forcing)
		writeFile.write(line + "\n")
	writeFile.close()

if __name__ == "__main__":
	if(len(sys.argv) >= 3):
		# to debug
		inputFilePath = sys.argv[1]
		outputFilePath = sys.argv[2]
		forcing = sys.argv[3]
		
		replace(inputFilePath,outputFilePath,forcing)

