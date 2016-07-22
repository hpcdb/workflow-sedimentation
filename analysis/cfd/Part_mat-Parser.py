import sys

def replace(inputFilePath,outputFilePath,visc,dens,kxx,kyy,kzz):
	readFile = open(inputFilePath,'r')
	writeFile = open(outputFilePath,'w+')

	for line in readFile:
		if "%=VISC%" in line:
			line = line.replace("%=VISC%",visc)
		if "%=DENS%" in line:
			line = line.replace("%=DENS%",dens)
		if "%=KXX%" in line:
			line = line.replace("%=KXX%",kxx)
		if "%=KYY%" in line:
			line = line.replace("%=KYY%",kyy)
		if "%=KZZ%" in line:
			line = line.replace("%=KZZ%",kzz)
		
		writeFile.write(line + "\n")
	writeFile.close()

if __name__ == "__main__":
	if(len(sys.argv) >= 3):
		# to debug
		inputFilePath = sys.argv[1]
		outputFilePath = sys.argv[2]
		visc = sys.argv[3]
		dens = sys.argv[4]
		kxx = sys.argv[5]
		kyy = sys.argv[6]
		kzz = sys.argv[7]
		
		replace(inputFilePath,outputFilePath,visc,dens,kxx,kyy,kzz)

