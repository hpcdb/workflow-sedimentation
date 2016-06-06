from time import sleep
import sys
import os
import shutil

def joinData(inner,outer,outputFilePath):
	ireader = open(inner,'r')
	oreader = open(outer,'r')
	print "ok"

	iatts = []
	oatts = []
	headerLine = ""
	for iline in ireader:
		for oline in oreader:
			iatts = iline.split(";")
			oatts = oline.split(";")
			headerLine = iline[len('TIME')+1:-1] + ";" + oline
			break
		break

	iindex = None
	oindex = None
	for i in range(len(iatts)):
		for o in range(len(oatts)):
			if(iatts[i] == oatts[o] and iatts[i] == "TIME"):
				iindex = i
				oindex = o
				break
		if(iindex != None):
			break

	print iindex
	print oindex

	writeFile = None
  	if not os.path.exists(outputFilePath):
  		writeFile = open(outputFilePath,'w+')
  		writeFile.write(headerLine)
	else:
		writeFile = open(outputFilePath,'a')

	ireader = open(inner,'r')
	ifirst = True
	for iline in ireader:
		if(ifirst):
			ifirst = False
		else:	
			ivals = iline.split(";")
			oreader = open(outer,'r')
			ofirst = True
			for oline in oreader:
				if(ofirst):
					ofirst = False
				else:
					ovals = oline.split(";")
					# print str(ivals[iindex]) + " - " + str(ovals[oindex])
					if(ivals[iindex] == ovals[oindex]):
						wline = iline[len(ivals[iindex])+1:-1] + ";" + oline
						writeFile.write(wline)

	ireader.close()
	oreader.close()
	writeFile.close()

if __name__ == "__main__":
	if(len(sys.argv) >= 3):
		# to debug
		inner = sys.argv[1]
		outer = sys.argv[2]
		outputFilePath = sys.argv[3]

		print "###################"
		print "Join Data"
		print "STP_FILE_PATH = " + inner
		print "PARAVIEW_FILE_PATH = " + outer
		print "OUTPUT_FILE_PATH = " + outputFilePath

		joinData(inner,outer,outputFilePath)



