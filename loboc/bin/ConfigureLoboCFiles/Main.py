import sys

from ConfigurationFile import ConfigurationFile

def main(argv):
    if(len(argv) > 0):
        configFile = ConfigurationFile(argv[0])
        configFile.read()
        configFile.writeConfigurationFiles()
    else:
        printUsage()

def printUsage():
    """"Print correct usage of arguments in this project."""
    print("python ConfigureLoboCFiles/Main.py nodes.txt")

if __name__ == "__main__": main(sys.argv[1:])