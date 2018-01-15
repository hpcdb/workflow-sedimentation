import sys

from Monitor import Monitor

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def printOptionalArguments():
    """Print optional arguments to the Python application"""
    print("python Main.py [-h <dfa_hostname> -p <dfa_port>]")

if __name__ == "__main__":
    # variables
    hostname = "localhost"
    port = 22000
    # handle input arguments
    printOptionalArguments()
    myargs = getopts(sys.argv)
    if '-h' in myargs:
        hostname = myargs['-h']
    elif '-p' in myargs:
        port = int(myargs['-p'])
    # monitor
    print("Starting Performance Monitor...")
    monitor = Monitor(hostname,port)
    monitor.run()
    print("Shutting down Performance Monitor...")