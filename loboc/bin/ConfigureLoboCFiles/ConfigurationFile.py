class ConfigurationFile():
    """A configuration file with properties of computational nodes in a cluster environment. 
    More practical, this code considers the configuration file of LoboC."""

    filePath = None
    databaseNode = None
    computingNodes = []

    # constants
    databaseFileName = "database.conf"
    computingFileName = "nodes.conf"

    def __init__(self, filePath):
        """Define a new configuration file"""
        self.filePath  = filePath;

    def organizeNodes(self, nodes):
        """Organize computational nodes to the database and computing processes."""
        for node in nodes:
            if(self.databaseNode == None):
                self.databaseNode = node
	    else:
            	self.computingNodes += [node]

    def read(self):
        """Read computational nodes from configuration file"""
        nodes = []
        with open(self.filePath) as f:
            for line in f:
                if line not in nodes:
                    nodes += [line]
        self.organizeNodes(nodes)

    def writeConfigurationFiles(self):
        """Write configuration files based on database and computing nodes."""
        # database
        file = open(self.databaseFileName,"w")
        file.write(self.databaseNode)
        file.close()
        # computing nodes
        file = open(self.computingFileName, "w")
        file.write("".join(str(x) for x in self.computingNodes))
        file.close()

