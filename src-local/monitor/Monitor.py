import requests
import time
import os.path

from requests.exceptions import ConnectionError

class Monitor:
    """It defines an object for monitoring performance data from DfAnalyzer"""
    filePath = 'performance.log'
    tokenPath = 'finish_performance.tkn'

    def __init__(self, hostname="localhost", port=None):
        self.hostname = hostname
        self.port = port
        self.updateURL()
        # create file
        self.file = open(self.filePath, 'w+')

    def updateURL(self):
        """Update URL with the current hostname"""
        self.url = "http://" + self.hostname
        if self.port:
            self.url += ":" + str(self.port)
        return self.url

    def run(self):
        """It starts the daemon for monitoring the performance metrics of DfAnalyzer"""
        first = True
        while not os.path.exists(self.tokenPath):
            # send HTTP request
            response = self.sendHTTPRequest(first);
            # handle HTTP response
            self.handleResponse(response)
            if first and response:
                first = False
            # wait for a while
            time.sleep(1)
        self.file.close()

    def sendHTTPRequest(self, header=False):
        """Send HTTP request to DfAnalyzer"""
        # REST API - HTTP request with GET method
        request = self.url + "/performance/" + str(header).lower()
        try:
            response = requests.get(request)
            return response
        except ConnectionError as e:
            return None

    def handleResponse(self, response):
        """It handles HTTP responses from DfAnalyzer by writing a log with performance data."""
        if response:
            # print(response.text)
            self.file.write(response.text + "\n")
            self.file.flush()