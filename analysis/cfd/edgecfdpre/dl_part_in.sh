#!/bin/bash
PATH=`pwd`

time /usr/bin/python %=WFDIR%/../bin/Part_in-Extractor.py false $PATH %=MESH%_1.part.in dl_in.data [IPRINT,DT,ETAF,FORCING]