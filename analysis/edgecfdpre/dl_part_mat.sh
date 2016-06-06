#!/bin/bash
PATH=`pwd`

time /usr/bin/python %=WFDIR%/../bin/Part_mat-Extractor.py false $PATH %=MESH%_1.part.mat dl_mat.data [VISC,DENS,KXX,KYY,KZZ]

