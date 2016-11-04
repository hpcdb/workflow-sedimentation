#!/bin/bash
monetdb-stop-db
monetdb-destroy-db -f
monetdb-generate-db 
monetdb-run-script-db database-script.sql
monetdb-start-db