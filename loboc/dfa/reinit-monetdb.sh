#!/bin/bash
monetdb-stop-all
monetdb-start-all
monetdb-stop-db
monetdb-destroy-db -f
monetdb-generate-db 
monetdb-run-script-db database-script.sql
monetdb-start-db
monetdb-stop-all
