#!/bin/bash
monetdb-stop-db
monetdb-destroy-db -f
monetdb-generate-db 
monetdb-run-sql-scripts database-script.sql
monetdb-start-db 