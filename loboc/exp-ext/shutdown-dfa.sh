#!/bin/bash
db_host=$1
curl -X POST http://$db_host:22000/shutdown