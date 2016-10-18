echo "Backing up provenance database..."
rm prov-db.dump
mclient -p 50000 -d dataflow_analyzer --dump > prov-db.dump
