echo "Backing up provenance database..."
rm prov-db.dump
mclient -p 54321 -d dataflow_analyzer --dump > prov-db.dump
