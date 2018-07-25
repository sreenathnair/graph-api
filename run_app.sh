
if [ $# -ne 1 ]; then
    echo "Usage: ./run_app.sh <source neo4j (prod/local)>"
    exit
fi

source_db=$1

if [ $source_db == "local" ]; then
    export DEST_URL=bolt://localhost:7687
    export NEO4J_USERNAME=neo4j
    export NEO4J_PASSWORD=root
elif [ $source_db == "prod" ]; then
    export DEST_URL=bolt://wp-np2-3d.ebi.ac.uk:7687
    export NEO4J_USERNAME=neo4j
    export NEO4J_PASSWORD=pdbe_neo
else
    echo "Illegal source Neo4J instance"
    exit
fi

source venv/bin/activate
python run.py