from py2neo import Graph, Node, Relationship
import os

#url = os.environ.get('GRAPHENEDB_URL', 'http://miranda.ebi.ac.uk:7474')
#url = os.environ.get('GRAPHENEDB_URL', 'http://localhost:7474')
#dest = 'bolt://miranda.ebi.ac.uk:7687'
#dest = 'bolt://localhost:7687'
dest = os.environ.get('DEST_URL')
username = os.environ.get('NEO4J_USERNAME')
password = os.environ.get('NEO4J_PASSWORD')

graph = Graph(dest, username=username, password=password)

