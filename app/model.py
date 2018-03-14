from py2neo import Graph, Node, Relationship
import os

#url = os.environ.get('GRAPHENEDB_URL', 'http://miranda.ebi.ac.uk:7474')
url = os.environ.get('GRAPHENEDB_URL', 'http://localhost:7474')
username = os.environ.get('NEO4J_USERNAME')
password = os.environ.get('NEO4J_PASSWORD')

graph = Graph(url + '/db/data/', username=username, password=password)

