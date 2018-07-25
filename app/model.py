#!/usr/bin/env python

"""
model.py: Retrieves an instance of Neo4J
"""

from py2neo import Graph, Node, Relationship
import os

__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"

dest = os.environ.get('DEST_URL')
username = os.environ.get('NEO4J_USERNAME')
password = os.environ.get('NEO4J_PASSWORD')

graph = Graph(dest, username=username, password=password)


def get_neo4j_instance():
    return graph