import os
from py2neo import Graph

from app.compound import get_compound_atoms
from app.model import get_neo4j_instance

print()
print('Running tests')

graph = get_neo4j_instance()
print('Using Neo4J instance -> {}'.format(graph))


def test_get_compound_atoms():
    response, response_status = get_compound_atoms('TDP', graph)
    assert not response_status == 500

    