import os
from py2neo import Graph

from app.compound import get_compound_atoms, get_compound_bonds, get_compound_in_pdb, get_compound_co_factors, get_compound_co_factors_het
from app.model import get_neo4j_instance

print()
print('Running tests')

graph = get_neo4j_instance()
print('Using Neo4J instance -> {}'.format(graph))


def test_get_compound_atoms():
    response, response_status = get_compound_atoms('TDP', graph)
    assert not response_status == 500

def test_get_compound_bonds():
    response, response_status = get_compound_bonds('TDP', graph)
    assert not response_status == 500

def test_get_compound_in_pdb():
    response, response_status = get_compound_in_pdb('TDP', graph)
    assert not response_status == 500

def test_get_compound_co_factors():
    response, response_status = get_compound_co_factors(graph)
    assert not response_status == 500

def test_get_compound_co_factors_het():
    response, response_status = get_compound_co_factors_het('TDP', graph)
    assert not response_status == 500