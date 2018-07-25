#!/usr/bin/env python

"""
amino_acid_codes.py: This have dictionary of amino acid three letter codes and respective one letter code and name
"""

__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"

amino_acid_codes = {

    'ALA': ('A', 'Alanine'),
    'ARG':	('R', 'Arginine'),
    'ASN':	('N', 'Asparagine'),
    'ASP':	('D', 'Aspartic acid'),
    'CYS': ('C', 'Cysteine'),
    'GLN': ('Q', 'Glutamine'),
    'GLU': ('E', 'Glutamic acid'),
    'GLY': ('G', 'Glycine'),
    'HIS': ('H', 'Histidine'),
    'ILE': ('I', 'Isoleucine'),
    'LEU': ('L', 'Leucine'),
    'LYS': ('K', 'Lysine'),
    'MET': ('M', 'Methionine'),
    'PHE': ('F', 'Phenylalanine'),
    'PRO': ('P', 'Proline'),
    'PYL': ('O', 'Pyrrolysine'),
    'SER': ('S', 'Serine'),
    'SEC': ('U', 'Selenocysteine'),
    'THR': ('T', 'Threonine'),
    'TRP': ('W', 'Tryptophan'),
    'TYR': ('Y', 'Tyrosine'),
    'VAL': ('V', 'Valine'),
    'ASX': ('B', 'Aspartic acid or Asparagine'),
    'GLX': ('Z', 'Glutamic acid or Glutamine'),
    'XAA': ('X', 'Any amino acid'),
    'XLE': ('J', 'Leucine or Isoleucine'),
    'TERM':	('', 'termination codon'),
    'CCS': ('CCS', 'CCS'),
    'CME': ('CME', 'CME'),
    'MLZ': ('MLZ', 'MLZ'),
    'CSO': ('CSO', 'CSO'),
    'TPO': ('TPO', 'TPO'),
    'SCH': ('SCH', 'SCH')

}
