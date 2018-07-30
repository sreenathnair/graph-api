#!/usr/bin/env python

"""
app.py: Main entry point of API calls. Configures cache settings and invokes respective methods to retrieve data from Neo4J
"""

from flask import Flask, jsonify
from .model import graph
from .amino_acid_codes import amino_acid_codes
import re
from .utils import find_ranges
import more_itertools as mit
import json
from .sifts import get_uniprot, get_interpro, get_cath, get_scop, get_go, get_ec, get_pfam, get_uniprot_to_pfam, get_uniprot_segments, get_isoforms, get_best_structures, get_homologene
from  .sifts import get_ensembl, get_best_structures_residue_range, get_uniref90
from .residue import get_mappings_for_residue_uniprot, get_mappings_for_residue_cath, get_mappings_for_residue_interpro, get_mappings_for_residue_pfam, get_mappings_for_residue_scop
from .residue import get_mappings_for_residue_binding_site, get_basic_residue_details
from .compound import get_compound_atoms, get_compound_bonds, get_compound_in_pdb, get_compound_co_factors, get_compound_co_factors_het
from .validation import get_validation_protein_ramachandran_sidechain_outliers, get_validation_rama_sidechain_listing, get_validation_rna_pucker_suite_outliers
from .validation import get_validation_global_percentiles, get_validation_summary_quality_scores, get_validation_key_validation_stats, get_validation_xray_refine_data_stats
from .validation import get_validation_residuewise_outlier_summary, get_validation_protein_rna_dna_geometry_outlier_residues
from werkzeug.contrib.cache import SimpleCache
from .pdb import get_binding_sites_for_entry, get_binding_sites_for_uniprot, get_secondary_structures
from .uniprot import get_unipdb, get_unipdb_residue
from flask import Response

__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"


app = Flask(__name__)
app.config['JSONIFY_PRETTYPRINT_REGULAR'] = False

cache = SimpleCache()
CACHE_ENABLED = False
cache_timeout = 5 * 60

@app.route('/')
def default():

    return 'Default'


@app.route('/api/mappings/<string:entry_id>')
def get_mappings_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt, Pfam, InterPro, CATH, SCOP, IntEnz, GO, Ensembl and HMMER accessions (and vice versa).
    """

    cache_result = cache.get('get_mappings_api:{}'.format(entry_id))
    response_status = 200

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        
        unp_result, unp_status = get_uniprot(entry_id, graph)
        pfam_result, pfam_status = get_pfam(entry_id, graph)
        cath_result, cath_status = get_cath(entry_id, graph)
        interpro_result, interpro_status = get_interpro(entry_id, graph)
        scop_result, scop_status = get_scop(entry_id, graph)
        go_result, go_status = get_go(entry_id, graph)
        ec_result, ec_status = get_ec(entry_id, graph)

        all_status = [unp_status, pfam_status, cath_status, interpro_status, scop_status, go_status, ec_status]

        if(all_status == [404 for x in range(0, 7)]):
            return jsonify({}), 404

        cache_result = {
            entry_id: {
                "UniProt": unp_result,
                "Pfam": pfam_result,
                "CATH": cath_result,
                "InterPro": interpro_result,
                "HMMER": "",
                "SCOP": scop_result,
                "GO": go_result,
                "EC": ec_result
            }
        }
        cache.set('get_mappings_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/uniprot/<string:entry_id>')
def get_uniprot_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt.
    """

    cache_result = cache.get('get_uniprot_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):

        response, response_status = get_uniprot(entry_id, graph)
        cache_result = {
            entry_id: {
                "UniProt": response
            }
        }
        
        cache.set('get_uniprot_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/pfam/<string:entry_id>')
def get_pfam_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to Pfam.
    """

    cache_result = cache.get('get_pfam_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):

        response, response_status = get_pfam(entry_id, graph)
        cache_result = {
            entry_id: {
                "Pfam": response
            }
        }
        
        cache.set('get_pfam_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/interpro/<string:entry_id>')
def get_interpro_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to InterPro.
    """

    cache_result = cache.get('get_interpro_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_interpro(entry_id, graph)
        cache_result = {
            entry_id: {
                "InterPro": response
            }
        }
        cache.set('get_interpro_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/cath/<string:entry_id>')
def get_cath_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to CATH.
    """

    cache_result = cache.get('get_cath_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_cath(entry_id, graph)
        cache_result = {
            entry_id: {
                "CATH": response
            }
        }
        cache.set('get_cath_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/scop/<string:entry_id>')
def get_scop_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to SCOP.
    """

    cache_result = cache.get('get_scop_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_scop(entry_id, graph)
        cache_result = {
            entry_id: {
                "SCOP": response
            }
        }
        cache.set('get_scop_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

    
@app.route('/api/mappings/go/<string:entry_id>')
def get_go_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to GO.
    """

    cache_result = cache.get('get_go_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_go(entry_id, graph)
        cache_result = {
            entry_id: {
                "GO": response
            }
        }
        cache.set('get_go_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/ec/<string:entry_id>')
def get_ec_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to EC.
    """

    cache_result = cache.get('get_ec_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        ec_result, response_status = get_ec(entry_id, graph)
        cache_result = {
            entry_id: {
                "EC": ec_result
            }
        }
        cache.set('get_ec_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/sequence_domains/<string:entry_id>')
def get_sequence_domains_api(entry_id):
    """
    Get mappings from protein chains to both Pfam and InterPro as assigned by the SIFTS process.
    """

    cache_result = cache.get('get_sequence_domains_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        interpro_response, interpro_response_status = get_interpro(entry_id, graph)
        pfam_response, pfam_response_status = get_pfam(entry_id, graph)

        if(interpro_response_status != 200 and pfam_response_status != 200):
            response_status = interpro_response_status
        else:
            response_status = 200

        cache_result = {
            entry_id: {
                "InterPro": interpro_response,
                "Pfam": pfam_response
            }
        }
        cache.set('get_sequence_domains_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/structural_domains/<string:entry_id>')
def get_structural_domains_api(entry_id):
    """
    Get mappings from protein chains to both SCOP and CATH as assigned by the SIFTS process.
    """

    cache_result = cache.get('get_structural_domains_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        cath_response, cath_response_status = get_cath(entry_id, graph)
        scop_response, scop_response_status = get_scop(entry_id, graph)

        if(cath_response_status != 200 and scop_response_status != 200):
            response_status = cath_response_status
        else:
            response_status = 200

        cache_result = {
            entry_id: {
                "CATH": cath_response,
                "SCOP": scop_response
            }
        }
        cache.set('get_structural_domains_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/best_structures/<string:accession>')
def get_best_structures_api(accession):
    """
    Get the list of PDB structures mapping to a UniProt accession sorted by coverage of the protein and, if the same, resolution.
    """

    cache_result = cache.get('get_best_structures_api:{}'.format(accession))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_best_structures(accession, graph)
        cache_result = {
            accession: response
        }
        cache.set('get_best_structures_api:{}'.format(accession), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/best_structures/<string:accession>/<string:unp_start>/<string:unp_end>')
def get_best_structures_residue_range_api(accession, unp_start, unp_end):
    """
    Get the list of PDB structures mapping to a UniProt residue range sorted by coverage of the protein and, if the same, resolution.
    """

    cache_result = cache.get('get_best_structures_residue_range_api:{}:{}:{}'.format(accession, unp_start, unp_end))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_best_structures_residue_range(accession, unp_start, unp_end, graph)
        cache_result = {
            accession: response
        }
        cache.set('get_best_structures_residue_range_api:{}:{}:{}'.format(accession, unp_start, unp_end), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/homologene/<string:entry_id>/<string:entity_id>')
def get_homologene_api(entry_id, entity_id):
    """
    Get homologene polypeptides for a given PDB entity
    """

    cache_result = cache.get('get_homologene_api:{}:{}'.format(entry_id, entity_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_homologene(entry_id, entity_id, graph)
        cache_result = {
            entry_id +'_' +entity_id : response
        }
        cache.set('get_homologene_api:{}:{}'.format(entry_id, entity_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/ensembl/<string:entry_id>')
def get_ensembl_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to Ensembl
    """

    cache_result = cache.get('get_ensembl_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_ensembl(entry_id, graph)
        cache_result = {
            entry_id : response
        }
        cache.set('get_ensembl_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/uniprot_to_pfam/<string:accession>')
def get_uniprot_to_pfam_api(accession):
    """
    Get mappings (as assigned by the SIFTS process) from a UniProt accession to a Pfam accession with details of the Pfam protein family
    """

    cache_result = cache.get('get_uniprot_to_pfam_api:{}'.format(accession))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_uniprot_to_pfam(accession, graph)
        cache_result = {
            accession: {
                "Pfam": response
            }
        }
        cache.set('get_uniprot_to_pfam_api:{}'.format(accession), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/uniprot_segments/<string:entry_id>')
def get_uniprot_segments_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt. The segments are calculated using the same rules as with the standard call but 
    with the UniProt sequence as reference. The outcome is a set of segments that reflects discontinuities in the UniProt sequence.
    """

    cache_result = cache.get('get_uniprot_segments_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_uniprot_segments(entry_id, graph)
        cache_result = {
            entry_id: {
                "UniProt": response
            }
        }
        cache.set('get_uniprot_segments_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/isoforms/<string:entry_id>')
def get_best_isoforms_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt. It returns the best isoform found in UniProt.
    """

    cache_result = cache.get('get_best_isoforms_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_isoforms(entry_id, graph, 'B')
        cache_result = {
            entry_id: response
        }
        cache.set('get_best_isoforms_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/all_isoforms/<string:entry_id>')
def get_all_isoforms_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt. It returns all the mappings to the UniProt isoforms of the canonical accession.
    """

    cache_result = cache.get('get_all_isoforms_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_isoforms(entry_id, graph, 'A')
        cache_result = {
            entry_id: response
        }
        cache.set('get_all_isoforms_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/uniref90/<string:entry_id>')
def get_uniref90_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt. It returns mappings to all the UniRef90 members of the UniProt accession.
    """

    cache_result = cache.get('get_uniref90_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_uniref90(entry_id, graph, 'A')
        cache_result = {
            entry_id: response
        }
        cache.set('get_uniref90_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/homologene_uniref90/<string:entry_id>')
def get_homologene_uniref90_api(entry_id):
    """
    Get mappings (as assigned by the SIFTS process) from PDB structures to UniProt. It returns mappings to all the UniRef90 members of the UniProt accession which also belong 
    to the same Homologene clusters.
    """

    cache_result = cache.get('get_homologene_uniref90_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_uniref90(entry_id, graph, 'H')
        cache_result = {
            entry_id: response
        }
        cache.set('get_homologene_uniref90_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/residue_mapping/<string:entry_id>/<string:entity_id>/<string:residue_number>')
def get_mappings_for_residue(entry_id, entity_id, residue_number):
    """
    Get mappings (as assigned by the SIFTS process) for a PDB Residue to UniProt, Pfam, InterPro, CATH, SCOP, IntEnz, GO, Ensembl and HMMER accessions (and vice versa).
    """

    cache_result = cache.get('get_mappings_for_residue:{}:{}:{}'.format(entry_id, entity_id, residue_number))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "entity_id": int(entity_id),
                "chains": []
            }
        }

        unp_result, unp_status = get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number, graph)
        pfam_result, pfam_status = get_mappings_for_residue_pfam(entry_id, entity_id, residue_number, graph)
        interpro_result, interpro_status = get_mappings_for_residue_interpro(entry_id, entity_id, residue_number, graph)
        cath_result, cath_status = get_mappings_for_residue_cath(entry_id, entity_id, residue_number, graph)
        scop_result, scop_status = get_mappings_for_residue_scop(entry_id, entity_id, residue_number, graph)
        binding_sites_result, binding_sites_status = get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number, True, graph)

        if(unp_status != 200 and pfam_status != 200 and interpro_status != 200 and cath_status != 200 and scop_status != 200 and binding_sites_status != 200):
            response_status = unp_status
        else:
            response_status = 200

        chains, response_status = get_basic_residue_details(entry_id, entity_id, residue_number, graph)

        modified_chains = []
        for chain in chains:
            
            modified_residues = []

            for residue in chain["residues"]:
                
                residue["features"] = {}
                residue["features"]["UniProt"] = unp_result
                residue["features"]["Pfam"] = pfam_result
                residue["features"]["InterPro"] = interpro_result
                residue["features"]["CATH"] = cath_result
                residue["features"]["SCOP"] = scop_result
                residue["features"]["binding_sites"] = binding_sites_result

                modified_residues.append(residue)

            chain["residues"] = modified_residues
            modified_chains.append(chain)

        cache_result[entry_id]["chains"] = modified_chains

        cache.set('get_mappings_for_residue:{}:{}:{}'.format(entry_id, entity_id, residue_number), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/residue_mapping/<string:entry_id>/<string:entity_id>/<string:residue_start>/<string:residue_end>')
def get_mappings_for_residue_range(entry_id, entity_id, residue_start, residue_end):
    """
    Get mappings (as assigned by the SIFTS process) for a PDB Residue range to UniProt, Pfam, InterPro, CATH, SCOP, IntEnz, GO, Ensembl and HMMER accessions (and vice versa).
    """

    cache_result = cache.get('get_mappings_for_residue_range:{}:{}:{}:{}'.format(entry_id, entity_id, residue_start, residue_end))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "entity_id": int(entity_id),
                "chains": []
            }
        }

        dict_chains = {}

        for residue_number in range(int(residue_start), int(residue_end) + 1):
            
            residue_number = str(residue_number)

            unp_result, unp_status = get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number, graph)
            pfam_result, pfam_status = get_mappings_for_residue_pfam(entry_id, entity_id, residue_number, graph)
            interpro_result, interpro_status = get_mappings_for_residue_interpro(entry_id, entity_id, residue_number, graph)
            cath_result, cath_status = get_mappings_for_residue_cath(entry_id, entity_id, residue_number, graph)
            scop_result, scop_status = get_mappings_for_residue_scop(entry_id, entity_id, residue_number, graph)
            binding_sites_result, binding_sites_status = get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number, True, graph)

            if(unp_status != 200 and pfam_status != 200 and interpro_status != 200 and cath_status != 200 and scop_status != 200 and binding_sites_status != 200):
                response_status = unp_status
            else:
                response_status = 200

            chains, response_status = get_basic_residue_details(entry_id, entity_id, residue_number, graph)

            for chain in chains:
                
                auth_asym_id = chain["auth_asym_id"]
                struct_asym_id = chain["struct_asym_id"]
                chain_key = (auth_asym_id, struct_asym_id)

                if dict_chains.get(chain_key) is None:
                    dict_chains[chain_key] = {
                        "auth_asym_id": auth_asym_id,
                        "struct_asym_id": struct_asym_id,
                        "residues": []
                    }

                for residue in chain["residues"]:
                    
                    residue["features"] = {}
                    residue["features"]["UniProt"] = unp_result
                    residue["features"]["Pfam"] = pfam_result
                    residue["features"]["InterPro"] = interpro_result
                    residue["features"]["CATH"] = cath_result
                    residue["features"]["SCOP"] = scop_result
                    residue["features"]["binding_sites"] = binding_sites_result

                    dict_chains[chain_key]["residues"].append(residue)

            
        for key in dict_chains.keys():
            cache_result[entry_id]["chains"].append(dict_chains[key])

        cache.set('get_mappings_for_residue_range:{}:{}:{}:{}'.format(entry_id, entity_id, residue_start, residue_end), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/binding_sites/<string:entry_id>')
def get_binding_sites_for_entry_api(entry_id):
    """
    This call provides details on binding sites in the entry as per STRUCT_SITE records in PDB files (or mmcif equivalent thereof), such as ligand, residues in the site, 
    description of the site, etc.
    """

    cache_result = cache.get('get_binding_sites_for_entry_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_binding_sites_for_entry(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_binding_sites_for_entry_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status
    

@app.route('/api/mappings/uniprot/binding_sites/<string:uniprot_accession>')
def get_binding_sites_for_uniprot_api(uniprot_accession):
    """
    This call provides details on binding sites for a UniProt accession from related PDB entries as per STRUCT_SITE records in PDB files (or mmcif equivalent thereof), 
    such as ligand, residues in the site, description of the site, etc.
    """

    cache_result = cache.get('get_binding_sites_for_uniprot_api:{}'.format(uniprot_accession))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        unp_binding_result, response_status = get_binding_sites_for_uniprot(uniprot_accession, graph)
        cache_result = {
            uniprot_accession: unp_binding_result
        }
        cache.set('get_binding_sites_for_uniprot_api:{}'.format(uniprot_accession), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/pdb/entry/secondary_structure/<string:entry_id>')
def get_secondary_structures_api(entry_id):
    """
    This call provides details about residue ranges of regular secondary structure (alpha helices and beta strands) found in protein chains of the entry. For strands, 
    sheet id can be used to identify a beta sheet.
    """

    cache_result = cache.get('get_secondary_structures_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_secondary_structures(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_secondary_structures_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/unipdb/<string:uniprot_accession>')
def get_unipdb_api(uniprot_accession):

    cache_result = cache.get('get_unipdb_api:{}'.format(uniprot_accession))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_unipdb(uniprot_accession, graph)
        cache_result = {
            uniprot_accession: response
        }
        cache.set('get_unipdb_api:{}'.format(uniprot_accession), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/mappings/unipdb/<string:uniprot_accession>/<string:unp_res>')
def get_unipdb_residue_api(uniprot_accession, unp_res):
    """
    Get mappings (as assigned by the SIFTS process) for a UniProt Residue to sequence and structural domains.
    """

    cache_result = cache.get('get_unipdb_residue_api:{}:{}'.format(uniprot_accession, unp_res))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        unipdb_result, response_status = get_unipdb_residue(uniprot_accession, unp_res, graph)
        cache_result = {
            uniprot_accession: unipdb_result
        }
        cache.set('get_unipdb_residue_api:{}:{}'.format(uniprot_accession, unp_res), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

    
@app.route('/api/pdb/compound/atoms/<string:chem_comp_id>')
def get_compound_atoms_api(chem_comp_id):
    """
    This set of calls provides information about atoms in a chemical groups defined in the PDB Chemical Component Dictionary. For each atoms, properties such as name, element symbol, 
    ideal coordinates, stereochemistry, aromaticity (when applicable), etc. are available.
    """

    cache_result = cache.get('get_compound_atoms_api:{}'.format(chem_comp_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_compound_atoms(chem_comp_id, graph)
        cache_result = {
            chem_comp_id: response
        }
        cache.set('get_compound_atoms_api:{}'.format(chem_comp_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/pdb/compound/bonds/<string:chem_comp_id>')
def get_compound_bonds_api(chem_comp_id):
    """
    This set of calls provides information about bonds in a chemical groups defined in the PDB Chemical Component Dictionary. For each bond, properties such as atom names, bond type, 
    stereochemistry and aromaticity (when applicable) etc. are available.
    """

    cache_result = cache.get('get_compound_bonds_api:{}'.format(chem_comp_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_compound_bonds(chem_comp_id, graph)
        cache_result = {
            chem_comp_id: response
        }
        cache.set('get_compound_bonds_api:{}'.format(chem_comp_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/pdb/compound/in_pdb/<string:chem_comp_id>')
def get_compound_in_pdb_api(chem_comp_id):
    """
    This set of calls returns a list of PDB entries that contain the compound defined in the PDB Chemical Component Dictionary.
    """

    cache_result = cache.get('get_compound_in_pdb_api:{}'.format(chem_comp_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_compound_in_pdb(chem_comp_id, graph)
        cache_result = {
            chem_comp_id: response
        }
        cache.set('get_compound_in_pdb_api:{}'.format(chem_comp_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/pdb/compound/cofactors')
def get_compound_co_factors_api():
    """
    This call provides a summary of the cofactor annotation in the PDB.
    """

    cache_result = cache.get('get_compound_co_factors_api')
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        cache_result = get_compound_co_factors(graph)
        cache.set('get_compound_co_factors_api', cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/pdb/compound/cofactors/het/<string:het_code>')
def get_compound_co_factors_het_api(het_code):
    """
    This call provides hetcodes similar to the query parameter with the same functional annotation.
    """

    cache_result = cache.get('get_compound_co_factors_het_api:{}'.format(het_code))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        co_factor_result, response_status = get_compound_co_factors_het(het_code, graph)
        cache_result = {
            het_code: [] if co_factor_result is None else [co_factor_result]
        }
        cache.set('get_compound_co_factors_het_api:{}'.format(het_code), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/validation/protein-ramachandran-sidechain-outliers/entry/<string:entry_id>')
def get_validation_protein_ramachandran_sidechain_outliers_api(entry_id):
    """
    This call returns backbone and sidechain outliers in protien chains, as calculated by Molprobity as part of wwPDB validation pipeline.
    """

    cache_result = cache.get('get_validation_protein_ramachandran_sidechain_outliers_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_protein_ramachandran_sidechain_outliers(
                entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_validation_protein_ramachandran_sidechain_outliers_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/rama_sidechain_listing/entry/<string:entry_id>')
def get_validation_rama_sidechain_listing_api(entry_id):
    """
    This call returns Ramachandran status (favoured, outlier, etc.), phi-psi values, sidechain status (rotamer name or outlier) as reported by Molprobity component of the 
    wwPDB validation pipeline.
    """

    cache_result = cache.get('get_validation_rama_sidechain_listing_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_rama_sidechain_listing(entry_id, graph)

        if(response_status != 200):
            cache_result = {}
        else:
            cache_result = {
                entry_id: response
            }
        cache.set('get_validation_rama_sidechain_listing_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/RNA_pucker_suite_outliers/entry/<string:entry_id>')
def get_validation_rna_pucker_suite_outliers_api(entry_id):
    """
    This call returns RNA backbone outliers, i.e. non-rotameric suites and unusual puckers, as calculated by Molprobity as part of wwPDB validation pipeline.
    """
    
    cache_result = cache.get('get_validation_rna_pucker_suite_outliers_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_rna_pucker_suite_outliers(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_validation_rna_pucker_suite_outliers_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/global-percentiles/entry/<string:entry_id>')
def get_validation_global_percentiles_api(entry_id):
    """
    Metrics here are the ones recommended by validation task force. Global is against whole PDB archive and relative is against entries of comparable resolution.
    """

    cache_result = cache.get('get_validation_global_percentiles_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        val_result, response_status = get_validation_global_percentiles(entry_id, graph)
        cache_result = {
            entry_id: val_result
        }
        cache.set('get_validation_global_percentiles_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/summary_quality_scores/entry/<string:entry_id>')
def get_validation_summary_quality_scores_api(entry_id):
    """
    These scores are harmonic means of absolute percentiles of geometric metrics (e.g. ramachandran, clashscore, sidechains), reflections-based metrics (Rfree, RSRZ) and both these 
    kinds of metrics taken together. Wherever a constitutent percentile is 0, the harmonic mean is defined to be 0. When constituent percentiles are all unavailable, the harmonic mean 
    is null.
    """

    cache_result = cache.get('get_validation_summary_quality_scores_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        val_result, response_status = get_validation_summary_quality_scores(entry_id, graph)
        cache_result = {
            entry_id: val_result
        }
        cache.set('get_validation_summary_quality_scores_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/key_validation_stats/entry/<string:entry_id>')
def get_validation_key_validation_stats_api(entry_id):
    """
    This is still a very high level summary, but covers metrics of interest not included in percentiles, or a little more detail than just percentile.
    """

    cache_result = cache.get('get_validation_key_validation_stats_api:{}'.format(entry_id))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        val_result, response_status = get_validation_key_validation_stats(entry_id, graph)
        cache_result = {
            entry_id: val_result
        }
        cache.set('get_validation_key_validation_stats_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/xray_refine_data_stats/entry/<string:entry_id>')
def get_validation_xray_refine_data_stats_api(entry_id):

    cache_result = cache.get('get_validation_xray_refine_data_stats_api:{}'.format(entry_id))
    response_status = None

    if (not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_xray_refine_data_stats(entry_id, graph)

        if(response_status != 200):
            cache_result = {}
        else:
            cache_result = {
                entry_id: response
            }

        cache.set('get_validation_xray_refine_data_stats_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status
    

@app.route('/api/validation/residuewise_outlier_summary/entry/<string:entry_id>')
def get_validation_residuewise_outlier_summary_api(entry_id):
    """
    A residue can have many types of geometric or experimental-data-based outliers. This call lists all kinds of outliers found in a residue. For residues with no recorded outlier, 
    there is no information returned.
    """

    cache_result = cache.get('get_validation_residuewise_outlier_summary_api:{}'.format(entry_id))
    response_status = None

    if (not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_residuewise_outlier_summary(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_validation_residuewise_outlier_summary_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status


@app.route('/api/validation/protein-RNA-DNA-geometry-outlier-residues/entry/<string:entry_id>')
def get_validation_protein_rna_dna_geometry_outlier_residues_api(entry_id):
    """
    Lists residues in protein, DNA, RNA chains that contain various types of geometry outliers.
    """

    cache_result = cache.get('get_validation_protein_rna_dna_geometry_outlier_residues_api:{}'.format(entry_id))
    response_status = None

    if (not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_protein_rna_dna_geometry_outlier_residues(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_validation_protein_rna_dna_geometry_outlier_residues_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    else:
        response_status = 200
    
    return jsonify(cache_result), response_status
