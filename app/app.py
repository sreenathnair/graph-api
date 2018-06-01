from flask import Flask, jsonify
from .model import graph
from .amino_acid_codes import amino_acid_codes
import re
from .utils import find_ranges
import more_itertools as mit
import json
from .sifts import get_uniprot, get_interpro, get_cath, get_scop, get_go, get_ec, get_pfam, get_uniprot_to_pfam, get_uniprot_segments, get_isoforms
from .residue import get_mappings_for_residue_uniprot, get_mappings_for_residue_cath, get_mappings_for_residue_interpro, get_mappings_for_residue_pfam, get_mappings_for_residue_scop
from .residue import get_mappings_for_residue_binding_site
from .compound import get_compound_atoms, get_compound_bonds, get_compound_in_pdb, get_compound_co_factors, get_compound_co_factors_het
from .validation import get_validation_protein_ramachandran_sidechain_outliers, get_validation_rama_sidechain_listing, get_validation_rna_pucker_suite_outliers
from .validation import get_validation_global_percentiles, get_validation_summary_quality_scores, get_validation_key_validation_stats, get_validation_xray_refine_data_stats
from .validation import get_validation_residuewise_outlier_summary
from werkzeug.contrib.cache import SimpleCache
from .pdb import get_binding_sites_for_entry, get_binding_sites_for_uniprot, get_secondary_structures
from .uniprot import get_unipdb, get_unipdb_residue
from flask import Response

app = Flask(__name__)

cache = SimpleCache()
CACHE_ENABLED = False
cache_timeout = 5 * 60

@app.route('/')
def default():

    return 'Default'


@app.route('/api/mappings/<string:entry_id>')
def get_mappings_api(entry_id):

    cache_result = cache.get('get_mappings_api')

    if(cache_result is None):

        cache_result = {
            entry_id: {
                "UniProt": get_uniprot(entry_id, graph),
                "Pfam": get_pfam(entry_id, graph),
                "CATH": get_cath(entry_id, graph),
                "InterPro": get_interpro(entry_id, graph),
                "HMMER": "",
                "SCOP": get_scop(entry_id, graph),
                "GO": get_go(entry_id, graph),
                "EC": get_ec(entry_id, graph)
            }
        }
        cache.set('get_mappings_api', cache_result, timeout=cache_timeout)

    return jsonify(cache_result)

@app.route('/api/mappings/uniprot/<string:entry_id>')
def get_uniprot_api(entry_id):

    cache_result = cache.get('get_uniprot_api:{}'.format(entry_id))

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

    return jsonify(cache_result), response_status


@app.route('/api/mappings/interpro/<string:entry_id>')
def get_interpro_api(entry_id):

    cache_result = cache.get('get_interpro_api:{}'.format(entry_id))

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "InterPro": get_interpro(entry_id, graph)
            }
        }
        cache.set('get_interpro_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)

    return jsonify(cache_result)


@app.route('/api/mappings/cath/<string:entry_id>')
def get_cath_api(entry_id):

    cache_result = cache.get('get_cath_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "CATH": get_cath(entry_id, graph)
            }
        }
        cache.set('get_cath_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/scop/<string:entry_id>')
def get_scop_api(entry_id):

    cache_result = cache.get('get_scop_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "SCOP": get_scop(entry_id, graph)
            }
        }
        cache.set('get_scop_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

    
@app.route('/api/mappings/go/<string:entry_id>')
def get_go_api(entry_id):

    cache_result = cache.get('get_go_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "GO": get_go(entry_id, graph)
            }
        }
        cache.set('get_go_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/ec/<string:entry_id>')
def get_ec_api(entry_id):

    cache_result = cache.get('get_ec_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "EC": get_ec(entry_id, graph)
            }
        }
        cache.set('get_ec_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/sequence_domains/<string:entry_id>')
def get_sequence_domains_api(entry_id):

    cache_result = cache.get('get_sequence_domains_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "InterPro": get_interpro(entry_id, graph),
                "Pfam": get_pfam(entry_id, graph)
            }
        }
        cache.set('get_sequence_domains_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

@app.route('/api/mappings/structural_domains/<string:entry_id>')
def get_structural_domains_api(entry_id):

    cache_result = cache.get('get_structural_domains_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "CATH": get_cath(entry_id, graph),
                "SCOP": get_scop(entry_id, graph)
            }
        }
        cache.set('get_structural_domains_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

@app.route('/api/mappings/uniprot_to_pfam/<string:accession>')
def get_uniprot_to_pfam_api(accession):

    cache_result = cache.get('get_uniprot_to_pfam_api:{}'.format(accession))

    if(cache_result is None):
        cache_result = {
            accession: {
                "Pfam": get_uniprot_to_pfam(accession, graph)
            }
        }
        cache.set('get_uniprot_to_pfam_api:{}'.format(accession), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

@app.route('/api/mappings/uniprot_segments/<string:entry_id>')
def get_uniprot_segments_api(entry_id):

    cache_result = cache.get('get_uniprot_segments_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: {
                "UniProt": get_uniprot_segments(entry_id, graph)
            }
        }
        cache.set('get_uniprot_segments_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/isoforms/<string:entry_id>')
def get_best_isoforms_api(entry_id):

    cache_result = cache.get('get_best_isoforms_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_isoforms(entry_id, graph, 'B')
        }
        cache.set('get_best_isoforms_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

@app.route('/api/mappings/all_isoforms/<string:entry_id>')
def get_all_isoforms_api(entry_id):

    cache_result = cache.get('get_all_isoforms_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_isoforms(entry_id, graph, 'A')
        }
        cache.set('get_all_isoforms_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/residue_mapping/<string:entry_id>/<string:entity_id>/<string:residue_number>')
def get_mappings_for_residue(entry_id, entity_id, residue_number):

    cache_result = cache.get('get_mappings_for_residue:{}:{}:{}'.format(entry_id, entity_id, residue_number))

    if(cache_result is None):
        cache_result = {
            entry_id: {}
        }

        cache_result[entry_id]["UniProt"] = get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number, graph)
        cache_result[entry_id]["Pfam"] = get_mappings_for_residue_pfam(entry_id, entity_id, residue_number, graph)
        cache_result[entry_id]["InterPro"] = get_mappings_for_residue_interpro(entry_id, entity_id, residue_number, graph)
        cache_result[entry_id]["CATH"] = get_mappings_for_residue_cath(entry_id, entity_id, residue_number, graph)
        cache_result[entry_id]["SCOP"] = get_mappings_for_residue_scop(entry_id, entity_id, residue_number, graph)
        cache_result[entry_id]["binding_sites"] = get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number, True, graph)

        cache.set('get_mappings_for_residue:{}:{}:{}'.format(entry_id, entity_id, residue_number), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

@app.route('/api/mappings/binding_sites/<string:entry_id>')
def get_binding_sites_for_entry_api(entry_id):

    cache_result = cache.get('get_binding_sites_for_entry_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_binding_sites_for_entry(entry_id, graph)
        }
        cache.set('get_binding_sites_for_entry_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)
    

@app.route('/api/mappings/uniprot/binding_sites/<string:uniprot_accession>')
def get_binding_sites_for_uniprot_api(uniprot_accession):

    cache_result = cache.get('get_binding_sites_for_uniprot_api:{}'.format(uniprot_accession))

    if(cache_result is None):
        cache_result = {
            uniprot_accession: get_binding_sites_for_uniprot(uniprot_accession, graph)
        }
        cache.set('get_binding_sites_for_uniprot_api:{}'.format(uniprot_accession), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/pdb/entry/secondary_structure/<string:entry_id>')
def get_secondary_structures_api(entry_id):

    cache_result = cache.get('get_secondary_structures_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_secondary_structures(entry_id, graph)
        }
        cache.set('get_secondary_structures_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/unipdb/<string:uniprot_accession>')
def get_unipdb_api(uniprot_accession):

    cache_result = cache.get('get_unipdb_api:{}'.format(uniprot_accession))

    if(cache_result is None):
        cache_result = {
            uniprot_accession: get_unipdb(uniprot_accession, graph)
        }
        cache.set('get_unipdb_api:{}'.format(uniprot_accession), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/mappings/unipdb/<string:uniprot_accession>/<string:unp_res>')
def get_unipdb_residue_api(uniprot_accession, unp_res):

    cache_result = cache.get('get_unipdb_residue_api:{}:{}'.format(uniprot_accession, unp_res))

    if(cache_result is None):
        cache_result = {
            uniprot_accession: get_unipdb_residue(uniprot_accession, unp_res, graph)
        }
        cache.set('get_unipdb_residue_api:{}:{}'.format(uniprot_accession, unp_res), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

    
@app.route('/api/pdb/compound/atoms/<string:chem_comp_id>')
def get_compound_atoms_api(chem_comp_id):

    cache_result = cache.get('get_compound_atoms_api:{}'.format(chem_comp_id))

    if(cache_result is None):
        cache_result = {
            chem_comp_id: get_compound_atoms(chem_comp_id, graph)
        }
        cache.set('get_compound_atoms_api:{}'.format(chem_comp_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/pdb/compound/bonds/<string:chem_comp_id>')
def get_compound_bonds_api(chem_comp_id):

    cache_result = cache.get('get_compound_bonds_api:{}'.format(chem_comp_id))

    if(cache_result is None):
        cache_result = {
            chem_comp_id: get_compound_bonds(chem_comp_id, graph)
        }
        cache.set('get_compound_bonds_api:{}'.format(chem_comp_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/pdb/compound/in_pdb/<string:chem_comp_id>')
def get_compound_in_pdb_api(chem_comp_id):

    cache_result = cache.get('get_compound_in_pdb_api:{}'.format(chem_comp_id))

    if(cache_result is None):
        cache_result = {
            chem_comp_id: get_compound_in_pdb(chem_comp_id, graph)
        }
        cache.set('get_compound_in_pdb_api:{}'.format(chem_comp_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)

@app.route('/api/validation/protein-ramachandran-sidechain-outliers/entry/<string:entry_id>')
def get_validation_protein_ramachandran_sidechain_outliers_api(entry_id):

    cache_result = cache.get('get_validation_protein_ramachandran_sidechain_outliers_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_validation_protein_ramachandran_sidechain_outliers(
                entry_id, graph)
        }
        cache.set('get_validation_protein_ramachandran_sidechain_outliers_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/validation/rama_sidechain_listing/entry/<string:entry_id>')
def get_validation_rama_sidechain_listing_api(entry_id):

    cache_result = cache.get('get_validation_rama_sidechain_listing_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_validation_rama_sidechain_listing(entry_id, graph)
        }
        cache.set('get_validation_rama_sidechain_listing_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/validation/RNA_pucker_suite_outliers/entry/<string:entry_id>')
def get_validation_rna_pucker_suite_outliers_api(entry_id):
    
    cache_result = cache.get('get_validation_rna_pucker_suite_outliers_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_validation_rna_pucker_suite_outliers(entry_id, graph)
        }
        cache.set('get_validation_rna_pucker_suite_outliers_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/pdb/compound/cofactors')
def get_compound_co_factors_api():

    cache_result = cache.get('get_compound_co_factors_api')

    if(cache_result is None):
        cache_result = get_compound_co_factors(graph)
        cache.set('get_compound_co_factors_api', cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/pdb/compound/cofactors/het/<string:het_code>')
def get_compound_co_factors_het_api(het_code):

    cache_result = cache.get('get_compound_co_factors_het_api:{}'.format(het_code))

    if(cache_result is None):
        cache_result = {
            het_code: [get_compound_co_factors_het(het_code, graph)]
        }
        cache.set('get_compound_co_factors_het_api:{}'.format(het_code), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/validation/global-percentiles/entry/<string:entry_id>')
def get_validation_global_percentiles_api(entry_id):

    cache_result = cache.get('get_validation_global_percentiles_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_validation_global_percentiles(entry_id, graph)
        }
        cache.set('get_validation_global_percentiles_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/validation/summary_quality_scores/entry/<string:entry_id>')
def get_validation_summary_quality_scores_api(entry_id):

    cache_result = cache.get('get_validation_summary_quality_scores_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_validation_summary_quality_scores(entry_id, graph)
        }
        cache.set('get_validation_summary_quality_scores_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/validation/key_validation_stats/entry/<string:entry_id>')
def get_validation_key_validation_stats_api(entry_id):

    cache_result = cache.get('get_validation_key_validation_stats_api:{}'.format(entry_id))

    if(cache_result is None):
        cache_result = {
            entry_id: get_validation_key_validation_stats(entry_id, graph)
        }
        cache.set('get_validation_key_validation_stats_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result)


@app.route('/api/validation/xray_refine_data_stats/entry/<string:entry_id>')
def get_validation_xray_refine_data_stats_api(entry_id):

    cache_result = cache.get('get_validation_xray_refine_data_stats_api:{}'.format(entry_id))

    if (not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_xray_refine_data_stats(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_validation_xray_refine_data_stats_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result), response_status
    

@app.route('/api/validation/residuewise_outlier_summary/entry/<string:entry_id>')
def get_validation_residuewise_outlier_summary_api(entry_id):

    cache_result = cache.get('get_validation_residuewise_outlier_summary_api:{}'.format(entry_id))

    if (not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        response, response_status = get_validation_residuewise_outlier_summary(entry_id, graph)
        cache_result = {
            entry_id: response
        }
        cache.set('get_validation_residuewise_outlier_summary_api:{}'.format(entry_id), cache_result, timeout=cache_timeout)
    
    return jsonify(cache_result), response_status