from flask import Flask, jsonify
from .model import graph
from .amino_acid_codes import amino_acid_codes
import re
from .utils import find_ranges
import more_itertools as mit
import json
from .sifts import get_uniprot, get_interpro, get_cath, get_scop, get_go, get_ec, get_pfam, get_uniprot_to_pfam, get_uniprot_segments, get_isoforms, get_best_structures
from .residue import get_mappings_for_residue_uniprot, get_mappings_for_residue_cath, get_mappings_for_residue_interpro, get_mappings_for_residue_pfam, get_mappings_for_residue_scop
from .residue import get_mappings_for_residue_binding_site
from .compound import get_compound_atoms, get_compound_bonds, get_compound_in_pdb, get_compound_co_factors, get_compound_co_factors_het
from .validation import get_validation_protein_ramachandran_sidechain_outliers, get_validation_rama_sidechain_listing, get_validation_rna_pucker_suite_outliers
from .validation import get_validation_global_percentiles, get_validation_summary_quality_scores, get_validation_key_validation_stats, get_validation_xray_refine_data_stats
from .validation import get_validation_residuewise_outlier_summary, get_validation_protein_rna_dna_geometry_outlier_residues
from werkzeug.contrib.cache import SimpleCache
from .pdb import get_binding_sites_for_entry, get_binding_sites_for_uniprot, get_secondary_structures
from .uniprot import get_unipdb, get_unipdb_residue
from flask import Response

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

@app.route('/api/mappings/uniprot_to_pfam/<string:accession>')
def get_uniprot_to_pfam_api(accession):

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


@app.route('/api/mappings/residue_mapping/<string:entry_id>/<string:entity_id>/<string:residue_number>')
def get_mappings_for_residue(entry_id, entity_id, residue_number):

    cache_result = cache.get('get_mappings_for_residue:{}:{}:{}'.format(entry_id, entity_id, residue_number))
    response_status = None

    if(not CACHE_ENABLED):
        cache_result = None

    if(cache_result is None):
        cache_result = {
            entry_id: {}
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

        cache_result[entry_id]["UniProt"] = unp_result
        cache_result[entry_id]["Pfam"] = pfam_result
        cache_result[entry_id]["InterPro"] = interpro_result
        cache_result[entry_id]["CATH"] = cath_result
        cache_result[entry_id]["SCOP"] = scop_result
        cache_result[entry_id]["binding_sites"] = binding_sites_result

        cache.set('get_mappings_for_residue:{}:{}:{}'.format(entry_id, entity_id, residue_number), cache_result, timeout=cache_timeout)
    else:
        response_status = 200

    return jsonify(cache_result), response_status

@app.route('/api/mappings/binding_sites/<string:entry_id>')
def get_binding_sites_for_entry_api(entry_id):

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

@app.route('/api/validation/protein-ramachandran-sidechain-outliers/entry/<string:entry_id>')
def get_validation_protein_ramachandran_sidechain_outliers_api(entry_id):

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


@app.route('/api/pdb/compound/cofactors')
def get_compound_co_factors_api():

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


@app.route('/api/validation/global-percentiles/entry/<string:entry_id>')
def get_validation_global_percentiles_api(entry_id):

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
