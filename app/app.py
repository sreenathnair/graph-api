from flask import Flask, jsonify
from .model import graph
from .amino_acid_codes import amino_acid_codes
import re
from .utils import find_ranges
import more_itertools as mit
import json
from .sifts import get_uniprot, get_interpro, get_cath, get_scop, get_go, get_ec, get_pfam, get_uniprot_to_pfam, get_uniprot_segments
from .residue import get_mappings_for_residue_uniprot, get_mappings_for_residue_cath, get_mappings_for_residue_interpro, get_mappings_for_residue_pfam, get_mappings_for_residue_scop
from .residue import get_mappings_for_residue_binding_site
from .compound import get_compound_atoms, get_compound_bonds, get_compound_in_pdb
from .validation import get_validation_protein_ramachandran_sidechain_outliers

app = Flask(__name__)

@app.route('/')
def default():

    return 'Default'

@app.route('/api/mappings/uniprot/<string:entry_id>')
def get_uniprot_api(entry_id):

    return jsonify({
        entry_id: {
            "UniProt": get_uniprot(entry_id, graph)
        }
    })


@app.route('/api/mappings/interpro/<string:entry_id>')
def get_interpro_api(entry_id):

    return jsonify({
        entry_id: {
            "InterPro": get_interpro(entry_id, graph)
        }
    })

@app.route('/api/mappings/cath/<string:entry_id>')
def get_cath_api(entry_id):

    return jsonify({
        entry_id: {
            "CATH": get_cath(entry_id, graph)
        }
    })
    

@app.route('/api/mappings/scop/<string:entry_id>')
def get_scop_api(entry_id):

    return jsonify({
        entry_id: {
            "SCOP": get_scop(entry_id, graph)
        }
    })


@app.route('/api/mappings/go/<string:entry_id>')
def get_go_api(entry_id):

    return jsonify({
        entry_id: {
            "GO": get_go(entry_id, graph)
        }
    })
    

@app.route('/api/mappings/ec/<string:entry_id>')
def get_ec_api(entry_id):

    return jsonify({
        entry_id: {
            "EC": get_ec(entry_id, graph)
        }
    })

@app.route('/api/mappings/sequence_domains/<string:entry_id>')
def get_sequence_domains_api(entry_id):

    return jsonify({
        entry_id: {
            "InterPro": get_interpro(entry_id, graph),
            "Pfam": get_pfam(entry_id, graph)
        }
    })
    

@app.route('/api/mappings/structural_domains/<string:entry_id>')
def get_structural_domains_api(entry_id):

    return jsonify({
        entry_id: {
            "CATH": get_cath(entry_id, graph),
            "SCOP": get_scop(entry_id, graph)
        }
    })

@app.route('/api/mappings/uniprot_to_pfam/<string:accession>')
def get_uniprot_to_pfam_api(accession):

    return jsonify({
        accession: {
            "Pfam": get_uniprot_to_pfam(accession, graph)
        }
    })

@app.route('/api/mappings/uniprot_segments/<string:entry_id>')
def get_uniprot_segments_api(entry_id):

    return jsonify({
        entry_id: {
            "UniProt": get_uniprot_segments(entry_id, graph)
        }
    })

@app.route('/api/mappings/residue_mapping/<string:entry_id>/<string:entity_id>/<string:residue_number>')
def get_mappings_for_residue(entry_id, entity_id, residue_number):

    final_result = { 
        entry_id: {} 
    }

    final_result[entry_id]["UniProt"] = get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number, graph)
    final_result[entry_id]["Pfam"] = get_mappings_for_residue_pfam(entry_id, entity_id, residue_number, graph)
    final_result[entry_id]["InterPro"] = get_mappings_for_residue_interpro(entry_id, entity_id, residue_number, graph)
    final_result[entry_id]["CATH"] = get_mappings_for_residue_cath(entry_id, entity_id, residue_number, graph)
    final_result[entry_id]["SCOP"] = get_mappings_for_residue_scop(entry_id, entity_id, residue_number, graph)
    final_result[entry_id]["binding_sites"] = get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number, True, graph)

    return jsonify(final_result)


@app.route('/api/mappings/binding_sites/<string:entry_id>')
def get_binding_sites_for_entry(entry_id):

    site_dict = {}
    site_ligand_dict = {}
    site_boundligand_dict = {}
    site_residue_dict = {}
    final_result = []

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[relation:HAS_BINDING_SITE]->(site:Binding_Site)
    OPTIONAL MATCH (site)<-[ligand_site_relation:IS_IN_BINDING_SITE]-(ligand_chain:Chain)<-[ligand_relation:CONTAINS_CHAIN]-(ligand_entity:Entity)-[:IS_AN_INSTANCE_OF]->(ligand:Ligand)
    OPTIONAL MATCH (site)-[:BOUNDED_BY]->(boundligand_chain:Chain)<-[bound_relation:CONTAINS_CHAIN]-(boundligand_entity:Entity)-[:IS_AN_INSTANCE_OF]->(boundligand:Ligand)
    RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, boundligand.ID, boundligand.NAME, boundligand.FORMULA, boundligand_chain.AUTH_ASYM_ID, boundligand_chain.STRUCT_ASYM_ID, 
    bound_relation.AUTH_SEQ_ID, boundligand_entity.ID, bound_relation.RES_ID, ligand.ID, ligand.NAME, ligand.FORMULA, ligand_chain.AUTH_ASYM_ID, ligand_relation.AUTH_SEQ_ID, 
    ligand_relation.STRUCT_ASYM_ID, ligand_site_relation.SYMMETRY_SYMBOL, ligand_entity.ID, ligand_relation.RES_ID ORDER BY site.ID
    """

    mappings = list(graph.run(query, parameters= {
        'entry_id': str(entry_id)
    }))

    for mapping in mappings:
        
        (site_id, site_name, site_evidence, bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, bound_auth_seq_id, bound_entity_id,
        bound_residue_id, ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, 
        ligand_residue_id) = mapping

        boundligand_to_list = (bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, 
                                                              bound_auth_seq_id, bound_entity_id, bound_residue_id)

        ligand_to_list = (ligand, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id)

        if(site_dict.get(site_id) is None):
            site_dict[site_id] = (site_name, site_evidence)

        if(site_boundligand_dict.get(site_id) is None):
            site_boundligand_dict[site_id] = [boundligand_to_list]
        else:
            site_boundligand_dict[site_id].append(boundligand_to_list)

        # usage of optional match may cause ligands with null values, so ignore those
        if(ligand is not None):
            if(site_ligand_dict.get(site_id) is None):
                site_ligand_dict[site_id] = [ligand_to_list]
            else:
                site_ligand_dict[site_id].append(ligand_to_list)


    del mappings

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[relation:HAS_BINDING_SITE]->(site:Binding_Site)
    MATCH (entry)-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[res_relation:IS_IN_BINDING_SITE]->(site)
    WITH entity, site, pdb_res, res_relation
    RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, pdb_res.ID, entity.ID, pdb_res.CHEM_COMP_ID, res_relation.AUTH_ASYM_ID, 
    res_relation.AUTH_SEQ_ID, res_relation.STRUCT_ASYM_ID, res_relation.SYMMETRY_SYMBOL ORDER BY site.ID
    """

    mappings = list(graph.run(query, parameters= {
        'entry_id': str(entry_id)
    }))

    
    for mapping in mappings:

    
        (site_id, site_name, site_evidence, pdb_res, pdb_res_entity_id, pdb_res_chem_comp_id, pdb_res_auth_asym_id, pdb_res_auth_seq_id, pdb_res_struct_asym_id, 
        pdb_res_sym_symbol) = mapping

        residue_to_list = (pdb_res_chem_comp_id, pdb_res_auth_asym_id, pdb_res_auth_seq_id, pdb_res_struct_asym_id, pdb_res_sym_symbol, pdb_res_entity_id, pdb_res)

        if(site_ligand_dict.get(site_id) is None):
            site_ligand_dict[site_id] = [residue_to_list]
        else:
            site_ligand_dict[site_id].append(residue_to_list)

    for key in site_dict.keys():
        (site_name, evidence) = site_dict[key]
        temp = {
            "site_id": key,
            "evidence_code": evidence,
            "details": site_name,
            "site_residues": [],
            "ligand_residues": []
        }

        for result in list(set(site_boundligand_dict[key])):
            (bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, bound_auth_seq_id, bound_entity_id, bound_residue_id) = result
            
            temp["ligand_residues"].append({
                "entity_id": int(bound_entity_id),
                "residue_number": int(bound_residue_id),
                "author_insertion_code": "null",
                "chain_id": bound_auth_asym_id,
                "author_residue_number": int(bound_auth_seq_id),
                "chem_comp_id": bound_ligand,
                "struct_asym_id": bound_struct_asym_id
            })
        
        for result in list(set(site_ligand_dict[key])):
            (ligand, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id) = result
            
            temp["site_residues"].append({
                "entity_id": int(ligand_entity_id),
                "residue_number": int(ligand_residue_id),
                "author_insertion_code": "null",
                "chain_id": ligand_auth_asym_id,
                "author_residue_number": int(ligand_auth_seq_id),
                "chem_comp_id": ligand,
                "struct_asym_id": ligand_struct_asym_id,
                "symmetry_symbol": ligand_sym_symbol
            })
        
        final_result.append(temp)

    return jsonify({
        entry_id: final_result
    })

@app.route('/api/mappings/uniprot/binding_sites/<string:uniprot_accession>')
def get_binding_sites_for_uniprot(uniprot_accession):

    final_result = []

    query = """
    MATCH (unp:UniProt {ACCESSION: $accession})-[:HAS_UNP_RESIDUE]-(unp_res:UNP_Residue)<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)
    MATCH (pdb_res)-[:IS_IN_BINDING_SITE]->(site:Binding_Site)
    WITH pdb_res, unp_res
    MATCH (pdb_res)<-[:HAS_PDB_RESIDUE]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry)
    RETURN DISTINCT unp_res.ID, unp_res.ONE_LETTER_CODE, pdb_res.ID, pdb_res.CHEM_COMP_ID, entity.ID, entry.ID
    """

    mappings = list(graph.run(query, parameters= {
        'accession': str(uniprot_accession)
    }))

    for mapping in mappings:
        (unp_res_id, unp_res_code, pdb_res_id, pdb_res_code, entity_id, entry_id) = mapping

        residue_result = get_mappings_for_residue_binding_site(entry_id, entity_id, pdb_res_id, False, graph)

        temp_result = {
            "unp_res_id": unp_res_id,
            "unp_res_code": unp_res_code,
            "pdb_res_id": pdb_res_id,
            "pdb_res_code": pdb_res_code,
            "entry_id": entry_id,
            "entity_id": entity_id,
            "binding_sites": residue_result
        }

        final_result.append(temp_result)

    return jsonify({
        uniprot_accession: final_result
    })


@app.route('/api/mappings/secondary_structures/<string:pdb_id>')
def get_secondary_structures(pdb_id):

    query = """
    MATCH (n:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_CHAIN]->(chain:Chain)
    WHERE rel.OBSERVED='Y' AND (rel.IS_IN_HELIX='Y' OR exists(rel.SHEETS))
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, rel.HELIX_SEGMENT, rel.SHEETS, toInteger(pdb_res.ID), toInteger(rel.AUTH_SEQ_ID) ORDER by toInteger(pdb_res.ID)
    """

    mappings = list(graph.run(query, parameters= {
        'entry_id': str(pdb_id)
    }))

    dict_helices = {}
    dict_strands_res_id = {}
    dict_strands_res_auth_seq_id = {}
    dict_struct_asym_id = {}

    for mapping in mappings:
        (entity_id, auth_asym_id, struct_asym_id, helix_segment, sheets, pdb_residue, auth_seq_id) = mapping

        #setting dictionary of helices
        if(helix_segment != '0'):
            if(dict_helices.get((entity_id, struct_asym_id, helix_segment)) is None):
                dict_helices[(entity_id, struct_asym_id, helix_segment)] = [(pdb_residue, auth_seq_id)]
            else:
                dict_helices[(entity_id, struct_asym_id, helix_segment)].append((pdb_residue, auth_seq_id))

        #setting dictionary of strands
        if(sheets is not None):
            sheets = re.sub(r'[\[\]]','', sheets).split(',')
            
            for sheet in sheets:
                if(dict_strands_res_id.get((entity_id, struct_asym_id, sheet)) is None):
                    dict_strands_res_id[(entity_id, struct_asym_id, sheet)] = [pdb_residue]
                    dict_strands_res_auth_seq_id[(entity_id, struct_asym_id, sheet)] = [auth_seq_id]
                else:
                    dict_strands_res_id[(entity_id, struct_asym_id, sheet)].append(pdb_residue)
                    dict_strands_res_auth_seq_id[(entity_id, struct_asym_id, sheet)].append(auth_seq_id)

        dict_struct_asym_id[struct_asym_id] = auth_asym_id

    final_result = {
        pdb_id: {
            "molecules": []
        }
    }

    dict_helices_chain = {}
    dict_entity = {}
    
    for key in dict_helices.keys():
        (entity_id, struct_asym_id, helix_segment) = key

        pdb_start_id, pdb_start_auth = dict_helices.get(key)[0]
        pdb_end_id, pdb_end_auth = dict_helices.get(key)[-1]
        temp_segment = {
                "start": {
                    "author_residue_number": int(pdb_start_auth),
                    "author_insertion_code": None,
                    "residue_number": int(pdb_start_id)
                },
                "end": {
                    "author_residue_number": int(pdb_end_auth),
                    "author_insertion_code": None,
                    "residue_number": int(pdb_end_id)
                }
            }
        if(dict_helices_chain.get((entity_id, struct_asym_id)) is None):
            dict_helices_chain[(entity_id, struct_asym_id)] = [temp_segment]
        else:
            dict_helices_chain[(entity_id, struct_asym_id)].append(temp_segment)

        if(dict_entity.get(entity_id) is None):
            dict_entity[entity_id] = set(list(struct_asym_id))
        else:
            dict_entity[entity_id].add(struct_asym_id)

    dict_strand_chain = {}

    for key in dict_strands_res_id.keys():
        (entity_id, struct_asym_id, sheet) = key
        res_id_list = []
        res_auth_seq_list = []

        for group in mit.consecutive_groups(dict_strands_res_id[key]):
            element = list(group)
            pdb_start = element[0]
            pdb_end = element[-1]
            res_id_list.append((pdb_start, pdb_end))

        for group in mit.consecutive_groups(dict_strands_res_auth_seq_id[key]):
            element = list(group)
            pdb_start = element[0]
            pdb_end = element[-1]
            res_auth_seq_list.append((pdb_start, pdb_end))
        
        for incr in range(len(res_id_list)):
            temp_strand = {
                "start": {
                    "author_residue_number": int(res_auth_seq_list[incr][0]),
                    "author_insertion_code": None,
                    "residue_number": int(res_id_list[incr][0])
                },
                "end": {
                    "author_residue_number": int(res_auth_seq_list[incr][1]),
                    "author_insertion_code": None,
                    "residue_number": int(res_id_list[incr][1])
                },
                "sheet_id": sheet
            }
            incr += 1

            if(dict_strand_chain.get((entity_id, struct_asym_id)) is None):
                dict_strand_chain[(entity_id, struct_asym_id)] = [temp_strand]
            else:
                dict_strand_chain[(entity_id, struct_asym_id)].append(temp_strand)


    for entity in dict_entity.keys():
        temp_entity = {
            "entity_id": entity,
            "chains": []
        }
        for chain in dict_entity[entity]:
            temp_chain = {
                "chain_id": dict_struct_asym_id[chain],
                "struct_asym_id": chain,
                "secondary_structure": {
                    "helices": [],
                    "strands": []
                }
            }
            for helix in dict_helices_chain[(entity, chain)]:
                temp_chain["secondary_structure"]["helices"].append(helix)
            for strand in dict_strand_chain[(entity, chain)]:
                temp_chain["secondary_structure"]["strands"].append(strand)
            
            temp_entity["chains"].append(temp_chain)
        
        final_result[pdb_id]["molecules"].append(temp_entity)
    
    return jsonify(final_result)

   

@app.route('/api/mappings/unipdb/<string:uniprot_accession>')
def get_unipdb(uniprot_accession):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)<-[relation:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)
    <-[:HAS_PDB_RESIDUE]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry)
    WITH entry.ID AS entryId, entity.ID AS entityId, relation.CHAINS AS chains, toInteger(pdb_res.ID) AS pdbRes, toInteger(unp_res.ID) AS unpRes, pdb_res.CHEM_COMP_ID AS pdbResCode, unp_res.ONE_LETTER_CODE AS unpResCode 
    RETURN entryId, entityId, chains, pdbRes, unpRes, pdbResCode, unpResCode ORDER BY entryId, pdbRes
    """

    mappings = list(graph.run(query, parameters= {
        'accession': str(uniprot_accession)
    }))

    final_result = {
        uniprot_accession: {
            'Pfam': get_mappings_for_unp_residue_pfam(uniprot_accession),
            'mappings': []
        }
    }

    entry_entity_dict = {}


    for mapping in mappings:
        
        (entry, entity, chains, pdb_res, unp_res, pdb_res_code, unp_res_code) = mapping
        pdb_res_one_letter, pdb_res_desc = amino_acid_codes[pdb_res_code]
        del pdb_res_desc
        data = (chains, pdb_res, unp_res, pdb_res_one_letter, unp_res_code)

        # create array of mappings in dict for an (entry, entity) pair if not present
        if entry_entity_dict.get((entry, entity)) is None:
            entry_entity_dict[(entry, entity)] = [data]
        else:
            entry_entity_dict[(entry, entity)].append(data)


    final_map = {}

    for key in entry_entity_dict.keys():
        final_map[key] = []
        pdb_seq = ''
        unp_seq = ''
        prev_pdb_res = None
        prev_unp_res = None
        pdb_start = None
        unp_start = None
        incr = 0
        start = True

        while incr < len(entry_entity_dict[key]):
            (chains, pdb_res, unp_res, pdb_res_code, unp_res_code) = entry_entity_dict[key][incr]
           
            if start or prev_pdb_res == pdb_res - 1:
                if start:
                    if prev_pdb_res != None:
                        pdb_start = prev_pdb_res
                        unp_start = prev_unp_res
                    else:
                        pdb_start = pdb_res
                        unp_start = unp_res
                pdb_seq += pdb_res_code
                unp_seq += unp_res_code
                start = False
            else:
                # check for a new segment
                final_map[key].append((chains, pdb_seq, unp_seq, pdb_start, prev_pdb_res, unp_start, prev_unp_res)) # pdb_res and unp_end will be ending residue numbers
                pdb_seq = pdb_res_code
                unp_seq = unp_res_code
                start = True
            
            # the very last residue
            if incr == len(entry_entity_dict[key]) - 1:
                final_map[key].append((chains, pdb_seq, unp_seq, pdb_start, pdb_res, unp_start, unp_res)) # pdb_res and unp_end will be ending residue numbers

            prev_pdb_res = pdb_res
            prev_unp_res = unp_res
            incr += 1

    
    for key in final_map.keys():
        entry_id, entity_id = key

        mappings = []

        for mapping in final_map[key]:
            (chains, pdb_seq, unp_seq, pdb_start, pdb_end, unp_start, unp_end) = mapping
            mappings.append({
                'chains': chains,
                'pdb_sequence': pdb_seq,
                'unp_sequence': unp_seq,
                'pdb_start': pdb_start,
                'pdb_end': pdb_end,
                'unp_start': unp_start,
                'unp_end': unp_end
            })

        final_result[uniprot_accession]['mappings'].append({    
            'entry_id': entry_id,
            'entity_id': entity_id,
            'segments': mappings
        })

    return jsonify(final_result)


def get_mappings_for_unp_residue_pfam(uniprot_accession):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)-[:IS_IN_PFAM]->(pfam:Pfam)
    WITH pfam.PFAM_ACCESSION AS pfamAccession, pfam.NAME AS pfamName, pfam.DESCRIPTION AS pfamDesc, toInteger(unp_res.ID) AS unpRes
    RETURN pfamAccession, pfamName, pfamDesc, min(unpRes) AS unpStart, max(unpRes) AS unpEnd
    """

    result = list(graph.run(query, parameters= {
        'accession': str(uniprot_accession)
    }))

    final_result = {}

    for pfam in result:
        
        (accession, name, desc, unp_start, unp_end) = pfam
        final_result[accession] = {
                "identifier": name,
                "name": name,
                "description": desc,
                "unp_start": unp_start,
                "unp_end": unp_end
            }

    return final_result


@app.route('/api/mappings/unipdb/<string:uniprot_accession>/<string:unp_res>')
def get_unipdb_residue(uniprot_accession, unp_res):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue {ID:$unpResidue})<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry)
    RETURN entry.ID AS entryId, entity.ID AS entityId, pdb_res.ID AS pdbRes
    """

    mappings = list(graph.run(query, parameters= {
        'accession': str(uniprot_accession), 'unpResidue': str(unp_res)
    }))

    final_result = {
        uniprot_accession: {
        }
    }

    pfam_dict = {}
    interpro_dict = {}
    cath_dict = {}
    scop_dict = {}

    for mapping in mappings:
        
        (entry_id, entity_id, pdb_res) = mapping
        temp_map = get_mappings_for_residue_pfam(entry_id, entity_id, pdb_res, graph)

        for pfam in temp_map.keys():
            if pfam_dict.get(pfam) is None:
                pfam_dict[pfam] = temp_map[pfam]

        temp_map = get_mappings_for_residue_interpro(entry_id, entity_id, pdb_res, graph)

        for interpro in temp_map.keys():
            if interpro_dict.get(interpro) is None:
                interpro_dict[interpro] = temp_map[interpro]

        temp_map = get_mappings_for_residue_cath(entry_id, entity_id, pdb_res, graph)

        for cath in temp_map.keys():
            if cath_dict.get(cath) is None:
                cath_dict[cath] = temp_map[cath]

        temp_map = get_mappings_for_residue_scop(entry_id, entity_id, pdb_res, graph)

        for scop in temp_map.keys():
            if scop_dict.get(scop) is None:
                scop_dict[scop] = temp_map[scop]

    final_result[uniprot_accession]['Pfam'] = pfam_dict
    final_result[uniprot_accession]['InterPro'] = interpro_dict
    final_result[uniprot_accession]['CATH'] = cath_dict
    final_result[uniprot_accession]['SCOP'] = scop_dict

    return jsonify(final_result)


@app.route('/api/pdb/compound/atoms/<string:chem_comp_id>')
def get_compound_atoms_api(chem_comp_id):

    return jsonify({
        chem_comp_id: get_compound_atoms(chem_comp_id, graph)
    })

@app.route('/api/pdb/compound/bonds/<string:chem_comp_id>')
def get_compound_bonds_api(chem_comp_id):

    return jsonify({
        chem_comp_id: get_compound_bonds(chem_comp_id, graph)
    })


@app.route('/api/pdb/compound/in_pdb/<string:chem_comp_id>')
def get_compound_in_pdb_api(chem_comp_id):

    return jsonify({
        chem_comp_id: get_compound_in_pdb(chem_comp_id, graph)
    })


@app.route('/api/validation/protein-ramachandran-sidechain-outliers/entry/<string:entry_id>')
def get_validation_protein_ramachandran_sidechain_outliers_api(entry_id):

    return jsonify({
        entry_id: get_validation_protein_ramachandran_sidechain_outliers(entry_id, graph)
    })