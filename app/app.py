from flask import Flask, jsonify
from .model import graph
from .amino_acid_codes import amino_acid_codes

app = Flask(__name__)

@app.route('/')
def default():

    return 'Default'

@app.route('/api/mappings/uniprot/<string:entry_id>')
def get_uniprot(entry_id):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[r:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt)
    RETURN toInteger(entity.ID) as entityId, unp.ACCESSION, unp.NAME, toInteger(pdb_res.ID) as pdbResId, toInteger(unp_res.ID) as unpResId, r.CHAINS order by toInteger(pdb_res.ID)
    """

    unp_map = {}
    unp_desc = {}
    result = list(graph.run(query, entry_id=entry_id))

    #print(result)

    prev = None
   
    incr = 0
    while incr < len(result):
        
        r = result[incr]

        (accession, name, pdb_res, unp_res, entity_id, chains) = (r['unp.ACCESSION'], r['unp.NAME'], r['pdbResId'], r['unpResId'], r['entityId'], r['r.CHAINS'])

        # query returns list of chains as string, making a list of chains
        chains = chains.translate({ord(c):'' for c in "[]' "}).split(',')

        rec = (pdb_res, unp_res, entity_id, chains)
        
        current = accession

        if prev != current:
            if unp_map.get(accession) is not None:
                unp_map[accession].append(rec)
            else:
                unp_map[accession] = [rec]
                unp_desc[accession] = name
            if incr != 0:
                r_prev = result[incr - 1]
                (prev_accession, prev_name, prev_pdb_res, prev_unp_res, prev_entity_id, prev_chains) = (r_prev['unp.ACCESSION'], r_prev['unp.NAME'], r_prev['pdbResId'], r_prev['unpResId'], r_prev['entityId'], r_prev['r.CHAINS'])
                prev_rec = (prev_pdb_res, prev_unp_res, prev_entity_id, prev_chains)
                unp_map[prev].append(prev_rec)
        # last record
        elif incr == len(result) - 1:
            unp_map[current].append(rec)

        prev = current
        incr += 1

    api_result = {
        entry_id: {
            "UniProt": {
                
            }
        }
    }

    for accession in unp_map.keys():

        api_result[entry_id]['UniProt'][accession] = {
            "identifier": unp_desc[accession],
            "name": unp_desc[accession],
            "mappings": [
            ]
        }

        incr = 0
        while incr < len(unp_map[accession]) - 1:
            mapping = unp_map[accession][incr]
            end_mapping = unp_map[accession][incr + 1]
            unp_end = end_mapping[1]
            end = end_mapping[0]
            del end_mapping

            for chain in list(mapping[3]):
              
                api_result[entry_id]['UniProt'][accession]["mappings"].append({
                    "entity_id": mapping[2],
                    "end": {
                        "author_residue_number": "",
                        "author_insertion_code": "",
                        "residue_number": end
                    },
                    "chain_id": chain,
                    "start": {
                        "author_residue_number": "",
                        "author_insertion_code": "",
                        "residue_number": mapping[0]
                    },
                    "unp_start": mapping[1],
                    "unp_end": unp_end,
                    "struct_asym_id": chain
                })

            incr += 2
    
    return jsonify(api_result)



@app.route('/api/mappings/residue_mapping/<string:entry_id>/<string:entity_id>/<string:residue_number>')
def get_mappings_for_residue(entry_id, entity_id, residue_number):

    final_result = { entry_id: {

    } }

    final_result[entry_id]["UniProt"] = get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number)
    final_result[entry_id]["Pfam"] = get_mappings_for_residue_pfam(entry_id, entity_id, residue_number)
    final_result[entry_id]["InterPro"] = get_mappings_for_residue_interpro(entry_id, entity_id, residue_number)
    final_result[entry_id]["CATH"] = get_mappings_for_residue_cath(entry_id, entity_id, residue_number)
    final_result[entry_id]["SCOP"] = get_mappings_for_residue_scop(entry_id, entity_id, residue_number)
    final_result[entry_id]["binding_sites"] = get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number)

    return jsonify(final_result)


def get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt)
    RETURN unp.ACCESSION as unp_accession, unp.NAME as unp_name
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    final_result = {}

    for unp in result:
        
        final_result[unp['unp_accession']] = {
                "identifier": unp['unp_name'],
                "name": unp['unp_name']
            }

    return final_result
    

def get_mappings_for_residue_pfam(entry_id, entity_id, residue_number):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)-[:IS_IN_PFAM]->(pfam:Pfam)
    RETURN pfam.ACCESSION as pfam_accession, pfam.NAME as pfam_name, pfam.DESCRIPTION as pfam_desc
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    final_result = {}

    for pfam in result:
        
        final_result[pfam['pfam_accession']] = {
                "identifier": pfam['pfam_name'],
                "name": pfam['pfam_name'],
                "description": pfam['pfam_desc']
            }

    return final_result

def get_mappings_for_residue_interpro(entry_id, entity_id, residue_number):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:IS_IN_INTERPRO]->(interpro:Interpro)
    RETURN interpro.ACCESSION as interpro_accession, interpro.NAME as interpro_name
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    final_result = {}

    for interpro in result:
        final_result[interpro['interpro_accession']] = {
                "identifier": interpro['interpro_name'],
                "name": interpro['interpro_name']
            }
        
    return final_result

def get_mappings_for_residue_cath(entry_id, entity_id, residue_number):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:IS_IN_CATH_DOMAIN]->(cath:CATH)
    RETURN cath.ARCH as arch, cath.CATHCODE as cathcode, cath.CLASS as class, cath.DOMAIN as domain, cath.HOMOL as homol, cath.TOPOL as topol, cath.NAME as name
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    final_result = {}

    for cath in result:
        final_result[cath['cathcode']] = {
                "homology": cath['homol'],
                "topology": cath['topol'],
                "architecture": cath['arch'],
                "identifier": cath['topol'],
                "class": cath['class'],
                "name": cath['name']
            }
        
    return final_result

def get_mappings_for_residue_scop(entry_id, entity_id, residue_number):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[relation:IS_IN_SCOP_DOMAIN]->(scop:SCOP)
    OPTIONAL MATCH (superfamily:SCOP {SUNID: relation.SUPERFAMILY_ID})
    OPTIONAL MATCH (fold:SCOP {SUNID: relation.FOLD_ID})
    OPTIONAL MATCH (class:SCOP {SUNID: relation.CLASS_ID})
    RETURN  distinct scop.SUNID as sunid, scop.DESCRIPTION as desc, superfamily.SUNID as super_sunid, superfamily.DESCRIPTION as super_desc, 
    fold.SUNID as fold_sunid, fold.DESCRIPTION as fold_desc, class.SUNID as class_sunid, class.DESCRIPTION as class_desc, scop.SCCS as sccs, 
    relation.SCOP_ID as scop_id, relation.AUTH_ASYM_ID as chain_id
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    final_result = {}
    sunid_mapping = {}

    for scop in result:
        
        if sunid_mapping.get(scop['sunid']) is None:
            sunid_mapping[scop['sunid']] = []

        sunid_mapping[scop['sunid']].append((scop['scop_id'], scop['chain_id']))

        final_result[scop['sunid']] = {
                "superfamily": {
                    "sunid": scop['super_sunid'],
                    "description": scop['super_desc']
                },
                "sccs": scop['sccs'],
                "fold": {
                    "sunid": scop['fold_sunid'],
                    "description": scop['fold_desc']
                },
                "identifier": scop['desc'],
                "class": {
                    "sunid": scop['class_sunid'],
                    "description": scop['class_desc']
                },
                "description": scop['desc'],
                "mappings": []
            }

    for sunid in sunid_mapping.keys():
        for mapping in sunid_mapping[sunid]:
            scop_id, chain_id = mapping

            final_result[sunid]['mappings'].append({
                "scop_id": scop_id,
                "chain_id": chain_id,
                "struct_asym_id": chain_id
            })

    return final_result



def get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {ID:$entity_id})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residue})-[res_relation:IS_IN_BINDING_SITE]->
    (site:Binding_Site)
    WITH site, pdb_res
    MATCH (ligand:Ligand)-[ligand_relation:IS_IN_BINDING_SITE]->(site)-[bound_relation:BOUNDED_BY]->(boundLigand:Ligand)
    RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, boundLigand.ID, boundLigand.NAME, boundLigand.FORMULA, bound_relation.AUTH_ASYM_ID, bound_relation.STRUCT_ASYM_ID, 
    bound_relation.AUTH_SEQ_ID, bound_relation.ENTITY_ID, bound_relation.RESIDUE_ID, ligand.ID, ligand.NAME, ligand.FORMULA, ligand_relation.AUTH_ASYM_ID, ligand_relation.AUTH_SEQ_ID, 
    ligand_relation.STRUCT_ASYM_ID, ligand_relation.SYMMETRY_SYMBOL, ligand_relation.ENTITY_ID, ligand_relation.RESIDUE_ID, pdb_res.ID, pdb_res.CHEM_COMP_ID ORDER BY site.ID
    """

    mappings = list(graph.run(query, parameters= {
        'entry_id': str(entry_id), 'entity_id': str(entity_id), 'residue': str(residue_number)
    }))

    site_dict = {}
    site_ligand_dict = {}
    site_boundligand_dict = {}
    residue_chem_comp_id = None
    final_result = []

    for mapping in mappings:
        
        (site_id, site_name, site_evidence, bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, bound_auth_seq_id, bound_entity_id,
        bound_residue_id, ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id,
        pdb_res, pdb_res_chem_comp_id) = mapping

        if(residue_chem_comp_id is None):
            residue_chem_comp_id = pdb_res_chem_comp_id

        if(site_dict.get(site_id) is None):
            site_dict[site_id] = (site_name, site_evidence)

        if(site_boundligand_dict.get(site_id) is None):
            site_boundligand_dict[site_id] = [(bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, 
                                                              bound_auth_seq_id, bound_entity_id, bound_residue_id)]
        else:
            site_boundligand_dict[site_id].append((bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, 
                                                              bound_auth_seq_id, bound_entity_id, bound_residue_id))

        if(site_ligand_dict.get(site_id) is None):
            site_ligand_dict[site_id] = [(ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id)]
        else:
            site_ligand_dict[site_id].append((ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id))


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
        for result in site_ligand_dict[key]:
            (ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id) = result
            
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

    return final_result

@app.route('/api/mappings/binding_sites/<string:entry_id>')
def get_binding_sites_for_entry(entry_id):

    site_dict = {}
    site_ligand_dict = {}
    site_boundligand_dict = {}
    site_residue_dict = {}
    final_result = []

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[relation:HAS_BINDING_SITE]->(site:Binding_Site)
    OPTIONAL MATCH (ligand:Ligand)-[ligand_relation:IS_IN_BINDING_SITE]->(site)
    OPTIONAL MATCH (site)-[bound_relation:BOUNDED_BY]->(boundLigand:Ligand)
    RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, boundLigand.ID, boundLigand.NAME, boundLigand.FORMULA, bound_relation.AUTH_ASYM_ID, bound_relation.STRUCT_ASYM_ID, 
    bound_relation.AUTH_SEQ_ID, bound_relation.ENTITY_ID, bound_relation.RESIDUE_ID, ligand.ID, ligand.NAME, ligand.FORMULA, ligand_relation.AUTH_ASYM_ID, ligand_relation.AUTH_SEQ_ID, 
    ligand_relation.STRUCT_ASYM_ID, ligand_relation.SYMMETRY_SYMBOL, ligand_relation.ENTITY_ID, ligand_relation.RESIDUE_ID ORDER BY site.ID
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
    WITH site, pdb_res, res_relation
    RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, pdb_res.ID, res_relation.ENTITY_ID, pdb_res.CHEM_COMP_ID, res_relation.AUTH_ASYM_ID, 
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
        print(unp_res_id, unp_res_code, pdb_res_id, pdb_res_code, entity_id, entry_id)

        residue_result = get_mappings_for_residue_binding_site(entry_id, entity_id, pdb_res_id)

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
    WITH pfam.ACCESSION AS pfamAccession, pfam.NAME AS pfamName, pfam.DESCRIPTION AS pfamDesc, toInteger(unp_res.ID) AS unpRes
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
        temp_map = get_mappings_for_residue_pfam(entry_id, entity_id, pdb_res)

        for pfam in temp_map.keys():
            if pfam_dict.get(pfam) is None:
                pfam_dict[pfam] = temp_map[pfam]

        temp_map = get_mappings_for_residue_interpro(entry_id, entity_id, pdb_res)

        for interpro in temp_map.keys():
            if interpro_dict.get(interpro) is None:
                interpro_dict[interpro] = temp_map[interpro]

        temp_map = get_mappings_for_residue_cath(entry_id, entity_id, pdb_res)

        for cath in temp_map.keys():
            if cath_dict.get(cath) is None:
                cath_dict[cath] = temp_map[cath]

        temp_map = get_mappings_for_residue_scop(entry_id, entity_id, pdb_res)

        for scop in temp_map.keys():
            if scop_dict.get(scop) is None:
                scop_dict[scop] = temp_map[scop]

    final_result[uniprot_accession]['Pfam'] = pfam_dict
    final_result[uniprot_accession]['InterPro'] = interpro_dict
    final_result[uniprot_accession]['CATH'] = cath_dict
    final_result[uniprot_accession]['SCOP'] = scop_dict

    return jsonify(final_result)

