#!/usr/bin/env python

"""
residue.py: This have list of function calls which retrieves data related to a PDB or UniProt residue
"""

__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"

def get_mappings_for_residue_uniprot(entry_id, entity_id, residue_number, graph):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt)
    RETURN unp.ACCESSION as unp_accession, unp.NAME as unp_name
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    if(len(result) == 0):
        return {}, 404

    final_result = {}

    for unp in result:
        
        final_result[unp['unp_accession']] = {
                "identifier": unp['unp_name'],
                "name": unp['unp_name']
            }

    return final_result, 200

def get_mappings_for_residue_pfam(entry_id, entity_id, residue_number, graph):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)-[:IS_IN_PFAM]->(pfam:Pfam)
    RETURN pfam.PFAM_ACCESSION as pfam_accession, pfam.NAME as pfam_name, pfam.DESCRIPTION as pfam_desc
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    if(len(result) == 0):
        return {}, 404

    final_result = {}

    for pfam in result:
        
        final_result[pfam['pfam_accession']] = {
                "identifier": pfam['pfam_name'],
                "name": pfam['pfam_name'],
                "description": pfam['pfam_desc']
            }

    return final_result, 200

def get_mappings_for_residue_interpro(entry_id, entity_id, residue_number, graph):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:IS_IN_INTERPRO]->(interpro:Interpro)
    RETURN interpro.INTERPRO_ACCESSION as interpro_accession, interpro.NAME as interpro_name
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    if(len(result) == 0):
        return {}, 404

    final_result = {}

    for interpro in result:
        final_result[interpro['interpro_accession']] = {
                "identifier": interpro['interpro_name'],
                "name": interpro['interpro_name']
            }
        
    return final_result, 200

def get_mappings_for_residue_cath(entry_id, entity_id, residue_number, graph):

    query = """
    MATCH (entry:Entry {ID:$entryId})-[:HAS_ENTITY]->(entity:Entity {ID:$entityId})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residueNumber})-[:IS_IN_CATH_DOMAIN]->(cath:CATH)
    RETURN cath.ARCH as arch, cath.CATHCODE as cathcode, cath.CLASS as class, cath.DOMAIN as domain, cath.HOMOL as homol, cath.TOPOL as topol, cath.NAME as name
    """

    result = list(graph.run(query, parameters= {
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': residue_number
    }))

    if(len(result) == 0):
        return {}, 404

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
        
    return final_result, 200

def get_mappings_for_residue_scop(entry_id, entity_id, residue_number, graph):

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
        'entryId': str(entry_id), 'entityId': str(entity_id), 'residueNumber': str(residue_number)
    }))

    if(len(result) == 0):
        return {}, 404

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

    return final_result, 200


def get_mappings_for_residue_binding_site(entry_id, entity_id, residue_number, site_residues, graph):

    query = None
    if(site_residues is True):
        query = """
        MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {ID:$entity_id})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residue})-[res_relation:IS_IN_BINDING_SITE]->
        (site:Binding_Site)
        WITH site, pdb_res
        OPTIONAL MATCH (ligand:Ligand)-[:IS_AN_INSTANCE_OF]->(ligand_entity:Entity)-[ligand_entity_relation:CONTAINS_CHAIN]->(ligand_chain:Chain)-[ligand_relation:IS_IN_BINDING_SITE]->(site)
        OPTIONAL MATCH (site)-[bound_relation:BOUNDED_BY]->(boundligand_chain:Chain)<-[boundligand_entity_relation:CONTAINS_CHAIN]-(boundligand_entity:Entity)-[:IS_AN_INSTANCE_OF]->(boundligand:Ligand)
        OPTIONAL MATCH (site)<-[res_all_relation:IS_IN_BINDING_SITE]-(pdb_res_all:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(pdb_res_all_entity:Entity) WHERE pdb_res_all.ID <> $residue
        RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, pdb_res_all.ID, pdb_res_all.CHEM_COMP_ID, res_all_relation.AUTH_ASYM_ID, res_all_relation.STRUCT_ASYM_ID, 
        res_all_relation.AUTH_SEQ_ID, pdb_res_all_entity.ID, res_all_relation.SYMMETRY_SYMBOL, boundligand.ID, boundligand.NAME, boundligand.FORMULA, boundligand_chain.AUTH_ASYM_ID, boundligand_chain.STRUCT_ASYM_ID, 
        boundligand_entity_relation.AUTH_SEQ_ID, boundligand_entity.ID, boundligand_entity_relation.RES_ID, ligand.ID, ligand.NAME, ligand.FORMULA, ligand_chain.AUTH_ASYM_ID, ligand_chain.AUTH_SEQ_ID, 
        ligand_chain.STRUCT_ASYM_ID, ligand_relation.SYMMETRY_SYMBOL, ligand_entity.ID, ligand_chain.RES_ID, pdb_res.ID, pdb_res.CHEM_COMP_ID ORDER BY site.ID
        """
    else:
        query = """
        MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {ID:$entity_id})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residue})-[res_relation:IS_IN_BINDING_SITE]->
        (site:Binding_Site)
        WITH site, pdb_res
        MATCH (site)-[bound_relation:BOUNDED_BY]->(boundChain:Chain)<-[entity_chain_rel:CONTAINS_CHAIN]-(entity:Entity)-[entity_lig_rel:IS_AN_INSTANCE_OF]->(ligand:Ligand)
        RETURN site.ID, site.DETAILS, site.EVIDENCE_CODE, ligand.ID, ligand.NAME, ligand.FORMULA, boundChain.AUTH_ASYM_ID, boundChain.STRUCT_ASYM_ID, 
        entity_chain_rel.AUTH_SEQ_ID, entity.ID, entity_chain_rel.RES_ID, pdb_res.ID, pdb_res.CHEM_COMP_ID ORDER BY site.ID
        """


    mappings = list(graph.run(query, parameters= {
        'entry_id': str(entry_id), 'entity_id': str(entity_id), 'residue': str(residue_number)
    }))

    if(len(mappings) == 0):
        return {}, 404

    site_dict = {}
    site_ligand_dict = {}
    site_boundligand_dict = {}
    site_pdb_res_dict = {}
    residue_chem_comp_id = None
    final_result = []

    for mapping in mappings:
        
        if(site_residues is True):
            (site_id, site_name, site_evidence, pdb_res_id, pdb_res_chem_comp_id, pdb_res_auth_asym_id, pdb_res_struct_asym_id, pdb_res_auth_seq_id, pdb_res_entity_id, pdb_res_symmetry_symbol, bound_ligand, 
            bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, bound_auth_seq_id, bound_entity_id, bound_residue_id, ligand, ligand_name, ligand_formula, 
            ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id, pdb_res, pdb_res_chem_comp_id) = mapping
        else:
            (site_id, site_name, site_evidence, bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, 
            bound_auth_seq_id, bound_entity_id, bound_residue_id, pdb_res, pdb_res_chem_comp_id) = mapping

        if(residue_chem_comp_id is None):
            residue_chem_comp_id = pdb_res_chem_comp_id

        if(site_dict.get(site_id) is None):
            site_dict[site_id] = (site_name, site_evidence)

        if(bound_ligand is not None):
            if(site_boundligand_dict.get(site_id) is None):
                site_boundligand_dict[site_id] = [(bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, 
                                                                bound_auth_seq_id, bound_entity_id, bound_residue_id)]
            else:
                site_boundligand_dict[site_id].append((bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, 
                                                                bound_auth_seq_id, bound_entity_id, bound_residue_id))
        if(site_residues is True):
            if(ligand is not None):
                if(site_ligand_dict.get(site_id) is None):
                    site_ligand_dict[site_id] = [(ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id)]
                else:
                    site_ligand_dict[site_id].append((ligand, ligand_name, ligand_formula, ligand_auth_asym_id, ligand_auth_seq_id, ligand_struct_asym_id, ligand_sym_symbol, ligand_entity_id, ligand_residue_id))

        if(site_residues is True):
            if(pdb_res_id is not None):
                if(site_pdb_res_dict.get(site_id) is None):
                    site_pdb_res_dict[site_id] = [(pdb_res_id, pdb_res_chem_comp_id, pdb_res_auth_asym_id, pdb_res_struct_asym_id, pdb_res_auth_seq_id, pdb_res_entity_id, pdb_res_symmetry_symbol)]
                else:
                    site_pdb_res_dict[site_id].append((pdb_res_id, pdb_res_chem_comp_id, pdb_res_auth_asym_id, pdb_res_struct_asym_id, pdb_res_auth_seq_id, pdb_res_entity_id, pdb_res_symmetry_symbol))

    for key in site_dict.keys():
        (site_name, evidence) = site_dict[key]
        if(site_residues is True):
            temp = {
                "site_id": key,
                "evidence_code": evidence,
                "details": site_name,
                "site_residues": [],
                "ligand_residues": []
            }
        else:
            temp = {
                "site_id": key,
                "evidence_code": evidence,
                "details": site_name,
                "ligand_residues": []
            }

        if(site_boundligand_dict.get(key) is not None):
            for result in list(set(site_boundligand_dict[key])):
                (bound_ligand, bound_ligand_name, bound_ligand_formula, bound_auth_asym_id, bound_struct_asym_id, bound_auth_seq_id, bound_entity_id, bound_residue_id) = result
                if(bound_residue_id is not None):
                    temp["ligand_residues"].append({
                        "entity_id": int(bound_entity_id),
                        "residue_number": int(bound_residue_id),
                        "author_insertion_code": "null",
                        "chain_id": bound_auth_asym_id,
                        "author_residue_number": int(bound_auth_seq_id),
                        "chem_comp_id": bound_ligand,
                        "struct_asym_id": bound_struct_asym_id
                    })
        
        if(site_residues is True):
            if(site_ligand_dict.get(key) is not None):
                for result in list(set(site_ligand_dict[key])):
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

        if(site_residues is True):
            if(site_pdb_res_dict.get(key) is not None):
                for result in list(set(site_pdb_res_dict[key])):
                    (pdb_res_id, pdb_res_chem_comp_id, pdb_res_auth_asym_id, pdb_res_struct_asym_id, pdb_res_auth_seq_id, pdb_res_entity_id, pdb_res_symmetry_symbol) = result

                    temp["site_residues"].append({
                        "entity_id": int(pdb_res_entity_id),
                        "residue_number": int(pdb_res_id),
                        "author_insertion_code": "null",
                        "chain_id": pdb_res_auth_asym_id,
                        "author_residue_number": int(pdb_res_auth_seq_id),
                        "chem_comp_id": pdb_res_chem_comp_id,
                        "struct_asym_id": pdb_res_struct_asym_id,
                        "symmetry_symbol": pdb_res_symmetry_symbol
                    })


        final_result.append(temp)

    return final_result, 200


def get_basic_residue_details(entry_id, entity_id, residue_number, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {ID:$entity_id})-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue {ID:$residue_number})-[chain_rel:IS_IN_CHAIN]->(chain:Chain)
    RETURN chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, chain_rel.AUTH_SEQ_ID, chain_rel.PDB_INS_CODE, chain_rel.OBSERVED
    """

    chains = []
    dict_chains = {}

    mappings = list(graph.run(query, parameters= {
        'entry_id': str(entry_id), 'entity_id': str(entity_id), 'residue_number': str(residue_number)
    }))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:
        (auth_asym_id, struct_asym_id, auth_seq_id, pdb_ins_code, observed) = mapping
        chain_key = (auth_asym_id, struct_asym_id)

        if dict_chains.get(chain_key) is None:
            dict_chains[chain_key] = {
                "auth_asym_id": auth_asym_id,
                "struct_asym_id": struct_asym_id,
                "residues": []
            }
        dict_chains[chain_key]["residues"].append({
            "residue_number": int(residue_number),
            "author_residue_number": int(auth_seq_id) if auth_seq_id is not None else None,
            "author_insertion_code": "" if pdb_ins_code is None else pdb_ins_code,
            "observed_in_chain": observed
        })

        
    for chain_key in dict_chains.keys():
        chains.append(dict_chains[chain_key])

    return chains, 200