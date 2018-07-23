
import more_itertools as mit
from .amino_acid_codes import amino_acid_codes


def get_uniprot(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[r:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt {TYPE:'CANONICAL'}),
    (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain)
    RETURN toInteger(entity.ID) as entityId, unp.ACCESSION, unp.NAME, toInteger(pdb_res.ID) as pdbResId, toInteger(unp_res.ID) as unpResId, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID,
    toInteger(chain_rel.AUTH_SEQ_ID) as auth_seq_id order by toInteger(pdb_res.ID)
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    api_result = {}
    list_of_accessions = []
    dict_unp_master = {}
    dict_residues = {}
    dict_chains = {}
    dict_entity = {}
    dict_residues_mappings = {}
    dict_struct_asym_id = {}

    for mapping in mappings:

        (entity_id, accession, name, pdb_res_id, unp_res_id, auth_asym_id, struct_asym_id, auth_seq_id) = mapping

        if(accession not in list_of_accessions):
            list_of_accessions.append(accession)
            dict_unp_master[accession] = name

        if(dict_entity.get(accession) is None):
            dict_entity[accession] = [entity_id]
        else:
            dict_entity[accession].append(entity_id)

        if(dict_chains.get(entity_id) is None):
            dict_chains[entity_id] = [auth_asym_id]
        else:
            dict_chains[entity_id].append(auth_asym_id)

        if(dict_struct_asym_id.get((entity_id, auth_asym_id)) is None):
            dict_struct_asym_id[(entity_id, auth_asym_id)] = struct_asym_id

        if(dict_residues.get((entity_id, auth_asym_id)) is None):
            dict_residues[(entity_id, auth_asym_id)] = [pdb_res_id]
        else:
            dict_residues[(entity_id, auth_asym_id)].append(pdb_res_id)

        dict_residues_mappings[(entity_id, auth_asym_id, pdb_res_id)] = (unp_res_id, auth_seq_id)
    
    for accession in list_of_accessions:
        
        temp_mappings = []
        for entity_id in set(dict_entity[accession]):
            for auth_asym_id in set(dict_chains[entity_id]):
                for group in mit.consecutive_groups(dict_residues[(entity_id, auth_asym_id)]):
                    group = list(group)
                    start_auth_seq_id = dict_residues_mappings[(entity_id, auth_asym_id, group[0])][1]
                    end_auth_seq_id = dict_residues_mappings[(entity_id, auth_asym_id, group[-1])][1]
                    element = {
                        "entity_id": int(entity_id),
                        "chain_id": auth_asym_id,
                        "struct_asym_id": dict_struct_asym_id[(entity_id, auth_asym_id)],
                        "start": {
                            "author_residue_number": None if start_auth_seq_id is None else int(start_auth_seq_id),
                            "author_insertion_code": "",
                            "residue_number": int(group[0])
                        },
                        "end": {
                            "author_residue_number": None if end_auth_seq_id is None else int(end_auth_seq_id),
                            "author_insertion_code": "",
                            "residue_number": int(group[-1])
                        },
                        "unp_start": int(dict_residues_mappings[(entity_id, auth_asym_id, group[0])][0]),
                        "unp_end": int(dict_residues_mappings[(entity_id, auth_asym_id, group[-1])][0])
                    }
                    temp_mappings.append(element)
        
        api_result[accession] = {
            "name": dict_unp_master[accession],
            "identifier": dict_unp_master[accession],
            "mappings": temp_mappings
        }
        del temp_mappings
                    
    return api_result, 200


def get_uniprot_segments(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[r:MAP_TO_UNIPROT_RESIDUE]->
    (unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt {TYPE:'CANONICAL'}),
    (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain)
    RETURN toInteger(entity.ID) as entityId, unp.ACCESSION, unp.NAME, toInteger(pdb_res.ID) as pdbResId, toInteger(unp_res.ID) as unpResId, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID,
    toInteger(chain_rel.AUTH_SEQ_ID) as auth_seq_id order by toInteger(unp_res.ID)
    """

    dict_unp = {}
    dict_pdb_res = {}
    dict_unp_res = {}
    dict_auth_seq_id = {}
    dict_auth_seq_id_nulls = {}
    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:
        (entity_id, unp_accession, unp_name, pdb_res_id, unp_res_id,
         auth_asym_id, struct_asym_id, auth_seq_id) = mapping

        # ugly fix, need to see for any issues in data consistency - auth_seq_id for unobserved residues is null which fails mit function
        if(auth_seq_id is None):
            auth_seq_id = pdb_res_id
            dict_auth_seq_id_nulls[(entity_id, unp_accession, unp_name, auth_asym_id, struct_asym_id, auth_seq_id)] = True

        key = (entity_id, unp_accession, unp_name,
               auth_asym_id, struct_asym_id)

        if(dict_unp.get(unp_accession) is None):
            dict_unp[unp_accession] = {
                "identifier": unp_name,
                "name": unp_name,
                "mappings": []
            }

        if(dict_pdb_res.get(key) is None):
            dict_pdb_res[key] = [pdb_res_id]
            dict_unp_res[key] = [unp_res_id]
            dict_auth_seq_id[key] = [auth_seq_id]
        else:
            dict_pdb_res[key].append(pdb_res_id)
            dict_unp_res[key].append(unp_res_id)
            dict_auth_seq_id[key].append(auth_seq_id)


    for key in dict_pdb_res.keys():

        (entity_id, unp_accession, unp_name, auth_asym_id, struct_asym_id) = key

        temp_list = []
        for group in mit.consecutive_groups(dict_pdb_res[key]):
            group = list(group)
            temp_list.append((group[0], group[-1]))
        dict_pdb_res[key] = temp_list

        # dont clear since it still holds the reference, so create new
        del temp_list
        temp_list = []

        for group in mit.consecutive_groups(dict_unp_res[key]):
            group = list(group)
            temp_list.append((group[0], group[-1]))
        dict_unp_res[key] = temp_list

        # dont clear since it still holds the reference, so create new
        del temp_list
        temp_list = []
        for group in mit.consecutive_groups(dict_auth_seq_id[key]):
            group = list(group)
            temp_list.append((group[0], group[-1]))
        dict_auth_seq_id[key] = temp_list

        # dont clear since it still holds the reference, so create new
        del temp_list
        temp_list = []

    api_result = {}
    
    for key in dict_pdb_res.keys():

        (entity_id, unp_accession, unp_name, auth_asym_id, struct_asym_id) = key

        api_result[unp_accession] = dict_unp[unp_accession]

        incr = 0

        while(incr < len(dict_pdb_res[key])):

            print('start->', dict_auth_seq_id[key][incr][0])
            auth_start_null = dict_auth_seq_id_nulls.get((entity_id, unp_accession, unp_name, auth_asym_id, struct_asym_id, dict_auth_seq_id[key][incr][0]))
            auth_end_null = dict_auth_seq_id_nulls.get((entity_id, unp_accession, unp_name, auth_asym_id, struct_asym_id, dict_auth_seq_id[key][incr][-1]))

            api_result[unp_accession]["mappings"].append({
                "start": {
                    "author_residue_number": None if auth_start_null is True else dict_auth_seq_id[key][incr][0],
                    "residue_number": dict_pdb_res[key][incr][0],
                    "author_insertion_code": ""
                },
                "end": {
                    "author_residue_number": None if auth_end_null is True else dict_auth_seq_id[key][incr][-1],
                    "residue_number": dict_pdb_res[key][incr][-1],
                    "author_insertion_code": ""
                },
                "unp_start": dict_unp_res[key][incr][0],
                "unp_end": dict_unp_res[key][incr][-1],
                "entity_id": entity_id,
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id
            })

            incr += 1

    return api_result, 200


def get_pfam(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)-[:IS_IN_PFAM]->(pfam:Pfam), 
    (pdb_res)-[chain_relation:IS_IN_CHAIN]->(chain:Chain)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, pdb_res.ID, chain_relation.AUTH_SEQ_ID, pfam.PFAM_ACCESSION, pfam.DESCRIPTION 
    ORDER BY chain.AUTH_ASYM_ID, toInteger(pdb_res.ID)
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    dict_pfam = {}
    dict_pfam_mappings = {}
    api_result = {}

    for mapping in mappings:

        (entity_id, auth_asym_id, struct_asym_id, res_id,
         auth_seq_id, pfam_accession, desc) = mapping

        if(dict_pfam.get(pfam_accession) is None):

            api_result[pfam_accession] = {
                "identifier": desc,
                "description": desc,
                "mappings": []
            }

            dict_pfam[pfam_accession] = desc

        key = (entity_id, auth_asym_id, struct_asym_id, pfam_accession)

        if(dict_pfam_mappings.get(key) is None):
            dict_pfam_mappings[key] = [(res_id, auth_seq_id)]
        else:
            dict_pfam_mappings[key].append((res_id, auth_seq_id))

    for key in dict_pfam_mappings.keys():

        (entity_id, auth_asym_id, struct_asym_id, pfam_accession) = key
        (start_res_id, start_auth_seq_id) = dict_pfam_mappings[key][0]
        (end_res_id, end_auth_seq_id) = dict_pfam_mappings[key][-1]

        api_result[pfam_accession]["mappings"].append({
            "entity_id": int(entity_id),
            "chain_id": auth_asym_id,
            "struct_asym_id": struct_asym_id,
            "end": {
                "author_residue_number": None if end_auth_seq_id is None else int(end_auth_seq_id),
                "author_insertion_code": "",
                "residue_number": int(end_res_id)
            },
            "start": {
                "author_residue_number": None if start_auth_seq_id is None else int(start_auth_seq_id),
                "author_insertion_code": "",
                "residue_number": int(start_res_id)
            }
        })

    return api_result, 200


def get_interpro(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_INTERPRO]->(interpro:Interpro),
    (pdb_res)-[chain_rel:IS_IN_CHAIN {OBSERVED: 'Y'}]->(chain:Chain) WHERE chain.AUTH_ASYM_ID = rel.AUTH_ASYM_ID
    RETURN toInteger(entity.ID), toInteger(pdb_res.ID), toInteger(chain_rel.AUTH_SEQ_ID), chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, rel.GROUP_ID, 
    interpro.INTERPRO_ACCESSION, interpro.NAME ORDER BY toInteger(pdb_res.ID)
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    dict_interpro_res = {}
    dict_interpro_auth_seq = {}

    api_result = {}

    for mapping in mappings:
        (entity_id, res_id, auth_seq_id, auth_asym_id, struct_asym_id,
         group_id, interpro_accession, interpro_name) = mapping

        key = (interpro_accession, interpro_name, entity_id,
               auth_asym_id, struct_asym_id, group_id)
        if(dict_interpro_res.get(key) is None):
            dict_interpro_res[key] = [res_id]
            dict_interpro_auth_seq[key] = [auth_seq_id]

            api_result[interpro_accession] = {
                "identifier": interpro_name,
                "name": interpro_name,
                "mappings": []
            }
        else:
            dict_interpro_res[key].append(res_id)
            dict_interpro_auth_seq[key].append(auth_seq_id)

    for key in dict_interpro_res.keys():
        (interpro_accession, interpro_name, entity_id,
         auth_asym_id, struct_asym_id, group_id) = key

        res_list = dict_interpro_res[key]
        auth_list = dict_interpro_auth_seq[key]

        temp_segment = {
            "entity_id": entity_id,
            "end": {
                "author_residue_number": auth_list[-1],
                "author_insertion_code": "",
                "residue_number": res_list[-1]
            },
            "start": {
                "author_residue_number": auth_list[0],
                "author_insertion_code": "",
                "residue_number": res_list[0]
            },
            "chain_id": auth_asym_id,
            "struct_asym_id": struct_asym_id
        }

        api_result[interpro_accession]["mappings"].append(
            temp_segment)

    return api_result, 200


def get_cath(entry_id, graph):

    query = """
    MATCH(entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_CATH_DOMAIN]->(cath:CATH)
    RETURN entity.ID, toInteger(pdb_res.ID), cath.CATHCODE, cath.DOMAIN, cath.HOMOL, cath.CLASS, cath.ARCH, cath.TOPOL, cath.NAME, rel.AUTH_ASYM_ID, 
    rel.STRUCT_ASYM_ID, rel.AUTH_START, rel.AUTH_END, rel.SEGMENT ORDER by toInteger(pdb_res.ID)
    """

    list_cath_code = []
    dict_cath = {}
    dict_domain = {}
    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    api_result = {}

    for mapping in mappings:
        (entity_id, res_id, cathcode, cath_domain, cath_homology, cath_class, cath_arch, cath_topology, cath_name, auth_asym_id,
         struct_asym_id, auth_start, auth_end, segment) = mapping

        if(cathcode not in list_cath_code):
            list_cath_code.append(cathcode)

            api_result[cathcode] = {
                "homology": cath_homology,
                "name": cath_name,
                "class": cath_class,
                "architecture": cath_arch,
                "identifier": cath_topology,
                "topology": cath_topology,
                "mappings": []
            }

        if(dict_domain.get(cath_domain) is None):
            dict_domain[cath_domain] = (auth_start, auth_end)

        key = (entity_id, cathcode, cath_domain,
               auth_asym_id, struct_asym_id, segment)

        if(dict_cath.get(key) is None):
            dict_cath[key] = [res_id]
        else:
            dict_cath[key].append(res_id)

    for key in dict_cath.keys():
        (entity_id, cathcode, cath_domain,
         auth_asym_id, struct_asym_id, segment) = key
        (auth_start, auth_end) = dict_domain[cath_domain]

        temp_segment = {
            "domain": cath_domain,
            "segment_id": int(segment),
            "entity_id": int(entity_id),
            "chain_id": auth_asym_id,
            "struct_asym_id": struct_asym_id,
            "start": {
                "author_residue_number": int(auth_start),
                "author_insertion_code": "",
                "residue_number": dict_cath[key][0]
            },
            "end": {
                "author_residue_number": int(auth_end),
                "author_insertion_code": "",
                "residue_number": dict_cath[key][-1]
            }
        }

        api_result[cathcode]["mappings"].append(temp_segment)

    return api_result, 200


def get_scop(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[relation:IS_IN_SCOP_DOMAIN]->(scop:SCOP)
    OPTIONAL MATCH (superfamily:SCOP {SUNID: relation.SUPERFAMILY_ID})
    OPTIONAL MATCH (fold:SCOP {SUNID: relation.FOLD_ID})
    OPTIONAL MATCH (class:SCOP {SUNID: relation.CLASS_ID})
    RETURN scop.SUNID as sunid, scop.DESCRIPTION as desc, superfamily.SUNID as super_sunid, superfamily.DESCRIPTION as super_desc, 
    fold.SUNID as fold_sunid, fold.DESCRIPTION as fold_desc, class.SUNID as class_sunid, class.DESCRIPTION as class_desc, scop.SCCS as sccs, 
    relation.SCOP_ID as scop_id, pdb_res.ID, entity.ID, relation.AUTH_ASYM_ID, relation.STRUCT_ASYM_ID, relation.SEGMENT_ID, relation.AUTH_START, 
    relation.AUTH_END ORDER BY toInteger(pdb_res.ID)
    """

    api_result = {}

    list_sunid = []
    dict_scop = {}
    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404
        
    dict_scop_auth = {}
    dict_scop_sunid = {}

    for mapping in mappings:
        (sunid, desc, super_sunid, super_desc, fold_sunid, fold_desc, class_sunid, class_desc, sccs, scop_id, res_id, entity_id, auth_asym_id, struct_asym_id,
         segment_id, auth_start, auth_end) = mapping

        if(sunid not in list_sunid):

            api_result[sunid] = {
                "superfamily": {
                    "sunid": int(super_sunid),
                    "description": super_desc
                },
                "sccs": sccs,
                "fold": {
                    "sunid": int(fold_sunid),
                    "description": fold_desc
                },
                "identifier": desc,
                "class": {
                    "sunid": int(class_sunid),
                    "description": class_desc
                },
                "description": desc,
                "mappings": []
            }

            list_sunid.append(sunid)

        key = (entity_id, auth_asym_id, struct_asym_id, scop_id, segment_id)

        if(dict_scop.get(key) is None):
            dict_scop[key] = [res_id]
        else:
            dict_scop[key].append(res_id)

        if(dict_scop_sunid.get(scop_id) is None):
            dict_scop_sunid[scop_id] = sunid

        # keep auth_start and auth_end for a domain, this is not to be stored in main dictionary
        if(dict_scop_auth.get(scop_id) is None):
            dict_scop_auth[scop_id] = (auth_start, auth_end)

    for key in dict_scop.keys():
        (entity_id, auth_asym_id, struct_asym_id, scop_id, segment_id) = key

        (auth_start, auth_end) = dict_scop_auth[scop_id]

        temp_segment = {
            "entity_id": int(entity_id),
            "end": {
                "author_residue_number": None if auth_end is None else int(auth_end),
                "author_insertion_code": "",
                "residue_number": int(dict_scop[key][-1])
            },
            "segment_id": int(segment_id),
            "chain_id": auth_asym_id,
            "scop_id": scop_id,
            "start": {
                "author_residue_number": None if auth_start is None else int(auth_start),
                "author_insertion_code": "",
                "residue_number": int(dict_scop[key][0])
            },
            "struct_asym_id": struct_asym_id
        }

        api_result[dict_scop_sunid[scop_id]]["mappings"].append(temp_segment)

    return api_result, 200


def get_go(entry_id, graph):

    bio_go_query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]->(chain:Chain)
    MATCH (entity)-[:HAS_GO]->(bio_go:GO_Biological_Process)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, bio_go.GO_ID, bio_go.DEFINITION, bio_go.NAME
    """

    bio_go_mappings = list(graph.run(bio_go_query, entry_id=entry_id))

    cell_go_query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]->(chain:Chain)
    MATCH (entity)-[:HAS_GO]->(cell_go:GO_Cellular_Component)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, cell_go.GO_ID, cell_go.DEFINITION, cell_go.NAME
    """

    cell_go_mappings = list(graph.run(cell_go_query, entry_id=entry_id))

    mol_go_query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity {POLYMER_TYPE:'P'})-[:CONTAINS_CHAIN]->(chain:Chain)
    MATCH (entity)-[:HAS_GO]->(mol_go:GO_Molecular_Function)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, mol_go.GO_ID, mol_go.DEFINITION, mol_go.NAME
    """

    mol_go_mappings = list(graph.run(mol_go_query, entry_id=entry_id))

    go_dict = {}

    if(len(bio_go_mappings) == 0 and len(cell_go_mappings) == 0 and len(mol_go_mappings) == 0):
        return {}, 404

    for mapping in bio_go_mappings:
        (entity_id, auth_asym_id, struct_asym_id, go_id, go_def, go_name) = mapping

        if(go_dict.get(go_id) is None):
            go_dict[go_id] = {
                "category": "Biological_process",
                "definition": go_def,
                "identifier": go_name,
                "name": go_name,
                "mappings": [{
                    "entity_id": int(entity_id),
                    "chain_id": auth_asym_id,
                    "struct_asym_id": struct_asym_id
                }]
            }
        else:
            go_dict[go_id]["mappings"].append({
                "entity_id": int(entity_id),
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id
            })

    for mapping in cell_go_mappings:
        (entity_id, auth_asym_id, struct_asym_id, go_id, go_def, go_name) = mapping

        if(go_dict.get(go_id) is None):
            go_dict[go_id] = {
                "category": "Cellular_component",
                "definition": go_def,
                "identifier": go_name,
                "name": go_name,
                "mappings": [{
                    "entity_id": int(entity_id),
                    "chain_id": auth_asym_id,
                    "struct_asym_id": struct_asym_id
                }]
            }
        else:
            go_dict[go_id]["mappings"].append({
                "entity_id": int(entity_id),
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id
            })

    for mapping in mol_go_mappings:
        (entity_id, auth_asym_id, struct_asym_id, go_id, go_def, go_name) = mapping

        if(go_dict.get(go_id) is None):
            go_dict[go_id] = {
                "category": "Molecular_function",
                "definition": go_def,
                "identifier": go_name,
                "name": go_name,
                "mappings": [{
                    "entity_id": int(entity_id),
                    "chain_id": auth_asym_id,
                    "struct_asym_id": struct_asym_id
                }]
            }
        else:
            go_dict[go_id]["mappings"].append({
                "entity_id": int(entity_id),
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id
            })

    api_result = {}

    for go_id in go_dict.keys():
        api_result[go_id] = go_dict[go_id]

    return api_result, 200


def get_ec(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:CONTAINS_CHAIN]->(chain:Chain), (entity)-[:HAS_EC]->(ec:EC)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, ec.EC_NUM, ec.ACCEPTED_NAME, ec.SYSTEMATIC_NAME, ec.REACTION, ec.SYNONYMS
    """

    dict_ec = {}
    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:

        (entity_id, auth_asym_id, struct_asym_id, ec_num,
         accepted_name, systematic_name, reaction, synonyms) = mapping
    
        # synonyms returns list as string, so converting them as list
        synonyms = synonyms.translate({ord(c): '' for c in "[]' "}).split(',')

        # remove EC from ec_num
        ec_num = ec_num.replace("EC ", "")

        if(dict_ec.get(ec_num) is None):
            dict_ec[ec_num] = {
                "reaction": reaction,
                "mappings": [{
                    "entity_id": int(entity_id),
                    "chain_id": auth_asym_id,
                    "struct_asym_id": struct_asym_id
                }],
                "systematic_name": systematic_name,
                "synonyms": synonyms,
                "identifier": accepted_name,
                "accepted_name": accepted_name
            }
        else:
            dict_ec[ec_num]["mappings"].append({
                "entity_id": int(entity_id),
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id
            })

    api_result = {}

    for key in dict_ec.keys():
        api_result[key] = dict_ec[key]

    return api_result, 200


def get_best_structures(accession, graph):

    query = """
    MATCH(uniprot:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(entity:Entity)-[:HAS_TAXONOMY]->(tax:Taxonomy), (entity)-[:CONTAINS_CHAIN]->(chain:Chain), (entity)<-[:HAS_ENTITY]-(entry:Entry)-[:EXPERIMENT]->(method:Method)
    WITH entry.ID AS entry_id, entity.ID AS entity_id, chain.AUTH_ASYM_ID AS chain_id, collect(DISTINCT pdb_res.ID) AS pdb_residues, collect(DISTINCT unp_res.ID) AS unp_residues, MIN(toInteger(pdb_res.ID)) AS pdb_start, MAX(toInteger(pdb_res.ID)) as pdb_end, MIN(toInteger(unp_res.ID)) as unp_start, MAX(toInteger(unp_res.ID)) as unp_end, tax.TAX_ID AS taxonomy_id, method.METHOD AS experiment, toFloat(entry.RESOLUTION) AS resolution, toFloat(entry.R_FACTOR) as r_factor
    RETURN entry_id, entity_id, chain_id, length(pdb_residues)/(unp_end-unp_start) AS coverage, pdb_start, pdb_end, unp_start, unp_end, taxonomy_id, resolution, r_factor, experiment ORDER BY coverage DESC, resolution
    """

    mappings = list(graph.run(query, accession=accession))

    if(len(mappings) == 0):
        return {}, 404

    api_result = []

    for mapping in mappings:
        (entry_id, entity_id, chain_id, coverage, pdb_start, pdb_end, unp_start, unp_end, taxonomy_id, resolution, r_factor, experiment) = mapping

        api_result.append({
            "end": pdb_end,
            "entity_id": int(entity_id),
            "chain_id": chain_id,
            "pdb_id": entry_id,
            "start": pdb_start,
            "unp_end": unp_end,
            "coverage": None if coverage is None else float("%.3f" % float(coverage)),
            "unp_start": unp_start,
            "resolution": None if resolution is None else float("%.1f" % float(resolution)),
            "experimental_method": experiment,
            "tax_id": int(taxonomy_id)
        })

    return api_result, 200


def get_best_structures_residue_range(accession, unp_start_in, unp_end_in, graph):

    query = """
    MATCH(uniprot:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(entity:Entity)-[:HAS_TAXONOMY]->(tax:Taxonomy), (entity)-[:CONTAINS_CHAIN]->(chain:Chain), (entity)<-[:HAS_ENTITY]-(entry:Entry)-[:EXPERIMENT]->(method:Method) WHERE toInteger(unp_res.ID) IN RANGE($unp_start, $unp_end)
    WITH entry.ID AS entry_id, entity.ID AS entity_id, chain.AUTH_ASYM_ID AS chain_id, collect(DISTINCT pdb_res.ID) AS pdb_residues, collect(DISTINCT unp_res.ID) AS unp_residues, MIN(toInteger(pdb_res.ID)) AS pdb_start, MAX(toInteger(pdb_res.ID)) as pdb_end, MIN(toInteger(unp_res.ID)) as unp_start, MAX(toInteger(unp_res.ID)) as unp_end, tax.TAX_ID AS taxonomy_id, method.METHOD AS experiment, toFloat(entry.RESOLUTION) AS resolution, toFloat(entry.R_FACTOR) as r_factor
    RETURN entry_id, entity_id, chain_id, length(pdb_residues)/(unp_end-unp_start) AS coverage, pdb_start, pdb_end, unp_start, unp_end, taxonomy_id, resolution, r_factor, experiment ORDER BY coverage DESC, resolution
    """

    mappings = list(graph.run(query, accession=accession, unp_start=int(unp_start_in), unp_end=int(unp_end_in)))

    if(len(mappings) == 0):
        return {}, 404

    api_result = []

    for mapping in mappings:
        (entry_id, entity_id, chain_id, coverage, pdb_start, pdb_end, unp_start, unp_end, taxonomy_id, resolution, r_factor, experiment) = mapping

        # skip entry if unp_start and unp_end not lies within range
        if int(unp_start_in) < int(unp_start) or int(unp_end_in) > int(unp_end):
            continue

        api_result.append({
            "end": pdb_end,
            "entity_id": int(entity_id),
            "chain_id": chain_id,
            "pdb_id": entry_id,
            "start": pdb_start,
            "unp_end": unp_end,
            "coverage": None if coverage is None else float("%.3f" % float(coverage)),
            "unp_start": unp_start,
            "resolution": None if resolution is None else float("%.1f" % float(resolution)),
            "experimental_method": experiment,
            "tax_id": int(taxonomy_id)
        })

    return api_result, 200


def get_homologene(entry_id, entity_id, graph):

    query = """
    MATCH (src_entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(src_entity:Entity {ID:$entity_id})-[:HAS_HOMOLOGOUS_UNIPROT]->(uniprot:UniProt)-[:IS_HOMOLOGOUS]->(homologene:Homologene)
    WITH homologene as homologene, uniprot.ACCESSION as accession
    MATCH (homologene)<-[:IS_HOMOLOGOUS]-(uniprot:UniProt {ACCESSION: accession})<-[:HAS_HOMOLOGOUS_UNIPROT]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry), (entity)-[:HAS_TAXONOMY]->(taxonomy:Taxonomy)
    RETURN homologene.HOMOLOGENE_ID, uniprot.ACCESSION, entity.ID, entity.DESCRIPTION, entry.ID, taxonomy.NCBI_SCIENTIFIC
    """

    mappings = list(graph.run(query, entry_id=entry_id, entity_id=entity_id))

    if(len(mappings) == 0):
        return {}, 404

    api_result = {
        "Homologene": {}
    }

    dict_homologene = {}

    for mapping in mappings:
        (homologene_id, accession, entity_id, entity_best_name, entry_id, organism_name) = mapping

        if dict_homologene.get(homologene_id) is None:
            dict_homologene[homologene_id] = {
                "identifier": int(homologene_id),
                "mappings": []
            }

        dict_homologene[homologene_id]["mappings"].append({
            "homologus_pdb_id": entry_id,
            "homologus_pdb_entity_id": int(entity_id),
            "name": entity_best_name,
            "organism_scientific_name": organism_name,
            "accession": accession
        })

    for homologene_id in dict_homologene.keys():

        api_result["Homologene"][homologene_id] = dict_homologene.get(homologene_id)

    return api_result, 200


def get_ensembl(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[ensembl_rel:HAS_ENSEMBL]->(ensembl:Ensembl_Gene),
    (pdb_res)-[:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt),
    (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain)
    RETURN ensembl.GENE_ID, toInteger(entity.ID), chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, toInteger(ensembl_rel.ORDINAL), ensembl_rel.COVERAGE, ensembl_rel.EXON_ID, 
    ensembl_rel.TRANSCRIPT_ID, ensembl_rel.TRANSLATION_ID, toInteger(ensembl_rel.GENOME_START), toInteger(ensembl_rel.GENOME_END), unp.ACCESSION, MIN(toInteger(pdb_res.ID)) AS pdb_res_start, 
    MAX(toInteger(pdb_res.ID)) AS pdb_res_end, MIN(toInteger(unp_res.ID)) AS unp_res_start, MAX(toInteger(unp_res.ID)) AS unp_res_end, MIN(toInteger(chain_rel.AUTH_SEQ_ID)) AS auth_seq_start, 
    MAX(toInteger(chain_rel.AUTH_SEQ_ID)) AS auth_seq_end
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    api_result = {
        "Ensembl": {}
    }

    dict_ensembl = {}

    for mapping in mappings:
        (gene_id, entity_id, auth_asym_id, struct_asym_id, ordinal, coverage, exon_id, transcript_id, translation_id, genome_start, genome_end, accession, pdb_start,
        pdb_end, unp_start, unp_end, auth_seq_start, auth_seq_end) = mapping

        if dict_ensembl.get(gene_id) is None:
            dict_ensembl[gene_id] = {
                "identifier": gene_id,
                "mappings": []
            }

        dict_ensembl[gene_id]["mappings"].append({
            "ordinal": ordinal,
            "entity_id": entity_id,
            "end": {
                "author_residue_number": auth_seq_end,
                "author_insertion_code": "",
                "residue_number": pdb_end
            },
            "accession": accession,
            "genome_start": genome_start,
            "chain_id": auth_asym_id,
            "start": {
                "author_residue_number": auth_seq_start,
                "author_insertion_code": "",
                "residue_number": pdb_start
            },
            "unp_end": unp_end,
            "transcript_id": transcript_id,
            "exon_id": exon_id,
            "coverage": float("%.3f" % float(coverage)),
            "unp_start": unp_start,
            "translation_id": translation_id,
            "genome_end": genome_end,
            "struct_asym_id": struct_asym_id
        })

    for gene_id in dict_ensembl.keys():

        api_result["Ensembl"][gene_id] = dict_ensembl.get(gene_id)

    return api_result, 200

def get_uniprot_to_pfam(accession, graph):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)-[:IS_IN_PFAM]->(pfam:Pfam)
    RETURN pfam.PFAM_ACCESSION, pfam.NAME, pfam.DESCRIPTION, toInteger(unp_res.ID) ORDER BY toInteger(unp_res.ID)
    """

    dict_pfam = {}
    api_result = {}

    mappings = list(graph.run(query, accession=accession))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:
        (pfam_accession, pfam_name, desc, res_id) = mapping

        if(dict_pfam.get(pfam_accession) is None):

            api_result[pfam_accession] = {
                "description": desc,
                "identifier": desc,
                "name": pfam_name,
                "mappings": []
            }
            dict_pfam[pfam_accession] = [res_id]

        else:
            dict_pfam[pfam_accession].append(res_id)

    for pfam_accession in dict_pfam.keys():

        for group in mit.consecutive_groups(dict_pfam[pfam_accession]):
            group = list(group)
            api_result[pfam_accession]["mappings"].append({
                "unp_start": group[0],
                "unp_end": group[-1]
            })

    return api_result, 200


def get_isoforms(entry_id, graph, mapping_type):

    best_query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[unp_rel:HAS_UNIPROT {BEST_MAPPING:'1'}]->(unp:UniProt)-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(entity)
    MATCH (entity)-[:CONTAINS_CHAIN]->(chain:Chain)<-[chain_rel:IS_IN_CHAIN]-(pdb_res)
    RETURN unp.ACCESSION, unp.NAME, entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, chain_rel.AUTH_SEQ_ID, chain_rel.OBSERVED_RATIO, pdb_res.ID, unp_res.ID, unp_rel.IDENTITY ORDER BY toInteger(pdb_res.ID)
    """

    all_query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[unp_rel:HAS_UNIPROT]->(unp:UniProt)-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(entity)
    MATCH (entity)-[:CONTAINS_CHAIN]->(chain:Chain)<-[chain_rel:IS_IN_CHAIN]-(pdb_res)
    RETURN unp.ACCESSION, unp.NAME, entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, chain_rel.AUTH_SEQ_ID, chain_rel.OBSERVED_RATIO, pdb_res.ID, unp_res.ID, unp_rel.IDENTITY ORDER BY toInteger(pdb_res.ID)
    """

    if(mapping_type == 'B'):
        mappings = list(graph.run(best_query, entry_id=entry_id))
    elif(mapping_type == 'A'):
        mappings = list(graph.run(all_query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    dict_unp_master = {}
    dict_mappings = {}

    api_result = {
        "UniProt": {}
    }

    for mapping in mappings:
        (accession, name, entity_id, auth_asym_id, struct_asym_id, auth_seq_id, observed_ratio, pdb_res_id, unp_res_id, identity) = mapping

        # AUTH_SEQ_ID is populated from RESIDUE not SIFTS_XREF_RESIDUE, AUTH_SEQ_ID is null if OBSERVED_RATIO is 0
        if(observed_ratio == '0.0'):
            auth_seq_id = None

        if(dict_unp_master.get(accession) is None):
            dict_unp_master[accession] = True
            api_result["UniProt"][accession] = {
                "identifier": name,
                "name": name,
                "mappings": []
            }
        
        key = (accession, entity_id, auth_asym_id, struct_asym_id, identity)

        if(dict_mappings.get(key) is None):
            dict_mappings[key] = [(auth_seq_id, pdb_res_id, unp_res_id)]
        else:
            dict_mappings[key].append((auth_seq_id, pdb_res_id, unp_res_id))

    for key in dict_mappings.keys():
        (accession, entity_id, auth_asym_id, struct_asym_id, identity) = key

        start = dict_mappings[key][0]
        end = dict_mappings[key][-1]

        api_result["UniProt"][accession]["mappings"].append({
            "entity_id": int(entity_id),
            "chain_id": auth_asym_id,
            "struct_asym_id": struct_asym_id,
            "pdb_start": int(start[1]),
            "pdb_end": int(end[1]),
            "start": {
                "author_residue_number": None if start[0] is None else int(start[0]),
                "author_insertion_code": "",
                "residue_number": int(start[1])
            },
            "end": {
                "author_residue_number": None if end[0] is None else int(end[0]),
                "author_insertion_code": "",
                "residue_number": int(end[1])
            },
            "identity": float("%.3f" % float(identity)),
            "unp_start": int(start[2]),
            "unp_end": int(end[2])
        })

    return api_result, 200


def get_uniref90(entry_id, graph, type_of_query):

    homologene_unf_query = """
    MATCH(entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)-[:MAP_TO_UNIREF_RESIDUE]->(unf_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unf:UniProt)-[:HAS_TAXONOMY]->(tax:Taxonomy), (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain), (unp_res)<-[:HAS_UNP_RESIDUE]-(unp:UniProt)-[:IS_HOMOLOGOUS]->(homologene:Homologene)<-[:IS_HOMOLOGOUS]-(unf)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, unf.ACCESSION, unf.NAME, tax.TAX_ID, tax.NCBI_SCIENTIFIC, tax.NCBI_COMMON, toInteger(pdb_res.ID), toInteger(chain_rel.AUTH_SEQ_ID), pdb_res.CHEM_COMP_ID, toInteger(unf_res.ID), unf_res.ONE_LETTER_CODE ORDER BY toInteger(pdb_res.ID)
    """

    all_unf_query = """
    MATCH(entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)-[:MAP_TO_UNIREF_RESIDUE]->(unf_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unf:UniProt)-[:HAS_TAXONOMY]->(tax:Taxonomy), (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, unf.ACCESSION, unf.NAME, tax.TAX_ID, tax.NCBI_SCIENTIFIC, tax.NCBI_COMMON, toInteger(pdb_res.ID), toInteger(chain_rel.AUTH_SEQ_ID), pdb_res.CHEM_COMP_ID, toInteger(unf_res.ID), unf_res.ONE_LETTER_CODE ORDER BY toInteger(pdb_res.ID)
    """

    all_unp_query = """
    MATCH(entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt)-[:HAS_TAXONOMY]->(tax:Taxonomy), (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, unp.ACCESSION, unp.NAME, tax.TAX_ID, tax.NCBI_SCIENTIFIC, tax.NCBI_COMMON, toInteger(pdb_res.ID), toInteger(chain_rel.AUTH_SEQ_ID), pdb_res.CHEM_COMP_ID, toInteger(unp_res.ID), unp_res.ONE_LETTER_CODE ORDER BY toInteger(pdb_res.ID)
    """

    homologene_unp_query = """
    MATCH(entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp:UniProt)-[:HAS_TAXONOMY]->(tax:Taxonomy), (pdb_res)-[chain_rel:IS_IN_CHAIN]->(chain:Chain), (unp)-[:IS_HOMOLOGOUS]->(homologene:Homologene)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, unp.ACCESSION, unp.NAME, tax.TAX_ID, tax.NCBI_SCIENTIFIC, tax.NCBI_COMMON, toInteger(pdb_res.ID), toInteger(chain_rel.AUTH_SEQ_ID), pdb_res.CHEM_COMP_ID, toInteger(unp_res.ID), unp_res.ONE_LETTER_CODE ORDER BY toInteger(pdb_res.ID)
    """

    if type_of_query == 'A':
        unf_mappings = list(graph.run(all_unf_query, entry_id=entry_id))
        unp_mappings = list(graph.run(all_unp_query, entry_id=entry_id))
    else:
        unf_mappings = list(graph.run(homologene_unf_query, entry_id=entry_id))
        unp_mappings = list(graph.run(homologene_unp_query, entry_id=entry_id))

    api_result = {
        "UniProt": {}
    }

    if(len(unf_mappings) == 0 and len(unp_mappings) == 0):
        return api_result, 404

    api_result = {
        "UniProt": {}
    }

    dict_details = {}

    mappings = unp_mappings + unf_mappings
   
    for mapping in mappings:
        
        (entity_id, auth_asym_id, struct_asym_id, unf_accession, unf_name, tax_id, scientific_name, common_name, pdb_res_id, auth_seq_id, chem_comp_id, unf_res_id, unf_one_letter_code) = mapping

        if api_result["UniProt"].get(unf_accession) is None:
            api_result["UniProt"][unf_accession] = {
                "identifier": unf_name,
                "name": unf_name,
                "organism_scientific_name": "{} {}".format(scientific_name, "" if common_name is None else "({})".format(common_name)),
                "tax_id": int(tax_id),
                "mappings": []
            }

        key = (unf_accession, entity_id, auth_asym_id, struct_asym_id)
        if dict_details.get((key)) is None:
            dict_details[key] = [(pdb_res_id, auth_seq_id, chem_comp_id, unf_res_id, unf_one_letter_code)]
        else:
            dict_details[key].append((pdb_res_id, auth_seq_id, chem_comp_id, unf_res_id, unf_one_letter_code))


    for key in dict_details.keys():
        (unf_accession, entity_id, auth_asym_id, struct_asym_id) = key
        details = dict_details.get(key)

        pdb_start = details[0][0]
        pdb_end = details[-1][0]
        auth_start = details[0][1]
        auth_end = details[-1][1]
        unf_start = details[0][3]
        unf_end = details[-1][3]
        count_identities = 0

        for detail in details:
            (pdb_res_id, auth_seq_id, chem_comp_id, unf_res_id, unf_one_letter_code) = detail

            (pdb_one_letter_code, chem_desc) = amino_acid_codes[chem_comp_id]

            if pdb_one_letter_code == unf_one_letter_code:
                count_identities += 1

        identity = count_identities/len(details)

        api_result["UniProt"][unf_accession]["mappings"].append({
            "entity_id": int(entity_id),
            "end": {
                "author_residue_number": None if auth_end is None else int(auth_end),
                "author_insertion_code": "",
                "residue_number": int(pdb_end)
            },
            "chain_id": auth_asym_id,
            "pdb_start": int(pdb_start),
            "start": {
                "author_residue_number": None if auth_start is None else int(auth_start),
                "author_insertion_code": "",
                "residue_number": int(pdb_start)
            },
            "unp_end": int(unf_end),
            "pdb_end": int(pdb_end),
            "struct_asym_id": struct_asym_id,
            "unp_start": int(unf_start),
            "identity": float("%.3f" % float(identity))
        })


    return api_result, 200
