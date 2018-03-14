from flask import Flask, jsonify
from .model import graph

app = Flask(__name__)

tasks = [
    {
        'id': 1,
        'title': u'Buy groceries',
        'description': u'Milk, Cheese, Pizza, Fruit, Tylenol', 
        'done': False
    },
    {
        'id': 2,
        'title': u'Learn Python',
        'description': u'Need to find a good Python tutorial on the web', 
        'done': False
    }
]

@app.route('/')
def hello():
    print(graph)
    return jsonify({'tasks': tasks})

@app.route('/api/mappings/uniprot/<string:entry_id>')
def get_uniprot(entry_id):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_UNIPROT]->(unp:UniProt)-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)
    MATCH (entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[r:MAP_TO_UNIPROT_RESIDUE]->(unp_res)
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



@app.route('/api/mappings/unipdb/<string:uniprot_accession>')
def get_unipdb(uniprot_accession):

    print(uniprot_accession)

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})<-[:HAS_UNIPROT]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry),
    (entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[relation:MAP_TO_UNIPROT_RESIDUE]->(unp_res:UNP_Residue)<-[:HAS_UNP_RESIDUE]-(unp)
    WITH entry.ID AS entryId, entity.ID AS entityId, relation.CHAINS AS chains, toInteger(pdb_res.ID) AS pdbRes, toInteger(unp_res.ID) AS unpRes
    RETURN entryId, entityId, chains, min(pdbRes) AS pdbStart, max(pdbRes) AS pdbEnd, min(unpRes) AS unpStart, max(unpRes) AS unpEnd
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

    for mapping in mappings:
        
        (entry, entity, chains, pdb_start, pdb_end, unp_start, unp_end) = mapping
        final_result[uniprot_accession]['mappings'].append({
            'entry_id': entry,
            'entity_id': entity,
            'chains': chains,
            'pdb_start': pdb_start,
            'pdb_end': pdb_end,
            'unp_start': unp_start,
            'unp_end': unp_end
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