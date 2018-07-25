#!/usr/bin/env python

"""
uniprot.py: This have list of function calls which retrieves data for a UniProt accession/residue
"""

from .amino_acid_codes import amino_acid_codes
from .residue import get_mappings_for_residue_pfam, get_mappings_for_residue_scop, get_mappings_for_residue_cath, get_mappings_for_residue_interpro


__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"


def get_unipdb(uniprot_accession, graph):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)<-[relation:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)
    <-[:HAS_PDB_RESIDUE]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry)
    WITH entry.ID AS entryId, entity.ID AS entityId, relation.CHAINS AS chains, toInteger(pdb_res.ID) AS pdbRes, toInteger(unp_res.ID) AS unpRes, pdb_res.CHEM_COMP_ID AS pdbResCode, unp_res.ONE_LETTER_CODE AS unpResCode 
    RETURN entryId, entityId, chains, pdbRes, unpRes, pdbResCode, unpResCode ORDER BY entryId, pdbRes
    """

    mappings = list(graph.run(query, parameters={
        'accession': str(uniprot_accession)
    }))
    
    if(len(mappings) == 0):
        return {}, 404

    final_result = {
        'Pfam': get_mappings_for_unp_residue_pfam(uniprot_accession, graph),
        'mappings': []
    }

    entry_entity_dict = {}

    for mapping in mappings:

        (entry, entity, chains, pdb_res, unp_res,
         pdb_res_code, unp_res_code) = mapping
    
        pdb_res_desc = None

        # modified/mutated residues will be different, handle them
        if(amino_acid_codes.get(pdb_res_code) is None):
            pdb_res_one_letter = pdb_res_code
        else:
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
            (chains, pdb_res, unp_res, pdb_res_code,
             unp_res_code) = entry_entity_dict[key][incr]

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
                # pdb_res and unp_end will be ending residue numbers
                final_map[key].append(
                    (chains, pdb_seq, unp_seq, pdb_start, prev_pdb_res, unp_start, prev_unp_res))
                pdb_seq = pdb_res_code
                unp_seq = unp_res_code
                start = True

            # the very last residue
            if incr == len(entry_entity_dict[key]) - 1:
                # pdb_res and unp_end will be ending residue numbers
                final_map[key].append(
                    (chains, pdb_seq, unp_seq, pdb_start, pdb_res, unp_start, unp_res))

            prev_pdb_res = pdb_res
            prev_unp_res = unp_res
            incr += 1

    for key in final_map.keys():
        entry_id, entity_id = key

        mappings = []

        for mapping in final_map[key]:
            (chains, pdb_seq, unp_seq, pdb_start,
             pdb_end, unp_start, unp_end) = mapping
            mappings.append({
                'chains': chains,
                'pdb_sequence': pdb_seq,
                'unp_sequence': unp_seq,
                'pdb_start': pdb_start,
                'pdb_end': pdb_end,
                'unp_start': unp_start,
                'unp_end': unp_end
            })

        final_result['mappings'].append({
            'entry_id': entry_id,
            'entity_id': int(entity_id),
            'segments': mappings
        })

    return final_result, 200


def get_mappings_for_unp_residue_pfam(uniprot_accession, graph):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue)-[:IS_IN_PFAM]->(pfam:Pfam)
    WITH pfam.PFAM_ACCESSION AS pfamAccession, pfam.NAME AS pfamName, pfam.DESCRIPTION AS pfamDesc, toInteger(unp_res.ID) AS unpRes
    RETURN pfamAccession, pfamName, pfamDesc, min(unpRes) AS unpStart, max(unpRes) AS unpEnd
    """

    result = list(graph.run(query, parameters={
        'accession': str(uniprot_accession)
    }))

    final_result = {}

    for pfam in result:

        (accession, name, desc, unp_start, unp_end) = pfam
        final_result = {
            "identifier": name,
            "name": name,
            "description": desc,
            "unp_start": unp_start,
            "unp_end": unp_end
        }

    return final_result


def get_unipdb_residue(uniprot_accession, unp_res, graph):

    query = """
    MATCH (unp:UniProt {ACCESSION:$accession})-[:HAS_UNP_RESIDUE]->(unp_res:UNP_Residue {ID:$unpResidue})<-[:MAP_TO_UNIPROT_RESIDUE]-(pdb_res:PDB_Residue)<-[:HAS_PDB_RESIDUE]-(entity:Entity)<-[:HAS_ENTITY]-(entry:Entry)
    RETURN entry.ID AS entryId, entity.ID AS entityId, pdb_res.ID AS pdbRes
    """

    mappings = list(graph.run(query, parameters={
        'accession': str(uniprot_accession), 'unpResidue': str(unp_res)
    }))

    if(len(mappings) == 0):
        return {}, 404

    final_result = {}

    pfam_dict = {}
    interpro_dict = {}
    cath_dict = {}
    scop_dict = {}

    for mapping in mappings:

        (entry_id, entity_id, pdb_res) = mapping
        temp_map, resp_status = get_mappings_for_residue_pfam(
            entry_id, entity_id, pdb_res, graph)

        for pfam in temp_map.keys():
            if pfam_dict.get(pfam) is None:
                pfam_dict[pfam] = temp_map[pfam]

        temp_map, resp_status = get_mappings_for_residue_interpro(
            entry_id, entity_id, pdb_res, graph)

        for interpro in temp_map.keys():
            if interpro_dict.get(interpro) is None:
                interpro_dict[interpro] = temp_map[interpro]

        temp_map, resp_status = get_mappings_for_residue_cath(
            entry_id, entity_id, pdb_res, graph)

        for cath in temp_map.keys():
            if cath_dict.get(cath) is None:
                cath_dict[cath] = temp_map[cath]

        temp_map, resp_status = get_mappings_for_residue_scop(
            entry_id, entity_id, pdb_res, graph)

        for scop in temp_map.keys():
            if scop_dict.get(scop) is None:
                scop_dict[scop] = temp_map[scop]

    final_result['Pfam'] = pfam_dict
    final_result['InterPro'] = interpro_dict
    final_result['CATH'] = cath_dict
    final_result['SCOP'] = scop_dict

    return final_result, 200