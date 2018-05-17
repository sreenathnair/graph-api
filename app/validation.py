

def get_validation_protein_ramachandran_sidechain_outliers(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_CHAIN {OBSERVED:'Y'}]->(chain:Chain)
    WHERE rel.RAMA='OUTLIER' OR rel.ROTA='OUTLIER'
    RETURN entity.ID, rel.RAMA, rel.ROTA, rel.MODEL, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, pdb_res.ID, rel.AUTH_SEQ_ID
    """

    list_rama = []
    list_rota = []
    mappings = list(graph.run(query, entry_id=entry_id))

    for mapping in mappings:

        (entity_id, rama, rota, model, auth_asym_id, struct_asym_id, res_id, auth_seq_id) = mapping

        if(rama == 'OUTLIER'):
            list_rama.append({
                "model_id": int(model),
                "entity_id": int(entity_id),
                "residue_number": int(res_id),
                "author_residue_number": int(auth_seq_id),
                "chain_id": auth_asym_id,
                "author_insertion_code": "",
                "alt_code": "",
                "struct_asym_id": struct_asym_id
            })

        if(rota == 'OUTLIER'):
            list_rota.append({
                "model_id": int(model),
                "entity_id": int(entity_id),
                "residue_number": int(res_id),
                "author_residue_number": int(auth_seq_id),
                "chain_id": auth_asym_id,
                "author_insertion_code": "",
                "alt_code": "",
                "struct_asym_id": struct_asym_id
            })

    return {
        "ramachandran_outliers": list_rama,
        "sidechain_outliers": list_rota
    }