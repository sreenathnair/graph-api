

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


def get_validation_rama_sidechain_listing(entry_id, graph):
    
    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_CHAIN {OBSERVED:'Y'}]->(chain:Chain)
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, rel.MODEL, rel.PSI, rel.CIS_PEPTIDE, pdb_res.ID, rel.AUTH_SEQ_ID, rel.RAMA, rel.PHI, pdb_res.CHEM_COMP_ID, rel.ROTA
    """

    dict_models = {}
    dict_chains = {}
    list_entities = []
    dict_residues = {}
    mappings = list(graph.run(query, entry_id=entry_id))

    for mapping in mappings:
        (entity_id, auth_asym_id, struct_asym_id, model, psi, cis_peptide, res_id, auth_seq_id, rama, phi, chem_comp_id, rota) = mapping

        key = (entity_id, auth_asym_id, struct_asym_id, model)

        if(entity_id not in list_entities):
            list_entities.append(entity_id)

        if(dict_chains.get(entity_id) is None):
            dict_chains[entity_id] = [(auth_asym_id, struct_asym_id)]
        elif((auth_asym_id, struct_asym_id) not in dict_chains[entity_id]):
            dict_chains[entity_id].append((auth_asym_id, struct_asym_id))

        model_key = (entity_id, auth_asym_id)

        if(dict_models.get(model_key) is None):
            dict_models[model_key] = [model]
        elif(model not in dict_models[model_key]):
            dict_models[model_key].append(model)

        if(dict_residues.get(key) is None):
            dict_residues[key] = [(psi, cis_peptide, res_id, auth_seq_id, rama, phi, chem_comp_id, rota)]
        else:
            dict_residues[key].append((psi, cis_peptide, res_id, auth_seq_id, rama, phi, chem_comp_id, rota))

    
    api_result = {
        "molecules": []
    }

    for entity in list_entities:

        entity_element = {
            "entity_id": int(entity),
            "chains": []
        }

        for (auth_asym_id, struct_asym_id) in dict_chains[entity]:
            
            chain_element = {
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id,
                "models": []
            }

            for model in dict_models[(entity, auth_asym_id)]:

                model_element = {
                    "model_id": int(model),
                    "residues": []
                }

                for (psi, cis_peptide, res_id, auth_seq_id, rama, phi, chem_comp_id, rota) in dict_residues[(entity, auth_asym_id, struct_asym_id, model)]:
                    
                    if(psi == None):
                        psi = ""

                    if(phi == None):
                        phi = ""

                    residue_element = {
                        "psi": None if psi == '' else float("%.1f" % float(psi)),
                        "cis_peptide": cis_peptide,
                        "residue_number": int(res_id),
                        "author_residue_number": int(auth_seq_id),
                        "rama": rama,
                        "phi": None if phi == '' else float("%.1f" % float(phi)),
                        "author_insertion_code": "",
                        "residue_name": chem_comp_id,
                        "rota": rota
                    }

                    model_element["residues"].append(residue_element)

                chain_element["models"].append(model_element)

            entity_element["chains"].append(chain_element)

        api_result["molecules"].append(entity_element)
    
    return api_result