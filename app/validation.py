from decimal import Decimal
import copy

RSRZ_OUTLIER_CUTOFF = 2
RNA_suite_not_nonRotameric = ["NotAvailable","Rotameric",None]

def get_validation_protein_ramachandran_sidechain_outliers(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_CHAIN {OBSERVED:'Y'}]->(chain:Chain)
    WHERE rel.RAMA='OUTLIER' OR rel.ROTA='OUTLIER'
    RETURN entity.ID, rel.RAMA, rel.ROTA, rel.MODEL, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, pdb_res.ID, rel.AUTH_SEQ_ID
    """

    list_rama = []
    list_rota = []
    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

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
    }, 200


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

    if(len(mappings) == 0):
        return {}, 404

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
    
    return api_result, 200


def get_validation_rna_pucker_suite_outliers(entry_id, graph):

    # sample entry : 3j8g
    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[rel:IS_IN_CHAIN]->(chain:Chain)
    WHERE rel.RNA_SUITE IN ['NonRotameric','Triaged/NotBinned'] OR rel.RNA_PUCKER='outlier'
    RETURN toInteger(entity.ID), chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, toInteger(pdb_res.ID), toInteger(rel.MODEL), rel.RNA_SUITE, rel.RNA_PUCKER, toInteger(rel.AUTH_SEQ_ID)
    """

    list_rna_pucker = []
    list_rna_suite = []
    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:
        (entity_id, auth_asym_id, struct_asym_id, res_id, model, rna_suite, rna_pucker, auth_seq_id) = mapping

        if(rna_suite in ['NonRotameric','Triaged/NotBinned']):
            list_rna_suite.append({
                "model_id": model,
                "entity_id": entity_id,
                "residue_number": res_id,
                "author_residue_number": auth_seq_id,
                "chain_id": auth_asym_id,
                "author_insertion_code": "",
                "alt_code": "",
                "struct_asym_id": struct_asym_id
            })

        if(rna_pucker == 'outlier'):
            list_rna_pucker.append({
                "model_id": model,
                "entity_id": entity_id,
                "residue_number": res_id,
                "author_residue_number": auth_seq_id,
                "chain_id": auth_asym_id,
                "author_insertion_code": "",
                "alt_code": "",
                "struct_asym_id": struct_asym_id
            })

    return {
        "pucker_outliers": list_rna_pucker,
        "suite_outliers": list_rna_suite
    }, 200

def get_validation_global_percentiles(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})
    RETURN entry.ABS_PERCENTILE_PERCENT_RSRZ,entry.REL_PERCENTILE_PERCENT_RSRZ,entry.PERCENT_RSRZ_OUTLIERS, entry.REL_PERCENTILE_PERCENT_ROTA,entry.ABS_PERCENTILE_PERCENT_ROTA,
    entry.PERCENT_ROTA_OUTLIERS,entry.REL_PERCENTILE_PERCENT_RAMA,entry.ABS_PERCENTILE_PERCENT_RAMA,entry.PERCENT_RAMA_OUTLIERS,entry.REL_PERCENTILE_CLASHSCORE,
    entry.ABS_PERCENTILE_CLASHSCORE,entry.CLASHSCORE,entry.REL_PERCENTILE_DCC_RFREE,entry.ABS_PERCENTILE_DCC_RFREE,entry.DCC_RFREE,entry.ABS_PERCENTILE_RNA_SUITENESS,
    entry.REL_PERCENTILE_RNA_SUITENESS,entry.RNA_SUITENESS
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    api_result = {}

    for mapping in mappings:
        (abs_rsrz,rel_rsrz,raw_rsrz,rel_rota,abs_rota,raw_rota,rel_rama,abs_rama,raw_rama,rel_clash,abs_clash,raw_clash,rel_dcc,abs_dcc,raw_dcc,abs_rna,rel_rna,raw_rna) = mapping

        if(rel_rsrz is not None):
            api_result["percent-RSRZ-outliers"] = {
                "relative": float("%.2f" % float(rel_rsrz)),
                "rawvalue": float("%.2f" % float(raw_rsrz)),
                "absolute": float("%.2f" % float(abs_rsrz))
            }
        if(rel_clash is not None):
            api_result["clashscore"] = {
                "relative": float("%.2f" % float(rel_clash)),
                "rawvalue": float("%.2f" % float(raw_clash)),
                "absolute": float("%.2f" % float(abs_clash))
            }
        if(rel_rota is not None):
            api_result["percent-rota-outliers"] = {
                "relative": float("%.2f" % float(rel_rota)),
                "rawvalue": float("%.2f" % float(raw_rota)),
                "absolute": float("%.2f" % float(abs_rota))
            }
        if(rel_rama is not None):
            api_result["percent-rama-outliers"] = {
                "relative": float("%.2f" % float(rel_rama)),
                "rawvalue": float("%.2f" % float(raw_rama)),
                "absolute": float("%.2f" % float(abs_rama))
            }
        if(rel_dcc is not None):
            api_result["DCC_Rfree"] = {
                "relative": float("%.2f" % float(rel_dcc)),
                "rawvalue": float("%.2f" % float(raw_dcc)),
                "absolute": float("%.2f" % float(abs_dcc))
            }
        if(rel_rna is not None):
            api_result["RNAsuiteness"] = {
                "relative": float("%.2f" % float(rel_rna)),
                "rawvalue": float("%.2f" % float(raw_rna)),
                "absolute": float("%.2f" % float(abs_rna))
            }

        # only 1 record is returned
        return api_result, 200


def get_validation_summary_quality_scores(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:EXPERIMENT]->(method:Method)
    RETURN entry.DATA_QUALITY,entry.OVERALL_QUALITY,entry.GEOMETRY_QUALITY,method.METHOD,method.METHOD_CLASS
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:

        (data_quality, overall_quality, geo_quality, method, method_class) = mapping
        
        # only 1 record is returned
        return {
            "overall_quality": None if(overall_quality is None) else float("%.2f" % float(overall_quality)),
            "geometry_quality": None if(geo_quality is None) else float("%.2f" % float(geo_quality)),
            "experiment_data_available": True if (method == 'X-ray diffraction' and method_class == 'x-ray') else "unknown",
            "data_quality": None if(data_quality is None) else float("%.2f" % float(data_quality))
        }, 200

def get_validation_key_validation_stats(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)
    OPTIONAL MATCH (pdb_res)-[:HAS_VALIDATION_DATA]->(val_bond:Val_Bond_Outlier)
    OPTIONAL MATCH (pdb_res)-[:HAS_VALIDATION_DATA]->(val_angle:Val_Angle_Outlier)
    RETURN entry.ANGLES_RMSZ,entry.NUM_ANGLES_RMSZ,entry.BONDS_RMSZ,entry.NUM_BONDS_RMSZ, count(val_bond) as num_val_bonds,count(val_angle) as num_val_angle
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    angle_rmsz, num_angle_rmsz, bonds_rmsz, num_bonds_rmsz, num_bond_outliers, num_angle_outliers = mappings[0]
    percent_bond_outliers = float((Decimal(num_bond_outliers) * Decimal(100.0) / Decimal(num_bonds_rmsz)))
    percent_angle_outliers = float((Decimal(num_angle_outliers) * Decimal(100.0) / Decimal(num_angle_rmsz)))

    api_result = {
        "bonds": {
            "rmsz":  float("%.2f" % float(bonds_rmsz)),
            "num_checked": int(num_bonds_rmsz),
            "percent_outliers": float("%.2f" % float(percent_bond_outliers)),
            "num_outliers": int(num_bond_outliers)
        },
        "angles": {
            "rmsz": float("%.2f" % float(angle_rmsz)),
            "num_checked": int(num_angle_rmsz),
            "percent_outliers": float("%.2f" % float(percent_angle_outliers)),
            "num_outliers": int(num_angle_outliers)
        }
    }

    query = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[chain_rel:IS_IN_CHAIN {OBSERVED:'Y'}]->(chain:Chain)
    RETURN chain_rel.RAMA AS rama,chain_rel.ROTA AS rota,chain_rel.RSRZ AS rsrz,chain_rel.RNA_SUITE AS rna_suite,chain_rel.RNA_PUCKER AS rna_pucker
    """

    mappings = list(graph.run(query, entry_id=entry_id))
    keys = ["rama", "rota", "rsrz", "rna_suite", "rna_pucker"]

    for key in keys:
        if(api_result.get(key) is None):
            api_result[key] = {
                "num_checked": 0,
                "num_outliers": 0
            }

    for mapping in mappings:

        for key in keys:
            value = mapping.get(key)
            
            if(value is not None):
                api_result[key]["num_checked"] += 1
                if(key in ["rama", "rota"] and value == "OUTLIER"):
                    api_result[key]["num_outliers"] += 1
                elif(key == "rsrz" and float(value) > RSRZ_OUTLIER_CUTOFF):
                    api_result[key]["num_outliers"] += 1
                elif(key == "rna_suite" and value not in RNA_suite_not_nonRotameric):
                    api_result[key]["num_outliers"] += 1
                elif(key == "rna_pucker" and value == "outlier"):
                    api_result[key]["num_outliers"] += 1
                    
    # correcting rna_pucker
    api_result["rna_pucker"]["num_checked"] = api_result["rna_suite"]["num_checked"]

    for key in keys:
        num_checked = Decimal(api_result[key]["num_checked"])
        num_outliers = Decimal(api_result[key]["num_outliers"])
        percent_outliers = None
        
        if(num_checked != 0):
            percent_outliers = (num_outliers * Decimal(100.0) / num_checked)

        api_result[key]["percent_outliers"] = None if percent_outliers is None else float("%.2f" % float(percent_outliers))

    # renaming keys
    api_result["protein_ramachandran"] = api_result["rama"]
    del api_result["rama"]

    api_result["RSRZ"] = api_result["rsrz"]
    del api_result["rsrz"]

    api_result["protein_sidechains"] = api_result["rota"]
    del api_result["rota"]


    return api_result, 200


def get_validation_xray_refine_data_stats(entry_id, graph):

    query = """
    MATCH (entry:Entry {ID:$entry_id})
    RETURN entry.DATA_COMPLETENESS,entry.NUM_FREE_REFLECTIONS,toInteger(entry.NUM_MILLER_INDICES),entry.TRANS_NCS,entry.DCC_RFREE, entry.DCC_REFINEMENT_PROGRAM,toInteger(entry.CENTRIC_OUTLIERS),
    entry.TWIN_L, entry.EDS_R,entry.PERCENT_FREE_REFLECTIONS,toInteger(entry.ACENTRIC_OUTLIERS),entry.DCC_R,entry.EDS_RES_LOW,entry.WILSON_B_ESTIMATE,entry.BULK_SOLVENT_K,entry.EDS_RES,
    entry.FO_FC_CORRELATION, entry.I_OVER_SIGMA,entry.TWIN_L2,entry.BULK_SOLVENT_B
    """

    mappings = list(graph.run(query, entry_id=entry_id))

    if(len(mappings) == 0):
        return {}, 404

    (data_completeness,num_free_reflections,num_miller_indices,trans_ncs,dcc_rfree,dcc_refinement_program,centric_outliers,twin_l,eds_r,percent_free_reflections,
    acentric_outliers,dcc_r,eds_res_low,wilson_b_estimate,bulk_solvent_k,eds_res,fo_fc_correlation,i_over_sigma,twin_l2,bulk_solvent_b) = mappings[0]

    return {
        "DataCompleteness": {
            "source": "EDS",
            "value": float("%.2f" % float(data_completeness))
        },
        "num-free-reflections": {
            "source": "EDS",
            "value": int(num_free_reflections)
        },
        "numMillerIndices": {
            "source": "Xtriage(Phenix)",
            "value": num_miller_indices
        },
        "TransNCS": {
            "source": "Xtriage(Phenix)",
            "value": trans_ncs
        },
        "DCC_refinement_program": {
            "source": "DCC",
            "value": dcc_refinement_program
        },
        "centric_outliers": {
            "source": "Xtriage(Phenix)",
            "value": centric_outliers
        },
        "TwinL": {
            "source": "Xtriage(Phenix)",
            "value": float("%.2f" % float(twin_l))
        },
        "EDS_R": {
            "source": "EDS",
            "value": float("%.2f" % float(eds_r))
        },
        "DCC_Rfree": {
            "source": "DCC",
            "value": float("%.2f" % float(dcc_rfree))
        },
        "percent-free-reflections": {
            "source": "EDS",
            "value": float("%.2f" % float(percent_free_reflections))
        },
        "acentric_outliers": {
            "source": "Xtriage(Phenix)",
            "value": acentric_outliers
        },
        "DCC_R": {
            "source": "DCC",
            "value": float("%.2f" % float(dcc_r))
        },
        "EDS_resolution_low": {
            "source": "EDS",
            "value": float("%.2f" % float(eds_res_low))
        },
        "WilsonBestimate": {
            "source": "Xtriage(Phenix)",
            "value": float("%.3f" % float(wilson_b_estimate))
        },
        "bulk_solvent_k": {
            "source": "EDS",
            "value": float("%.3f" % float(bulk_solvent_k))
        },
        "EDS_resolution": {
            "source": "EDS",
            "value": float("%.2f" % float(eds_res))
        },
        "Fo_Fc_correlation": {
            "source": "EDS",
            "value": float("%.3f" % float(fo_fc_correlation))
        },
        "IoverSigma": {
            "source": "Xtriage(Phenix)",
            "value": i_over_sigma
        },
        "TwinL2": {
            "source": "Xtriage(Phenix)",
            "value": float("%.3f" % float(twin_l2))
        },
        "bulk_solvent_b": {
            "source": "EDS",
            "value": float("%.3f" % float(bulk_solvent_b))
        }
    }, 200
    

def get_validation_residuewise_outlier_summary(entry_id, graph):

    query1 = """
    MATCH(entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[chain_rel:IS_IN_CHAIN {OBSERVED:'Y'}]->(chain:Chain)
    MATCH (pdb_res)-[val_rel:HAS_VALIDATION_DATA]->(val) WHERE chain.STRUCT_ASYM_ID=val_rel.STRUCT_ASYM_ID
    RETURN DISTINCT entity.ID, pdb_res.ID, chain_rel.AUTH_SEQ_ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, chain_rel.MODEL, labels(val)[0] ORDER BY toInteger(pdb_res.ID)
    """

    query2 = """
    MATCH (entry:Entry {ID:$entry_id})-[:HAS_ENTITY]->(entity:Entity)-[:HAS_PDB_RESIDUE]->(pdb_res:PDB_Residue)-[chain_rel:IS_IN_CHAIN]->(chain:Chain) WHERE chain_rel.RAMA='OUTLIER' OR 
    chain_rel.ROTA='OUTLIER' OR toFloat(chain_rel.RSRZ) > 2.0
    RETURN entity.ID, chain.AUTH_ASYM_ID, chain.STRUCT_ASYM_ID, chain_rel.MODEL, pdb_res.ID, chain_rel.AUTH_SEQ_ID, chain_rel.RAMA, chain_rel.ROTA, chain_rel.RSRZ
    ORDER BY toInteger(pdb_res.ID)
    """


    mappings = list(graph.run(query1, entry_id=entry_id)) + list(graph.run(query2, entry_id=entry_id))
    

    if(len(mappings) == 0):
        return {}, 404

    api_result = {
        "molecules": []
    }

    dict_entity = {}
    dict_chains = {}
    dict_model = {}
    dict_residues = {}

    entity_id, pdb_res_id, auth_seq_id, auth_asym_id, struct_asym_id, model, type_of_outlier, rama, rota, rsrz = (None for x in range(0,10))

    for mapping in mappings:

        print(mapping)

        types_of_outlier = []

        if(len(mapping) == 7):

            (entity_id, pdb_res_id, auth_seq_id, auth_asym_id, struct_asym_id, model, type_of_outlier) = mapping

            if(model is None): continue

            if(type_of_outlier == "Val_Clash"):
                types_of_outlier.append("clashes")
            elif(type_of_outlier == "Val_Chiral_Outlier"):
                types_of_outlier.append("chirals")
            elif(type_of_outlier == "Val_Bond_Outlier"):
                types_of_outlier.append("bond_lengths")
            elif(type_of_outlier == "Val_Angle_Outlier"):
                types_of_outlier.append("bond_angles")
            elif(type_of_outlier == "Val_Plane_Outlier"):
                types_of_outlier.append("planes")
            elif(type_of_outlier == "Val_Symm_Clash"):
                types_of_outlier.append("symm_clashes")

        elif(len(mapping) == 9):
            
            (entity_id, auth_asym_id, struct_asym_id, model, pdb_res_id, auth_seq_id, rama, rota, rsrz) = mapping

            if(rota == 'OUTLIER'):
                types_of_outlier.append("sidechain_outliers")
            if(rama == 'OUTLIER'):
                types_of_outlier.append("ramachandran_outliers")
            if(rsrz is not None and float(rsrz) > RSRZ_OUTLIER_CUTOFF):
                types_of_outlier.append("RSRZ")

        if(dict_entity.get(entity_id) is None):
            dict_entity[entity_id] = [(auth_asym_id, struct_asym_id)]
        else:
            dict_entity[entity_id].append((auth_asym_id, struct_asym_id))

        if(dict_chains.get((entity_id, auth_asym_id)) is None):
            dict_chains[(entity_id, auth_asym_id)] = [model]
        else:
            dict_chains[(entity_id, auth_asym_id)].append(model)

        if(dict_model.get((entity_id, auth_asym_id, model)) is None):
            dict_model[(entity_id, auth_asym_id, model)] = [(pdb_res_id, auth_seq_id)]
        else:
            dict_model[(entity_id, auth_asym_id, model)].append((pdb_res_id, auth_seq_id))

        if(dict_residues.get((entity_id, auth_asym_id, model, pdb_res_id)) is None):
            dict_residues[(entity_id, auth_asym_id, model, pdb_res_id)] = [x for x in types_of_outlier]
        else:
            dict_residues[(entity_id, auth_asym_id, model, pdb_res_id)] += [x for x in types_of_outlier]


    for entity_id in dict_entity.keys():
        entity_segment = {
            "entity_id": int(entity_id),
            "chains": []
        }
        for auth_asym_id, struct_asym_id in set(dict_entity[entity_id]):
            chain_element = {
                "chain_id": auth_asym_id,
                "struct_asym_id": struct_asym_id,
                "models": []
            }
            for model in set(dict_chains[(entity_id, auth_asym_id)]):
                model_element = {
                    "model_id": int(model),
                    "residues": []
                }
                for pdb_res_id, auth_seq_id in set(dict_model[(entity_id, auth_asym_id, model)]):
                    residue_element = {
                        "author_insertion_code": "",
                        "author_residue_number": int(auth_seq_id),
                        "alt_code": "",
                        "outlier_types": dict_residues[(entity_id, auth_asym_id, model, pdb_res_id)],
                        "residue_number": int(pdb_res_id)
                    }
                    model_element["residues"].append(residue_element)
                chain_element["models"].append(model_element)
            entity_segment["chains"].append(chain_element)
        api_result["molecules"].append(entity_segment)
        

    return api_result, 200

