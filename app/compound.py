

def get_compound_atoms(chem_comp_id, graph):

    query = """
    MATCH (ligand:Ligand {ID:$chem_comp_id})-[:HAS_ATOM]->(atom:Atom)
    RETURN atom.STEREO_CONFIG, atom.LEAVING_ATOM_FLAG, atom.ALT_ATOM_ID, atom.AROMATIC_FLAG, atom.ELEMENT_SYMBOL, atom.MODEL_CARTN_X_IDEAL, atom.MODEL_CARTN_Y_IDEAL, 
    atom.MODEL_CARTN_Z_IDEAL, atom.CHARGE, atom.ID
    """
    
    mappings = list(graph.run(query, chem_comp_id=chem_comp_id))

    if(len(mappings) == 0):
        return {}, 404

    result_list = []

    for mapping in mappings:
        (stereo_config, leaving_atom_flag, alt_atom_id, aromatic_flag, element_symbol, x_ideal, y_ideal, z_ideal, charge, atom_id) = mapping

        if(leaving_atom_flag == 'N'):
            leaving_atom_flag = False
        else:
            leaving_atom_flag = True

        if(aromatic_flag == 'N'):
            aromatic_flag = False
        else:
            aromatic_flag = True

        if(stereo_config == 'N'):
            stereo_config = False
        else:
            stereo_config = True

        result_list.append({
            "stereo": stereo_config,
            "leaving_atom": leaving_atom_flag,
            "pdb_name": alt_atom_id,
            "aromatic": aromatic_flag,
            "element": element_symbol,
            "ideal_y": None if y_ideal is None else float("%.3f" % float(y_ideal)),
            "ideal_x": None if x_ideal is None else float("%.3f" % float(x_ideal)),
            "charge": None if charge is None else float(charge),
            "ideal_z": None if z_ideal is None else float("%.3f" % float(z_ideal)),
            "atom_name": atom_id
        })

    
    return result_list, 200


def get_compound_bonds(chem_comp_id, graph):

    query = """
    MATCH (ligand:Ligand {ID:$chem_comp_id})-[:HAS_ATOM]->(atom1:Atom)-[relation:BONDS]->(atom2:Atom)
    RETURN relation.STEREO_CONFIG, relation.AROMATIC_FLAG, relation.ORDER_NUMBER, relation.VALUE_ORDER, relation.VALUE_DIST_IDEAL, atom1.ID, atom2.ID
    """

    mappings = list(graph.run(query, chem_comp_id=chem_comp_id))
    
    if(len(mappings) == 0):
        return {}, 404
    
    result_list = []

    for mapping in mappings:

        (stereo_config, aromatic_flag, bond_order, bond_type, ideal_length, atom1, atom2) = mapping

        if(stereo_config == 'N'):
            stereo_config = False
        else:
            stereo_config = True

        if(aromatic_flag == 'N'):
            aromatic_flag = False
        else:
            aromatic_flag = True

        result_list.append({
            "stereo": stereo_config,
            "atom_1": atom1,
            "atom_2": atom2,
            "bond_type": bond_type,
            "bond_order": float(bond_order),
            "ideal_length": float("%.3f" % float(ideal_length)),
            "aromatic": aromatic_flag
        })

    return result_list, 200


def get_compound_in_pdb(chem_comp_id, graph):

    query = """
    MATCH (entity:Entity {POLYMER_TYPE:'B', CHEM_COMP_LIST:$chem_comp_id})<-[:HAS_ENTITY]->(entry:Entry)
    RETURN collect(entry.ID) as entries
    """

    mappings = list(graph.run(query, chem_comp_id=chem_comp_id))


    if(len(mappings) != 0):
        return mappings[0][0], 200
    else:
        return {}, 404


def get_compound_co_factors(graph):

    query = """
    MATCH (co:CO_Factor_Class)-[:HAS_EC]->(ec:EC), (co)<-[:ACTS_AS_COFACTOR]-(ligand:Ligand)
    RETURN co.COFACTOR_ID, co.NAME, collect(DISTINCT split(ec.EC_NUM,' ')[1]) as ec, collect(DISTINCT ligand.ID) as cofactors
    """

    api_result = {}

    mappings = list(graph.run(query))

    if(len(mappings) == 0):
        return {}, 404

    for mapping in mappings:
        (cofactor_class, cofactor_class_name, list_ec, list_cofactors) = mapping

        api_result[cofactor_class_name] = [{
            "EC": list_ec,
            "cofactors": list_cofactors
        }]

    return api_result, 200

def get_compound_co_factors_het(het_code, graph):

    query = """
    MATCH (src_ligand:Ligand {ID:$het_code})-[:ACTS_AS_COFACTOR]->(co:CO_Factor_Class)
    OPTIONAL MATCH (co)<-[:ACTS_AS_COFACTOR]->(dest_ligand:Ligand)
    RETURN co.NAME, dest_ligand.ID, dest_ligand.NAME
    """

    api_result = {
        "acts_as": "cofactor",
        "chem_comp_ids": []
    }
    mappings = list(graph.run(query, het_code=het_code))

    if(len(mappings) == 0):
        return None, 404

    for mapping in mappings:
        (co_name, dest_ligand_id, dest_ligand_name) = mapping

        api_result["class"] = co_name
        api_result["chem_comp_ids"].append({
            "chem_comp_id": dest_ligand_id,
            "name": dest_ligand_name
        })

    return api_result, 200