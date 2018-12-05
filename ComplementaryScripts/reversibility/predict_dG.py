from equilibrator_api import ComponentContribution, Reaction
from cobra.io import read_sbml_model

def predict_dG(model, pH = 7, ionic_strength = 0.1):
    eq_api = ComponentContribution(pH=7.0, ionic_strength = 0.1)
    lst = []
    for r in model.reactions:
        kegg_string = build_kegg_string(r)
        if kegg_string:
            try:
                rxn = Reaction.parse_formula(kegg_string)
            except KeyError:
                print("eQuilibrator could not predict dG for {0}".format(r.id))
                continue
            # eq_api.dG0_prime(rxn)
            try:
                dgm = eq_api.dGm_prime(rxn)
                dg0 = eq_api.dG0_prime(rxn)
            except:
                print("eQuilibrator could not predict dG for {0}".format(r.id))
            else:
                lst.append([r.id, dg0, dgm])

    return lst

def build_kegg_string(r):
    mass_balance = r.check_mass_balance()
    try:
        r.annotation["kegg.reaction"]
    except KeyError:
        return None
    
    if len(mass_balance):
        print("{0} is not mass balanced".format(r.id))
        return None


    reactants = []
    products = []
    for m, coeff in r.metabolites.items():
        try:
            kegg_id = m.annotation["kegg.compound"]
        except KeyError:
            return None
        else:
            if isinstance(kegg_id, list):
                print("Double kegg annotation: ", kegg_id, " Use first only")
                kegg_id = kegg_id[0]
            if coeff > 0:
                products.append(kegg_id)
            else:
                reactants.append(kegg_id)



    reactant_string = " + ".join(reactants)
    product_string = " + ".join(products)
    return reactant_string + " <=> " + product_string

    
    
if __name__ == '__main__':
    scoGEM = read_sbml_model("../../ModelFiles/xml/scoGEM.xml")
    predict_dG(scoGEM)
