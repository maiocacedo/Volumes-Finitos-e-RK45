import sympy as sp

def symbol_references(in_list):
    slist = []

    for e in in_list:
        globals()[e] = sp.Symbol(e)
        slist.append(e)
    return slist