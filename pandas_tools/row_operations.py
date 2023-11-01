def prioritise_aa(row):
    aa_sub = row.MUTATION_AA
    if aa_sub == "p.?":
        return 1
    return 0
