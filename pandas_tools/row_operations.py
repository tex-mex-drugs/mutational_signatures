def prioritise_exon_aa(row):
    aa_sub = row.MUTATION_AA
    if aa_sub == "p.?":
        return 1
    return 0


# Get rid of mutations that aren't of interest, i.e. p.?, synonymous, and stop lost
def filter_mutation_aa(row):
    mutation_aa = row["MUTATION_AA"]
    # Filter out stop lost mutations
    if "ext" in mutation_aa:
        return False
    # Filter out anonymous mutations and intron variants
    if mutation_aa[-1] == "=" or mutation_aa[-1] == "?":
        return False
    return True
