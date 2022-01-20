
def og_parts_to_ogs(ogs, scores):
    """Convert to best-scores for each whole og rather than part og
    """
    ogs_full = [og.split(".")[0] for og in ogs]
    ogs_set = set(ogs_full)
    # index will get the first occurence (highest-scoring) for each full og
    # sort them to get them in score-order rather than arbitrary
    i_keep = list(sorted([ogs_full.index(og) for og in ogs_set]))
    ogs = [ogs_full[iog] for iog in i_keep]
    scores = [scores[iog] for iog in i_keep]
    return ogs, scores
