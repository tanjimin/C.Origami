import numpy as np
import matplotlib.pyplot as plt

def chr_score(matrix, res = 10000, radius = 500000, pseudocount_coeff = 30):
    pseudocount = matrix.mean() * pseudocount_coeff
    pixel_radius = int(radius / res)
    scores = []
    for loc_i, loc in enumerate(range(len(matrix))):
        scores.append(point_score(loc, pixel_radius, matrix, pseudocount))
    return scores

def point_score(locus, radius, matrix, pseudocount):
    l_edge = max(locus - radius, 0)
    r_edge = min(locus + radius, len(matrix))
    l_mask = matrix[l_edge : locus, l_edge : locus]
    r_mask = matrix[locus : r_edge, locus : r_edge]
    center_mask = matrix[l_edge : locus, locus : r_edge]
    score = (max(l_mask.mean(), r_mask.mean()) +  pseudocount) /\
            (center_mask.mean() + pseudocount)
    return score
