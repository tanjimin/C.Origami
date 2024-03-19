import numpy as np
from scipy.stats import pearsonr, spearmanr
import insulation_score as insu
from tqdm import tqdm

def mse(preds, targets):
    mse = ((preds - targets) ** 2).mean(axis = (1, 2))
    results = list(mse.astype(float))
    return results

def insulation_pearson(preds, targets):
    scores = []
    for pred, target in zip(preds, tqdm(targets)):
        pred_insu = np.array(insu.chr_score(pred))
        label_insu = np.array(insu.chr_score(target))
        nas = np.logical_or(np.isnan(pred_insu), np.isnan(label_insu))
        if nas.sum() == len(pred):
            scores.append(np.nan)
        else:
            metric, p_val = pearsonr(pred_insu[~nas], label_insu[~nas])
            scores.append(metric)
    results = scores
    return results

def observed_vs_expected_with_means(preds, targets, preds_mean, targets_mean):
    scores = []
    for pred, target in zip(preds - preds_mean, tqdm(targets - targets_mean)):
        metric, p_val = pearsonr(pred.reshape(-1), target.reshape(-1))
        scores.append(metric)
    results = scores
    return results

def observed_vs_expected(preds, targets):
    scores = []
    preds_mean = preds.mean(axis = 0, keepdims = True)
    targets_mean = targets.mean(axis = 0, keepdims = True)
    for pred, target in zip(preds - preds_mean, tqdm(targets - targets_mean)):
        metric, p_val = pearsonr(pred.reshape(-1), target.reshape(-1))
        scores.append(metric)
    results = scores
    return results

def distance_stratified_correlation(preds, targets):
    scores = []
    for pred, target in zip(preds, tqdm(targets)):
        distance_list = []
        for dis_i in range(len(pred)):
            pred_diag_i = np.diagonal(pred, offset = dis_i)
            target_diag_i = np.diagonal(target, offset = dis_i)
            if len(pred_diag_i) < 2: break
            metric, p_val = pearsonr(pred_diag_i, target_diag_i)
            distance_list.append(metric)
        scores.append(distance_list)
    results = scores
    return results
