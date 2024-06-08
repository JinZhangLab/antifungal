# antifungal/predict.py
import sys
import os

import numpy as np
import pandas as pd

from .ChemoinfoPy.ProteinTools import calc_Pro_Des_values, seq_preprocess
from .ChemoinfoPy.FeatureSelection import pls
import pickle
from importlib import import_module

antifungal_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if antifungal_path not in sys.path:
    sys.path.append(antifungal_path)

sys.modules['ChemoinfoPy'] = import_module('antifungal.ChemoinfoPy')

base_dir = os.path.dirname(os.path.abspath(__file__))

# transfer of the rdist to probability
def transfer_rdist_prob(x):
    """
    Converts the relative distance (rdist) of regression result to a probability value.

    Parameters:
    x (float): The relative distance value.

    Returns:
    float: The probability value derived from the relative distance.
    """
    return float(np.round((1 - np.tanh(x)) * 100, 1))


def calc_activity(X):
    """
    Calculates the antifungal activity and related probabilities for a given set of peptide descriptors.

    This function loads predictive models for various Candida species and uses them
    to predict the Minimum Inhibitory Concentration (MIC) and associated probabilities.

    Parameters:
    X (numpy.array): Array of peptide descriptors.

    Returns:
    dict: A dictionary containing the predicted classification, probabilities, MIC values, and their probabilities for each species.
    """
    Data_myseq = dict()
    with open(os.path.join(base_dir, 'models','cmodel.pkl'), 'rb') as file:
        cScaler = pickle.load(file)
        c_is_selected = pickle.load(file)
        cSvcModel = pickle.load(file)

    Xc = cScaler.transform(X)
    yhat_c = cSvcModel.predict(Xc[:, c_is_selected])
    yhat_c_proba = cSvcModel.predict_proba(Xc[:, c_is_selected])[:, 1]

    Data_myseq['antifungal'] = [bool(y) for y in yhat_c]
    Data_myseq['prob_antifungal'] = [np.round(value * 100, 1) for value in yhat_c_proba]
    
    with open(os.path.join(base_dir, 'models','rmodel_C_a.pkl'), 'rb') as file:
        rScaler_C_a = pickle.load(file)
        is_fs_VIP_C_a = pickle.load(file)
        plsModel_C_a = pickle.load(file)
        CovModel_C_a = pickle.load(file)
        max_Mdist_C_a = pickle.load(file)
        svrModel_C_a = pickle.load(file)

    Xc_C_a = rScaler_C_a.transform(X)
    pyhat_C_a = svrModel_C_a.predict(Xc_C_a[:, is_fs_VIP_C_a])
    yhat_C_a = 2 ** -pyhat_C_a
    Mdist_C_a = CovModel_C_a.mahalanobis(
        plsModel_C_a.project(Xc_C_a[:, is_fs_VIP_C_a]))

    rMdist_C_a = Mdist_C_a / max_Mdist_C_a
    Data_myseq['MIC_C_albicans'] = [round(value, 2) for value in yhat_C_a]
    Data_myseq['prob_MIC_C_albicans'] = [transfer_rdist_prob(value) for value in rMdist_C_a]
    with open(os.path.join(base_dir, 'models','rmodel_C_k.pkl'), 'rb') as file:
        rScaler_C_k = pickle.load(file)
        is_fs_VIP_C_k = pickle.load(file)
        plsModel_C_k = pickle.load(file)
        CovModel_C_k = pickle.load(file)
        max_Mdist_C_k = pickle.load(file)
        svrModel_C_k = pickle.load(file)

    Xc_C_k = rScaler_C_k.transform(X)
    pyhat_C_k = svrModel_C_k.predict(Xc_C_k[:, is_fs_VIP_C_k])
    yhat_C_k = 2 ** -pyhat_C_k
    Mdist_C_k = CovModel_C_k.mahalanobis(
        plsModel_C_k.project(Xc_C_k[:, is_fs_VIP_C_k]))
    rMdist_C_k = Mdist_C_k / max_Mdist_C_k
    Data_myseq['MIC_C_krusei'] = [round(value, 2) for value in yhat_C_k]
    Data_myseq['prob_MIC_C_krusei'] = [transfer_rdist_prob(value) for value in rMdist_C_k]
    with open(os.path.join(base_dir, 'models','rmodel_C_n.pkl'), 'rb') as file:
        rScaler_C_n = pickle.load(file)
        is_fs_VIP_C_n = pickle.load(file)
        plsModel_C_n = pickle.load(file)
        CovModel_C_n = pickle.load(file)
        max_Mdist_C_n = pickle.load(file)
        svrModel_C_n = pickle.load(file)

    Xc_C_n = rScaler_C_n.transform(X)
    pyhat_C_n = svrModel_C_n.predict(Xc_C_n[:, is_fs_VIP_C_n])
    yhat_C_n = 2 ** -pyhat_C_n
    Mdist_C_n = CovModel_C_n.mahalanobis(
        plsModel_C_n.project(Xc_C_n[:, is_fs_VIP_C_n]))
    rMdist_C_n = Mdist_C_n / max_Mdist_C_n
    Data_myseq['MIC_C_neoformans'] = [round(value, 2) for value in yhat_C_n]
    Data_myseq['prob_MIC_C_neoformans'] = [transfer_rdist_prob(value) for value in rMdist_C_n]
    with open(os.path.join(base_dir, 'models','rmodel_C_p.pkl'), 'rb') as file:
        rScaler_C_p = pickle.load(file)
        is_fs_VIP_C_p = pickle.load(file)
        plsModel_C_p = pickle.load(file)
        CovModel_C_p = pickle.load(file)
        max_Mdist_C_p = pickle.load(file)
        svrModel_C_p = pickle.load(file)

    Xc_C_p = rScaler_C_p.transform(X)
    pyhat_C_p = svrModel_C_p.predict(Xc_C_p[:, is_fs_VIP_C_p])
    yhat_C_p = 2 ** -pyhat_C_p
    Mdist_C_p = CovModel_C_p.mahalanobis(
        plsModel_C_p.project(Xc_C_p[:, is_fs_VIP_C_p]))
    rMdist_C_p = Mdist_C_p / max_Mdist_C_p
    Data_myseq['MIC_C_parapsilosis'] = [round(value, 2) for value in yhat_C_p]
    Data_myseq['prob_MIC_C_parapsilosis'] = [transfer_rdist_prob(value) for value in rMdist_C_p]
    
    afi = 2 ** (-1 * (pyhat_C_a + pyhat_C_k + pyhat_C_n + pyhat_C_p) / 4)
    afi_proba = [a * b * c * d * e for a, b, c, d, e in zip(Data_myseq['prob_antifungal'], Data_myseq['prob_MIC_C_albicans'], Data_myseq['prob_MIC_C_krusei'], Data_myseq['prob_MIC_C_neoformans'], Data_myseq['prob_MIC_C_parapsilosis'])]
    afi_proba = [value / 100 / 100 / 100 / 100 for value in afi_proba]

    Data_myseq["AFI"] = [round(float(value), 2) for value in afi]
    Data_myseq["prob_AFI"] = [round(float(value), 2) for value in afi_proba]
    return Data_myseq

def calc_descriptor(seq):
    """
    Calculates the peptide descriptors for a given sequence or list of sequences.

    This function processes the input sequence(s) and calculates the relevant descriptors.

    Parameters:
    seq (str or list of str): A single peptide sequence or a list of peptide sequences.

    Returns:
    des (numpy.array): Array of peptide descriptors.
    """
    # If the input only contains a single peptide, make it a list
    if not isinstance(seq, list):
        seq = [seq]

    des = [calc_Pro_Des_values(seq_preprocess(seqi)) for seqi in seq]
    des = np.array(des)
    return des


def predict_MIC(seq):
    """
    Predicts antifungal activity for a given sequence or list of sequences.

    This function processes the input sequence(s), calculates the relevant descriptors, and then uses
    these descriptors to predict antifungal activity using the `calc_activity` function.

    Parameters:
    seq (str or list of str): A single peptide sequence or a list of peptide sequences.

    Returns:
    dict: A dictionary containing predicted activities, probabilities, and the original peptide sequences.
    """
    des = calc_descriptor(seq)
    peptide = calc_activity(des)
    peptide['peptide_seq'] = seq
    return peptide


if __name__ == "__main__":
   seq = ['HIHIRHMWLLR','HIHIRHMWLLRR']
   pred = predict_MIC(seq)
   print(pred)
   