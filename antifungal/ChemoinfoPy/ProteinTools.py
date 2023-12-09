# -*- coding: utf-8 -*-
"""
calculatin decriptor of multiple protein based on aa sequence

powered by jin zhang
"""
import sys
from ..propy import PyPro
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor

n_jobs = -1

standAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L",
           "K", "M", "F", "P", "S", "T", "V", "W", "Y"]  # Standard amino acids


def seq_preprocess(mySeq):
    mySeq = str(mySeq)
    mySeq = mySeq.strip()
    mySeq = mySeq.upper()
    return mySeq


def replace_unstand_AA(mySeq):
    mySeq = mySeq.replace('X', 'G')
    mySeq = mySeq.replace('Z', 'K')
    mySeq = mySeq.replace('B', 'L')
    mySeq = mySeq.replace('U', 'S')
    return mySeq


def is_standAA(mySeq):
    mySeq = seq_preprocess(mySeq)
    num_AA = [mySeq.count(a) for a in standAA]
    if sum(num_AA) == len(mySeq):
        return True
    else:
        return False


def calc_Pro_Des(mySeq):
    mySeq = seq_preprocess(mySeq)
    data = dict()
    # calculating descriptor based on modlamp

    DesObject = GlobalDescriptor(mySeq)
    # GlobalDescriptor in modlamp inlude as follows:
    # ['Length', 'MW', 'ChargeDensity', 'pI', 'InstabilityInd', 'Aromaticity', 'AliphaticInd', 'BomanInd', 'HydRatio']
    DesObject.calculate_all()
    data.update(zip(DesObject.featurenames, DesObject.descriptor.tolist()[0]))

    # =============================================================================
    # PeptideDescriptor in modlamp inlude as follows:
    # AASI (An amino acid selectivity index scale for helical antimicrobial peptides, [1] D. Juretić, D. Vukicević, N. Ilić, N. Antcheva, A. Tossi, J. Chem. Inf. Model. 2009, 49, 2873–2882.)
    # ABHPRK (modlabs inhouse physicochemical feature scale (Acidic, Basic, Hydrophobic, Polar, aRomatic, Kink-inducer)
    # argos (Argos hydrophobicity amino acid scale, [2] Argos, P., Rao, J. K. M. & Hargrave, P. A., Eur. J. Biochem. 2005, 128, 565–575.)
    # bulkiness (Amino acid side chain bulkiness scale, [3] J. M. Zimmerman, N. Eliezer, R. Simha, J. Theor. Biol. 1968, 21, 170–201.)
    # charge_phys (Amino acid charge at pH 7.0 - Hystidine charge +0.1.)
    # charge_acid (Amino acid charge at acidic pH - Hystidine charge +1.0.)
    # cougar (modlabs inhouse selection of global peptide descriptors)
    # eisenberg (the Eisenberg hydrophobicity consensus amino acid scale, [4] D. Eisenberg, R. M. Weiss, T. C. Terwilliger, W. Wilcox, Faraday Symp. Chem. Soc. 1982, 17, 109.)
    # Ez (potential that assesses energies of insertion of amino acid side chains into lipid bilayers, [5] A. Senes, D. C. Chadi, P. B. Law, R. F. S. Walters, V. Nanda, W. F. DeGrado, J. Mol. Biol. 2007, 366, 436–448.)
    # flexibility (amino acid side chain flexibilitiy scale, [6] R. Bhaskaran, P. K. Ponnuswamy, Int. J. Pept. Protein Res. 1988, 32, 241–255.)
    # grantham (amino acid side chain composition, polarity and molecular volume, [8] Grantham, R. Science. 185, 862–864 (1974).)
    # gravy (GRAVY hydrophobicity amino acid scale, [9] J. Kyte, R. F. Doolittle, J. Mol. Biol. 1982, 157, 105–132.)
    # hopp-woods (Hopp-Woods amino acid hydrophobicity scale,*[10] T. P. Hopp, K. R. Woods, Proc. Natl. Acad. Sci. 1981, 78, 3824–3828.*)
    # ISAECI (Isotropic Surface Area (ISA) and Electronic Charge Index (ECI) of amino acid side chains, [11] E. R. Collantes, W. J. Dunn, J. Med. Chem. 1995, 38, 2705–2713.)
    # janin (Janin hydrophobicity amino acid scale, [12] J. L. Cornette, K. B. Cease, H. Margalit, J. L. Spouge, J. A. Berzofsky, C. DeLisi, J. Mol. Biol. 1987, 195, 659–685.)
    # kytedoolittle (Kyte & Doolittle hydrophobicity amino acid scale, [13] J. Kyte, R. F. Doolittle, J. Mol. Biol. 1982, 157, 105–132.)
    # levitt_alpha (Levitt amino acid alpha-helix propensity scale, extracted from http://web.expasy.org/protscale. [14] M. Levitt, Biochemistry 1978, 17, 4277-4285.)
    # MSS (A graph-theoretical index that reflects topological shape and size of amino acid side chains, [15] C. Raychaudhury, A. Banerjee, P. Bag, S. Roy, J. Chem. Inf. Comput. Sci. 1999, 39, 248–254.)
    # MSW (Amino acid scale based on a PCA of the molecular surface based WHIM descriptor (MS-WHIM), extended to natural amino acids, [16] A. Zaliani, E. Gancia, J. Chem. Inf. Comput. Sci 1999, 39, 525–533.)
    # pepArc (modlabs pharmacophoric feature scale, dimensions are: hydrophobicity, polarity, positive charge, negative charge, proline.)
    # pepcats (modlabs pharmacophoric feature based PEPCATS scale, [17] C. P. Koch, A. M. Perna, M. Pillong, N. K. Todoroff, P. Wrede, G. Folkers, J. A. Hiss, G. Schneider, PLoS Comput. Biol. 2013, 9, e1003088.)
    # polarity (Amino acid polarity scale, [18] J. M. Zimmerman, N. Eliezer, R. Simha, J. Theor. Biol. 1968, 21, 170–201.)
    # PPCALI (modlabs inhouse scale derived from a PCA of 143 amino acid property scales, [19] C. P. Koch, A. M. Perna, M. Pillong, N. K. Todoroff, P. Wrede, G. Folkers, J. A. Hiss, G. Schneider, PLoS Comput. Biol. 2013, 9, e1003088.)
    # refractivity (Relative amino acid refractivity values, [20] T. L. McMeekin, M. Wilensky, M. L. Groves, Biochem. Biophys. Res. Commun. 1962, 7, 151–156.)
    # t_scale (A PCA derived scale based on amino acid side chain properties calculated with 6 different probes of the GRID program, [21] M. Cocchi, E. Johansson, Quant. Struct. Act. Relationships 1993, 12, 1–8.)
    # TM_tend (Amino acid transmembrane propensity scale, extracted from http://web.expasy.org/protscale, [22] Zhao, G., London E. Protein Sci. 2006, 15, 1987-2001.)
    # z3 (The original three dimensional Z-scale, [23] S. Hellberg, M. Sjöström, B. Skagerberg, S. Wold, J. Med. Chem. 1987, 30, 1126–1135.)
    # z5 (The extended five dimensional Z-scale, [24] M. Sandberg, L. Eriksson, J. Jonsson, M. Sjöström, S. Wold, J. Med. Chem. 1998, 41, 2481–2491.)
    # =============================================================================

    # =============================================================================
    # calculating descriptor based on propy
    # The protein descriptors calculated by propy
    # AAC: amino acid composition descriptors (20)
    # DPC: dipeptide composition descriptors (400)
    # TPC: tri-peptide composition descriptors (8000)
    # MBauto: Normalized Moreau-Broto autocorrelation descriptors (depend on the given properties, the default is 240)
    # Moranauto: Moran autocorrelation descriptors(depend on the given properties, the default is 240)
    # Gearyauto: Geary autocorrelation descriptors(depend on the given properties, the default is 240)
    # CTD: Composition, Transition, Distribution descriptors (CTD) (21+21+105=147)
    # SOCN: sequence order coupling numbers (depend on the choice of maxlag, the default is 60)
    # QSO: quasi-sequence order descriptors (depend on the choice of maxlag, the default is 100)
    # PAAC: pseudo amino acid composition descriptors (depend on the choice of lamda, the default is 50)
    # APAAC: amphiphilic pseudo amino acid composition descriptors(depend on the choice of lamda, the default is 50)
    # Ref: https://pypi.org/project/propy3/
    # =============================================================================
    DesObject = PyPro.GetProDes(mySeq)
    # AAC: amino acid composition descriptors (20)
    data.update(DesObject.GetAAComp())

    # DPC: dipeptide composition descriptors (400)
    data.update(DesObject.GetDPComp())

    # TPC: tri-peptide composition descriptors (8000)
    data.update(DesObject.GetTPComp())

    # Composition Transition Distribution descriptors (147)
    data.update(DesObject.GetCTD())

    # Geary autocorrelation descriptors (240).
    data.update(DesObject.GetGearyAuto())

    # Moran autocorrelation descriptors (240).
    data.update(DesObject.GetMoranAuto())

    # MBauto: Normalized Moreau-Broto autocorrelation descriptors
    #        (depend on the given properties, the default is 240)
    data.update(DesObject.GetMoreauBrotoAuto())

    # PCCA Type I Pseudo amino acid composition descriptors (default is 30).
    data.update(DesObject.GetPAAC())

    # QSO Quasi sequence order descriptors default is 50.
    data.update(DesObject.GetQSO())

    # SOCN Sequence order coupling numbers default is 45.
    data.update(DesObject.GetSOCN())

    # SubSeq Obtain the sub sequences wit length 2*window+1, whose central point is ToAA.
    # data.update(DesObject.GetSubSeq())
    return data


def calc_Pro_Des_values(mySeq):
    return list(calc_Pro_Des(mySeq).values())


def calc_Pro_Des_keys(mySeq):
    return list(calc_Pro_Des(mySeq).keys())


'''def calc_protein_descriptors(sequence, global_descriptors=False, 
                             AA_comp=False, DP_comp=False, TP_comp=False, 
                             CTD=False, GearyAuto=False, MoranAuto=False, 
                             MoreauBrotoAuto=False, PAAC=False, QSO=False, 
                             SOCN=False):
    # =============================================================================
    # calculating descriptor based on propy
    # The protein descriptors calculated by propy
    # AAC: amino acid composition descriptors (20)
    # DPC: dipeptide composition descriptors (400)
    # TPC: tri-peptide composition descriptors (8000)
    # MBauto: Normalized Moreau-Broto autocorrelation descriptors (depend on the given properties, the default is 240)
    # Moranauto: Moran autocorrelation descriptors(depend on the given properties, the default is 240)
    # Gearyauto: Geary autocorrelation descriptors(depend on the given properties, the default is 240)
    # CTD: Composition, Transition, Distribution descriptors (CTD) (21+21+105=147)
    # SOCN: sequence order coupling numbers (depend on the choice of maxlag, the default is 60)
    # QSO: quasi-sequence order descriptors (depend on the choice of maxlag, the default is 100)
    # PAAC: pseudo amino acid composition descriptors (depend on the choice of lamda, the default is 50)
    # APAAC: amphiphilic pseudo amino acid composition descriptors(depend on the choice of lamda, the default is 50)
    # Ref: https://pypi.org/project/propy3/
    # =============================================================================

    data = {}
    if global_descriptors:
        DesObject = GlobalDescriptor(sequence)
        DesObject.calculate_all()
        data.update(zip(DesObject.featurenames, DesObject.descriptor.tolist()[0]))
    if AA_comp:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetAAComp())
    if DP_comp:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetDPComp())
    if TP_comp:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetTPComp())
    if CTD:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetCTD())
    if GearyAuto:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetGearyAuto())
    if MoranAuto:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetMoranAuto())
    if MoreauBrotoAuto:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetMoreauBrotoAuto())
    if PAAC:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetPAAC())
    if QSO:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetQSO())
    if SOCN:
        DesObject = PyPro.GetProDes(sequence)
        data.update(DesObject.GetSOCN())

    return data'''


def calculate_ten_features(sequence, const):
    DesObject = GlobalDescriptor(sequence)
    # ['MW', 'Charge', 'ChargeDensity', 'pI', 'InstabilityInd', 'Aromaticity', 'AliphaticInd', 'BomanInd', 'HydrophRatio']
    if const == 'MW':
        DesObject.calculate_MW(amide=True)
    elif const == 'Charge':
        DesObject.calculate_charge(ph=7.4, amide=True)
    elif const == 'ChargeDensity':
        DesObject.charge_density(ph=7.4, amide=True)
    elif const == 'pI':
        DesObject.isoelectric_point(amide=True)
    elif const == 'InstabilityInd':
        DesObject.instability_index()
    elif const == 'Aromaticity':
        DesObject.aromaticity()
    elif const == 'AliphaticInd':
        DesObject.aliphatic_index()
    elif const == 'BomanInd':
        DesObject.boman_index()
    elif const == 'HydrophRatio':
        DesObject.hydrophobic_ratio()
    else:
        raise ValueError(
            "The constraint is not valid, please choose from ['MW', 'Charge', 'ChargeDensity', 'pI', 'InstabilityInd', 'Aromaticity', 'AliphaticInd', 'BomanInd', 'HydrophRatio']")
    return DesObject.descriptor[0][0]