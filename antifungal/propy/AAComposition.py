# -*- coding: utf-8 -*-
"""
The module is used for computing the composition of amino acids, dipetide and
3-mers (tri-peptide) for a given protein sequence.

References
----------
.. [1] Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
   fold class predictions. Nucleic Acids Res, 22, 3616-3619.

.. [2] Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
   subcellular localization prediction. Bioinformatics, 17, 721-728.

.. [3] Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold
   class prediction: new methods of statistical classification. Proc Int Conf
   Intell Syst Mol Biol, 106-112.

Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.3.27
Email: oriental-cds@163.com
"""

# Core Library
import re
from typing import Any, Dict

AALetter = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]


def CalculateAAComposition(ProteinSequence):
    """
    Calculate the composition of Amino acids for a given protein sequence.

    Parameters
    ----------
    ProteinSequence: str
        a pure protein sequence

    Returns
    -------
    result : Dict
        contains the composition of 20 amino acids.

    Examples
    --------
    >>> result = CalculateAAComposition(protein)
    """
    LengthSequence = len(ProteinSequence)
    result = {}
    for i in AALetter:
        result[i] = round(float(ProteinSequence.count(i)) / LengthSequence * 100, 3)
    return result


def CalculateDipeptideComposition(ProteinSequence):
    """
    Calculate the composition of dipeptidefor a given protein sequence.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence

    Returns
    -------
    result : Dict
        contains the composition of 400 dipeptides

    Examples
    --------
    >>> result = CalculateDipeptideComposition(protein)
    """
    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        for j in AALetter:
            Dipeptide = i + j
            Result[Dipeptide] = round(
                float(ProteinSequence.count(Dipeptide)) / (LengthSequence - 1) * 100, 2
            )
    return Result


def Getkmers():
    """
    Get the amino acid list of 3-mers.

    Returns
    -------
    result : List
        contains 8000 tri-peptides

    Examples
    --------
    >>> result = Getkmers()
    """
    kmers = list()
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                kmers.append(i + j + k)
    return kmers


def GetSpectrumDict(proteinsequence):
    """
    Calcualte the spectrum descriptors of 3-mers for a given protein.

    Parameters
    ----------
    proteinsequence : a pure protein sequence

    Returns
    -------
    result : Dict
        contains the composition values of 8000 3-mers

    Examples
    --------
    >>> result = GetSpectrumDict(protein)
    """
    result = {}
    kmers = Getkmers()
    for i in kmers:
        result[i] = len(re.findall(i, proteinsequence))
    return result


def CalculateAADipeptideComposition(ProteinSequence):
    """
    Calculate the composition of AADs, dipeptide and 3-mers for a given protein
    sequence.

    Parameters
    ----------
    ProteinSequence : a pure protein sequence

    Returns
    -------
    result : Dict
        contains all composition values of AADs, dipeptide and 3-mers (8420).

    Examples
    --------
    >>> result = CalculateAADipeptideComposition(protein)
    """
    result: Dict[Any, Any] = {}
    result.update(CalculateAAComposition(ProteinSequence))
    result.update(CalculateDipeptideComposition(ProteinSequence))
    result.update(GetSpectrumDict(ProteinSequence))

    return result
