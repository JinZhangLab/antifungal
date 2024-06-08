# Antifungal
This repository hosts the source code of python package [Antifungal](https://pypi.org/project/antifungal/) , which server as the backend of our online platform for antifungal antivity prediction and rational designs of antifungal peptides. For online use, visit [Antifungipept](https://antifungipept.chemoinfolab.com).

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Antifungal activity prediction](#antifungal-activity-prediction)
  - [Rational design methods](#rational-design-methods)
- [Directory Structure](#directory-structure)
- [Major Release Notes](#major-release-notes)
- [Citation](#citation)
- [License](#license)


## Installation

To install the package, it is recommended to use a Python 3.8 environment. Follow these steps:

```bash
# Create a virtual environment (optional but recommended)
conda create --name antifungal python=3.8
conda  activate antifungal

# install the antifungal with pip
pip install antifungal
```


## Usage
This package provide the functionalities, like antifungal activity prediction and the rational design of antifungal peptides. 

### Antifungal activity prediction
Antifungal activity prediction provides the probability of a peptide being antifungal and the minimum inhibitory concentration (MIC) against four fungal species: *Candida albicans*, *Candida krusei*, *Cryptococcus neoformans*, and *Candida parapsilosis*. The antifungal index (AFI) is also calculated to evaluate the overall antifungal activity of the peptide.

```python
from antifungal.predict import predict_MIC

seq = ['HIHIRHMWLLR','HIHIRHMWLLRR']
pred = predict_MIC(seq)
print(pred)
# Expected output: 
# {
#     'antifungal': [True, True],
#     'prob_antifungal': [95.2, 97.9],
#     'MIC_C_albicans': [21.8, 17.34],
#     'prob_MIC_C_albicans': [99.8, 99.8],
#     'MIC_C_krusei': [7.13, 5.87],
#     'prob_MIC_C_krusei': [99.3, 99.4],
#     'MIC_C_neoformans': [24.4, 15.57],
#     'prob_MIC_C_neoformans': [99.3, 99.6],
#     'MIC_C_parapsilosis': [18.3, 17.05],
#     'prob_MIC_C_parapsilosis': [84.5, 82.6],
#     'AFI': [16.23, 12.82],
#     'prob_AFI': [79.16, 79.9],
#     'peptide_seq': ['HIHIRHMWLLR', 'HIHIRHMWLLRR']
# }
```

### Rational design methods

We categorize the rational design of antifungal peptides into three types based on sequence modifications: increasing, preserving, and reducing sequence length.

1. **Sequence Length Increasing Methods**
  - **Augment**: Adds amino acids to the N-terminus, C-terminus, or both ends of the peptide sequence.
  - **Insert**: Inserts amino acids at specific positions within the peptide sequence.
  - **Duplicate**: Duplicates existing amino acids within the peptide sequence.

2. **Sequence Length Preserving Methods**
  - **Single Point Mutation**: Replaces individual amino acids with other amino acids.
  - **Swap**: Swaps positions of adjacent amino acids within the sequence.
  - **Shift**: Shifts the entire sequence to the left or right by specified offsets.
  - **Global Optimize**: Optimizes the peptide sequence by performing multiple amino acid mutations simultaneously.

3. **Sequence Length Reducing Methods**
  - **Segment**: Segments the peptide into shorter sub-sequences.
  - **Delete**: Deletes specific amino acids from the peptide sequence.

All methods utilize common functions for generating candidate sequences and predicting antifungal activity.

```python
# Demonstrative design using segment method
from antifungal.design import segment

segment_instance = segment("HIHIRHMWLLRRR")
segment_predictions = segment_instance.get_candidate_sequences().predict()
print(segment_predictions)
# Expected output: 
# {
#     'antifungal': [True, True, True, True, True, True],
#     'prob_antifungal': [95.2, 97.9, 98.2, 97.7, 97.5, 99.0],
#     'MIC_C_albicans': [21.8, 17.34, 21.03, 19.68, 25.91, 24.36],
#     'prob_MIC_C_albicans': [99.8, 99.8, 99.8, 99.8, 99.8, 99.8],
#     'MIC_C_krusei': [7.13, 5.87, 5.85, 6.45, 5.96, 4.7],
#     'prob_MIC_C_krusei': [99.3, 99.4, 99.7, 99.0, 99.5, 98.9],
#     'MIC_C_neoformans': [24.4, 15.57, 10.01, 16.36, 10.46, 16.2],
#     'prob_MIC_C_neoformans': [99.3, 99.6, 99.8, 99.6, 99.9, 99.5],
#     'MIC_C_parapsilosis': [18.3, 17.05, 17.08, 18.28, 18.21, 19.01],
#     'prob_MIC_C_parapsilosis': [84.5, 82.6, 85.6, 83.7, 82.4, 85.2],
#     'AFI': [16.23, 12.82, 12.04, 13.96, 13.1, 13.7],
#     'prob_AFI': [79.16, 79.9, 83.47, 80.47, 79.7, 82.84],
#     'peptide_seq': ['HIHIRHMWLLR', 'HIHIRHMWLLRR', 'HIHIRHMWLLRRR', 'IHIRHMWLLRR', 'IHIRHMWLLRRR', 'HIRHMWLLRRR'],
#     'seq_name': ['segment_1_11', 'segment_1_12', 'segment_1_13', 'segment_2_12', 'segment_2_13', 'segment_3_13']
# }
```

**Globally optimize**  
The globally_optimize class is a length-preserving globall optimization method that utilizes an evolutionary algorithm to perform multiple amino acid mutations simultaneously, enabling the improvement of antifungal activity of peptides.
```python
# Demonstrative usage of globally optimize method (this will take a few minutes to hours depending on the computational resources)
from antifungal.design import globally_optimize
optimization_instance = globally_optimize("HIHIRHMWLLRRR")
optimized_seq, results = optimization_instance.optimize()
print(results)
# Expected outputs:
# {
#     "optimized_seq": "FICFRCMWFCRRL",
#     "antifungal_idx": [3.96]
# }
```

## Directory Structure
- **data/**: Contains data used for model development.
  - **training_data/**: Stores the datasets utilized in training the predictive models.
  - **screening_data**: Contains data from extensive screening studies detailed in the referenced article.
- **model/**: Houses the trained models for antifungal peptide prediction.

- **propy/**: A slightily modified version of the package [propy](https://pypi.org/project/propy/) for the supports of more verstaile peptide sequences.

- **ChemoinfoPy/**: Contains Python scripts for variable selection, peptide sequence preprocessing, descriptor calculation, and dataset partitioning for correction and validation purposes.


## Major Release Notes
- **Version 0.1.0 (Dec 9, 2023)**: The initial release of the Antifungal package, providing fundamental functionalities for antifungal peptide prediction.
- **Version 0.1.2 (Apr 10, 2024)**: Support rational design of antifungal peptides.
- **Version 0.1.3 (June 8, 2024)**: Expanded the rational design capabilities to include up to nine distinct methods for optimizing antifungal peptides.


## Citation
If you use this tool, please cite our Publications:
1. J. Zhang, et al. Large-Scale Screening of Antifungal Peptides Based on Quantitative Structure–Activity Relationship. ACS Med. Chem. Lett. 2022, 13, 1, 99–104. [link](https://pubs.acs.org/doi/10.1021/acsmedchemlett.1c00556). 
2. J. Zhang, et al. *In Silico* Design and Synthesis of Antifungal Peptides Guided by Quantitative Antifungal Activity. J. Chem. Inf. Model. 2024, 64 (10), 4277–4285. [link](https://doi.org/10.1021/acs.jcim.4c00142)


## License
This project is licensed under the [MIT License](https://choosealicense.com/licenses/mit/), which allows for broad usage and modification under specific terms.