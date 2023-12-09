# Antifungal Peptide Prediction Tool

This repository hosts the Antifungal Peptide Prediction Tool, a Python package for predicting and analyzing antifungal peptides. It integrates various functionalities including peptide sequence processing, descriptor calculation, and machine learning-based prediction models.

## Installation

To install the package, follow these steps:

```bash
# Create a virtual environment (optional but recommended)
pip install virtualenv
virtualenv env --python=python3.8
env\Scripts\activate # For Windows
source env/bin/activate # For Linux

# install antifungal package with pip
pip install antifungal
```


## Usage
The Antifungal Peptide Prediction Tool serves multiple purposes including the prediction of antifungal activities and the rational design of peptides. It enables users to segment peptide sequences, perform single-point mutation analysis, and globally optimize peptide sequences for enhanced properties.

### Example Usage for antifungal activity prediction
```python
from antifungal.predict import predict_MIC

seq = ['HIHIRHMWLLR','HIHIRHMWLLRR']
pred = predict_MIC(seq)
print(pred)
# Expected output: 
{
    'antifungal': [True, True],
    'prob_antifungal': [95.2, 97.9],
    'MIC_C_albicans': [21.8, 17.34],
    'prob_MIC_C_albicans': [99.8, 99.8],
    'MIC_C_krusei': [7.13, 5.87],
    'prob_MIC_C_krusei': [99.3, 99.4],
    'MIC_C_neoformans': [24.4, 15.57],
    'prob_MIC_C_neoformans': [99.3, 99.6],
    'MIC_C_parapsilosis': [18.3, 17.05],
    'prob_MIC_C_parapsilosis': [84.5, 82.6],
    'AFI': [16.23, 12.82],
    'prob_AFI': [79.16, 79.9],
    'peptide_seq': ['HIHIRHMWLLR', 'HIHIRHMWLLRR']
}
```


### Example Usage for antifungal peptide design
```python
from antifungal.design import segment, single_point_mutation, global_optimization

# Example for segment class
segment_instance = segment("HIHIRHMWLLRRR")
segment_predictions = segment_instance.get_segmented_sequences().predict()
print(segment_predictions)
# Expected output: 
{
    'antifungal': [True, True, True, True, True, True],
    'prob_antifungal': [95.2, 97.9, 98.2, 97.7, 97.5, 99.0],
    'MIC_C_albicans': [21.8, 17.34, 21.03, 19.68, 25.91, 24.36],
    'prob_MIC_C_albicans': [99.8, 99.8, 99.8, 99.8, 99.8, 99.8],
    'MIC_C_krusei': [7.13, 5.87, 5.85, 6.45, 5.96, 4.7],
    'prob_MIC_C_krusei': [99.3, 99.4, 99.7, 99.0, 99.5, 98.9],
    'MIC_C_neoformans': [24.4, 15.57, 10.01, 16.36, 10.46, 16.2],
    'prob_MIC_C_neoformans': [99.3, 99.6, 99.8, 99.6, 99.9, 99.5],
    'MIC_C_parapsilosis': [18.3, 17.05, 17.08, 18.28, 18.21, 19.01],
    'prob_MIC_C_parapsilosis': [84.5, 82.6, 85.6, 83.7, 82.4, 85.2],
    'AFI': [16.23, 12.82, 12.04, 13.96, 13.1, 13.7],
    'prob_AFI': [79.16, 79.9, 83.47, 80.47, 79.7, 82.84],
    'peptide_seq': ['HIHIRHMWLLR', 'HIHIRHMWLLRR', 'HIHIRHMWLLRRR', 'IHIRHMWLLRR', 'IHIRHMWLLRRR', 'HIRHMWLLRRR'],
    'seq_name': ['segment_1_11', 'segment_1_12', 'segment_1_13', 'segment_2_12', 'segment_2_13', 'segment_3_13']
}

# Example for single_point_mutation class
mutation_instance = single_point_mutation("HIHIRHMWLLRRR")
mutation_predictions = mutation_instance.get_mutated_sequences().predict()
print(mutation_predictions)
# Expected output for single_point_mutation:
{
    "antifungal": [true, true, ...], 
    "prob_antifungal": [100.0, 95.8, ...], 
    "MIC_C_albicans": [28.24, 23.48, ...], 
    "prob_MIC_C_albicans": [99.9, 99.8, ...], 
    "MIC_C_krusei": [8.37, 8.57, ...], 
    "prob_MIC_C_krusei": [99.9, 99.9, ...], 
    "MIC_C_neoformans": [6.58, 5.2, ...], 
    "prob_MIC_C_neoformans": [99.7, 99.2, ...], 
    "MIC_C_parapsilosis": [27.36, 23.58, ...], 
    "prob_MIC_C_parapsilosis": [86.9, 87.2, ...], 
    "AFI": [14.37, 12.54, ...],
    "prob_AFI": [86.47, 82.62, ...], 
    "peptide_seq": ["AIHIRHMWLLRRR", "CIHIRHMWLLRRR", ...], 
    "seq_name": ["mutate_1_A", "mutate_1_C", ...] 
}


# Example for global_optimization class, this will take a few minutes to hours to run depending on the number of iterations and sequence length
optimization_instance = global_optimization("HIHIRHMWLLRRR")
optimized_seq, results = optimization_instance.optimize()
print(results)
# Expected output for global_optimization:
{
    "optimized_seq": "FICFRCMWFCRRL",
    "antifungal_idx": [3.96]
}
```

## Directory Structure
- **data/**: Contains data used for model development.
  - **training_data/**: Stores the datasets utilized in training the predictive models.
  - **screening_data**: Contains data from extensive screening studies detailed in the referenced article.
- **model/**: Houses the trained models for antifungal peptide prediction.

- **propy/**: A modified version of the propy package, optimized for enhanced performance and bug fixes. The original package can be found at [propy](https://pypi.org/project/propy/).

- **ChemoinfoPy/**: Contains Python scripts for variable selection, peptide sequence preprocessing, descriptor calculation, and dataset partitioning for correction and validation purposes.


## Reference
For a comprehensive understanding and methodological details, refer to the paper "Large-Scale Screening of Antifungal Peptides Based on Quantitative Structure–Activity Relationship," published in ACS Med. Chem. Lett. (2022, 13, 1, 99–104)[link]](https://pubs.acs.org/doi/10.1021/acsmedchemlett.1c00556). For additional resources and online tools, visit the [Antifungal Webserver](https://www.chemoinfolab.com/antifungal).

## License
This project is licensed under the [MIT License](https://choosealicense.com/licenses/mit/), which allows for broad usage and modification under specific terms.