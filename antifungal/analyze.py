# antifungal/analyze.py
import warnings
import numpy as np
from .predict import predict_MIC
from sklearn.linear_model import LinearRegression
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_score, KFold
from scipy import stats
import matplotlib.pyplot as plt

class antifungal_contribution:
    def __init__(self, sequence):
        self.amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        self.sequence = sequence


    def get_amino_acid_contribution(self, n_sampling=1000, mutation_rate=0.3, n_bootstrap=1000, seed=42):
        sequence_length = len(self.sequence)
        np.random.seed(seed)

        if n_sampling < 10*sequence_length:
            n_sampling = 10*sequence_length
            warnings.warn('Number of samples is set to 10 times the sequence length.')
        
        X = np.ones((n_sampling, sequence_length))
        mutated_sequences = []

        for i in range(n_sampling):
            sampled_indices = np.random.choice(sequence_length, int(sequence_length * mutation_rate), replace=False)
            mut_seq = list(self.sequence)

            for index in sampled_indices:
                mut_seq[index] = np.random.choice([aa for aa in self.amino_acids if aa != self.sequence[index]])

            X[i, sampled_indices] = 0
            mutated_sequences.append(''.join(mut_seq))
        
        predictions = predict_MIC(mutated_sequences)
        AFI = np.array(predictions['AFI'])
        is_antifungal = np.array(predictions['antifungal'])
        
        AFI[~is_antifungal] = AFI[is_antifungal].max()
        y = -np.log2(AFI)

        # Determine the optimal number of components for PLS using 10-fold cross-validation
        best_score = -np.inf
        best_n_components = 0
        for n_components in range(1, min(sequence_length, 20) + 1):  # Limiting to 20 components or sequence length
            pls = PLSRegression(n_components=n_components)
            scores = cross_val_score(pls, X, y, cv=KFold(10, shuffle=True, random_state=seed))
            mean_score = scores.mean()
            if mean_score > best_score:
                best_score = mean_score
                best_n_components = n_components
        
        # Fit the final PLS model with the best number of components
        pls = PLSRegression(n_components=best_n_components)
        pls.fit(X, y)
        # contributions = pls.coef_.flatten()

        # Bootstrap for confidence intervals
        bootstrapped_coefs = []
        for _ in range(n_bootstrap):
            sample_indices = np.random.choice(np.arange(n_sampling), n_sampling, replace=True)
            X_resampled = X[sample_indices]
            y_resampled = y[sample_indices]
            pls.fit(X_resampled, y_resampled)
            bootstrapped_coefs.append(pls.coef_.flatten())
        
        bootstrapped_coefs = np.array(bootstrapped_coefs)

        contributions = np.mean(bootstrapped_coefs, axis=0)
        r2 = pls.score(X, y)

        self.bootstrapped_coefs = bootstrapped_coefs
        self.r2 = r2

        return contributions, r2
    
    def plot_contributions(self, sig_level = 0.05, ax = None, save_path = None):
        # Perform significance tests
        lower_bound = np.percentile(self.bootstrapped_coefs, 2.5, axis=0)
        upper_bound = np.percentile(self.bootstrapped_coefs, 97.5, axis=0)
        contributions = np.mean(self.bootstrapped_coefs, axis=0)

        significant_positive = (lower_bound > 0)
        significant_negative = (upper_bound < 0)
        not_significant = ~significant_positive & ~significant_negative

        if ax is None:
            fig, ax = plt.subplots()


        ax.bar(np.arange(len(contributions))[significant_positive], contributions[significant_positive], color='green', alpha=0.5, label='Positive')
        ax.bar(np.arange(len(contributions))[significant_negative], contributions[significant_negative], color='red', alpha=0.5, label='Negative')
        ax.bar(np.arange(len(contributions))[not_significant], contributions[not_significant], color='gray', alpha=0.5, label='Not significant')
        ax.errorbar(np.arange(len(contributions)), contributions, yerr=[contributions - lower_bound, upper_bound - contributions], fmt='none', color='black', capsize=5)


        ax.axhline(0, color='black', linestyle='--')
        ax.grid(axis='both', linestyle='--', alpha=0.5)
        ax.set_xticks(np.arange(len(contributions)))
        ax.set_xticklabels(list(self.sequence))
        ax.set_xlabel('Amino acid')
        ax.set_ylabel('Antifungal Contribution')
        ax.set_title(f'Antifungal Contribution tests, R2 = {self.r2:.2f}')

        ax.legend(loc='upper right')


        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')

        return ax


