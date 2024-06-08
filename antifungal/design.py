# antifungal/design.py
import numpy as np
import random
import pandas as pd
from .predict import predict_MIC
from .ChemoinfoPy.ProteinTools import calc_Pro_Des_values, calc_Pro_Des
from modlamp.descriptors import GlobalDescriptor
from scipy.optimize import minimize, NonlinearConstraint


# Design antifungal peptide by methods that preserve the sequence length, includinging single point muation, swap, shift and globally optimize. 
class single_point_mutate():
    """
    A class for generating and predicting the antifungal activity of all sequences with exhaustive single point mutations. This helps identify mutations that enhance antifungal activity.

    Attributes:
        standard_AA (list): List of standard amino acids.
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (DataFrame): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible mutated sequences.
        seq_name (list): List of names assigned to each mutated sequence, with the format 'mutate_{original_aa}_{position}_{mutated_aa}'.

    Methods:
        get_candidate_sequences(): Generates all possible single point mutations of the sequence.
        predict(): Predicts the antifungal activity of all mutated sequences.
    """
    def __init__(self, sequence):
        """
        Constructs all the necessary attributes for the single_point_mutation object.

        Args:
            original_sequence (str): The original peptide sequence for mutation analysis.
            original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        """
        self.standard_AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        if not isinstance(sequence, str):
            raise ValueError("Input sequence must be a string.")
        
        if not all(aa in self.standard_AA for aa in sequence):
            raise ValueError("Sequence contains invalid amino acids. Only the uppercase single letter expression of standard amino acids are allowed.")
        
        self.original_sequence = sequence
        self.original_sequence_activity = predict_MIC(sequence)

    def get_candidate_sequences(self):
        """
        Generates all possible mutated sequences through single point mutation of the original sequence.

        Returns:
            self: The instance itself with updated 'candidate_sequences' and 'seq_name' attributes.
        """
        candidate_sequences = []
        seq_name = []

        for index, amino_acid in enumerate(self.original_sequence):
            for mutation in self.standard_AA:
                if mutation != amino_acid:
                    mut_seq = list(self.original_sequence)
                    mut_seq[index] = mutation
                    candidate_sequences.append("".join(mut_seq))
                    seq_name.append(f"mutate_{amino_acid}_{index+1}_{mutation}")

        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self

    def predict(self):
        """
        Predicts the antifungal activity of all mutated sequences generated by 'get_candidate_sequences'.

        Returns:
            Dictonary: A dictonary containing predictions for each mutated sequence.
        """
        predictions = predict_MIC(self.candidate_sequences)
        predictions["seq_name"] = self.seq_name
        return predictions


class swap(single_point_mutate):
    """
    A class for generating and predicting the antifungal activity of all sequences with exhaustive amino acid swaps. This include all possible sequences that can be generated by swapping two adjoint amino acids in the original sequence.

    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible swapped sequences.
        seq_name (list): List of names assigned to each swapped sequence, with the format 'swap_{position1}_{AA1}_{position2}_{AA2}'.

    Methods:
        get_candidate_sequences(): Generates all possible amino acid swaps of the sequence.
        predict(): Predicts the antifungal activity of all swapped sequences.
    """

    def get_candidate_sequences(self):
        """
        Generates all possible swapped sequences through amino acid substitution of the original sequence.

        Returns:
            self: The instance itself with updated 'candidate_sequences', and 'seq_name' attributes.
        """
        candidate_sequences = []
        seq_name = []

        for index, AA in enumerate(self.original_sequence):
            if index < len(self.original_sequence) - 1:
                tmp_seq = list(self.original_sequence)
                tmp_seq[index], tmp_seq[index+1] = tmp_seq[index+1], tmp_seq[index]
                candidate_sequences.append("".join(tmp_seq))
                seq_name.append(f"swap_{index+1}_{AA}_{index+2}_{tmp_seq[index+1]}")

        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self



class shift(single_point_mutate):
    """
    A class for generating and predicting the antifungal activity of all sequences with exhaustive amino acid shifts. This includes all possible sequences that can be generated by shifting the sequence to the right or left by all possible numbers of offset. 

    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (DataFrame): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible shifted sequences.
        seq_name (list): List of names assigned to each shifted sequence, with the format 'shift_{direction}_{offset}'.

    Methods:
        get_candidate_sequences(): Generates all possible shifted sequences of the sequence.
        predict(): Predicts the antifungal activity of all shifted sequences.
    """

    def get_candidate_sequences(self):
        """
        Generates all possible shifted sequences by cyclically shifting the original sequence to the left or right.

        Returns:
            self: The instance itself with updated 'candidate_sequences', and 'seq_name' attributes.
        """
        candidate_sequences = []
        seq_name = []

        for direction in ['left', 'right']:
            for offset in range(1, len(self.original_sequence)):
                if direction == 'left':
                    shifted_seq = self.original_sequence[offset:] + self.original_sequence[:offset]
                else:
                    shifted_seq = self.original_sequence[-offset:] + self.original_sequence[:-offset]
                candidate_sequences.append(shifted_seq)
                seq_name.append(f"shift_{direction}_{offset}")

        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self

# Design anfungal peptide by methods that reduce sequence length, including segment and delete
class segment(single_point_mutate):
    """
    A class for for designing antifungal peptides by segmenting peptide sequences into various subpeptides 
    of different lengths. This allows for the prediction of antifungal activity to identify the most active fragments.

    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible mutated sequences.
        seq_name (list): List of names assigned to each mutated sequence, with the format 'segment_{start}_{end}'.

    Methods:
        get_candidate_sequences(): Segments the peptide into subpeptides.
        predict(): Predicts the antifungal activity of all subpeptides.
    """

    def __init__(self, sequence, min_size=11):
        """
        Constructs all the necessary attributes for the segment object.

        Args:
            sequence (str): The peptide sequence to be segmented.
            min_size (int, optional): Minimum size for each subpeptide segment. Default is 11.
        """
        super().__init__(sequence)
        self.min_size = min_size

    def get_candidate_sequences(self):
        """
        Segments the peptide into subpeptides with a minimum length of 'min_size'.

        Returns:
            self: The instance itself with updated 'subpeptides', 'start_end', and 'seq_name' attributes.
        """
        sequence = self.original_sequence
        min_size = self.min_size
        candidate_sequences = []
        seq_name = []
        for i in range(0, len(sequence) - min_size + 1):
            for j in range(i + min_size, len(sequence) + 1):
                if j == len(sequence):
                    candidate_sequences.append(sequence[i:])
                else:
                    candidate_sequences.append(sequence[i:j])
                seq_name.append(f"segment_{i+1}_{j}")
        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self


class delete(single_point_mutate):
    """
    A class for designing antifungal peptides by deleting one amino acid from the peptide sequence. 
    
    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible mutated sequences.
        seq_name (list): List of names assigned to each mutated sequence, with the format 'delete_{position}_{amino_acid}'.

    Methods:
        get_candidate_sequences(): Deletes each amino acid from the sequence.
        predict(): Predicts the antifungal activity of all deleted sequences.
    """

    def get_candidate_sequences(self):
        """
        Deletes each amino acid from the sequence.

        Returns:
            self: The instance itself with updated 'candidate_sequences', and 'seq_name' attributes.
        """
        sequence = self.original_sequence
        candidate_sequences = []
        seq_name = []
        for i in range(len(sequence)):
            candidate_sequences.append(sequence[:i] + sequence[i+1:])
            seq_name.append(f"delete_{i+1}_{sequence[i]}")
        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self


# Design antifungal peptide by methods that increse the sequence length, including insert, augment and duplicate.
class insert(single_point_mutate):
    """
    A class for designing antifungal peptides by inserting an amino acid into the peptide sequence.

    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible mutated sequences.
        seq_name (list): List of names assigned to each mutates, with the form of 'insert_{position}_{AA}'.

    Methods:
        get_candidate_sequences(): Inserts each amino acid into the sequence.
        predict(): Predicts the antifungal activity of all inserted sequences.
    """

    def get_candidate_sequences(self):
        """
        Inserts each amino acid into the sequence.

        Returns:
            self: The instance itself with updated 'candidate_sequences', and 'seq_name' attributes.
        """
        sequence = self.original_sequence
        candidate_sequences = []
        seq_name = []
        for i in range(len(sequence) + 1):
            for AA in self.standard_AA:
                candidate_sequences.append(sequence[:i] + AA + sequence[i:])
                seq_name.append(f"insert_{i+1}_{AA}")
        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self

class augment(single_point_mutate):
    """
    A class for designing antifungal peptides by adding an amino acid to the N-terminus, C-terminus, or both ends of the peptide sequence.

    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible augmented sequences.
        seq_name (list): List of names assigned to each augmented sequence, with the format 'augment_{position}_{AA}'.

    Methods:
        get_candidate_sequences(): Adds each amino acid to the sequence.
        predict(): Predicts the antifungal activity of all augmented sequences.
    """

    def get_candidate_sequences(self):
        """
        Adds each amino acid to the N-terminus, C-terminus, or both ends of the sequence.

        Returns:
            self: The instance itself with updated 'candidate_sequences' and 'seq_name' attributes.
        """
        sequence = self.original_sequence
        candidate_sequences = []
        seq_name = []

        # Adding to N-terminus
        for AA in self.standard_AA:
            candidate_sequences.append(AA + sequence)
            seq_name.append(f"augment_N_{AA}")

        # Adding to C-terminus
        for AA in self.standard_AA:
            candidate_sequences.append(sequence + AA)
            seq_name.append(f"augment_C_{AA}")

        # Adding to both N-terminus and C-terminus
        for AA in self.standard_AA:
            for BB in self.standard_AA:
                candidate_sequences.append(AA + sequence + BB)
                seq_name.append(f"augment_NC_{AA}_{BB}")

        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self


class duplicate(single_point_mutate):
    """
    A class for designing antifungal peptides by duplicating each amino acid in the peptide sequence.

    Attributes:
        original_sequence (str): The original peptide sequence.
        original_sequence_activity (dict): Predicted antifungal activity of the original sequence.
        candidate_sequences (list): List of all possible duplicated sequences.
        seq_name (list): List of names assigned to each duplicated sequence, with the format 'duplicate_{position}_{AA}'.

    Methods:
        get_candidate_sequences(): Duplicates each amino acid in the sequence.
        predict(): Predicts the antifungal activity of all duplicated sequences.
    """

    def get_candidate_sequences(self):
        """
        Duplicates each amino acid in the sequence.

        Returns:
            self: The instance itself with updated 'candidate_sequences' and 'seq_name' attributes.
        """
        sequence = self.original_sequence
        candidate_sequences = []
        seq_name = []

        for i in range(len(sequence)):
            duplicated_seq = sequence[:i] + sequence[i] + sequence[i] + sequence[i+1:]
            candidate_sequences.append(duplicated_seq)
            seq_name.append(f"duplicate_{i+1}_{sequence[i]}")

        self.candidate_sequences = candidate_sequences
        self.seq_name = seq_name
        return self


# Evolutionary algorithm based antifungal peptide design method, characterized by the simutaneous alteration of multiple amino acids in the peptide sequence.
class globally_optimize():
    """
    A class for performing global optimization on a peptide sequence to enhance its  antifungal activity.

    Attributes:
        sequence (str): The input peptide sequence.
        seed (int): Seed value for random number generation, used in optimization process.

    Methods:
        cost(seq_nums): Computes the cost (e.g., antifungal activity) of a sequence.
        constraint(seq_nums, const): Applies a specified constraint to a sequence.
        transfer_alphabeta_num(seq): Converts a sequence from alphabet to numerical representation.
        trasfer_num_alphabeta(numbers): Converts a sequence from numerical to alphabet representation.
        softmax(x): Applies the softmax function to an array.
        optimize(constraintName, constraint_lb_tol_percent, constraint_ub_tol_percent, maxiter): Performs global optimization with specified constraints and parameters.
    """

    def __init__(self, sequence, seed=0):
        """
        Initializes the global_optimization object with a sequence and a seed for random number generation.

        Args:
            sequence (str): The input peptide sequence.
            seed (int, optional): Seed value for random number generation. Defaults to 0.
        """
        self.sequence = sequence
        self.seed = seed

    def cost(self, seq_nums):
        """
        Calculates the cost (e.g., antifungal activity index) of a sequence in numerical form.

        Args:
            seq_nums (numpy.ndarray): Numerical representation of the peptide sequence.

        Returns:
            float: Cost of the sequence, representing its antifungal activity.
        """
        sequence = self.trasfer_num_alphabeta(seq_nums)
        sequence_activity = predict_MIC(sequence)
        return sequence_activity["AFI"]

    def constraint(self, seq_nums, const):
        """
        Applies a specified constraint to a sequence.

        Args:
            seq_nums (numpy.ndarray): Numerical representation of the sequence.
            const (str): The constraint to be applied (e.g., molecular weight 'MW', charge 'Charge').

        Returns:
            float: The value of the constraint for the sequence.

        Raises:
            ValueError: If an invalid constraint is specified.
        """
        seq = self.trasfer_num_alphabeta(seq_nums)
        DesObject = GlobalDescriptor(seq)
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

    def transfer_alphabeta_num(self, seq):
        """
        Converts a sequence from alphabet (amino acid) to numerical (one-hot encoding) representation.

        Args:
            seq (str): Alphabet representation of the sequence.

        Returns:
            numpy.ndarray: Numerical representation of the sequence.
        """
        # one hot encoding
        standard_AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        encoded_seq = []
        for i, seq_aa in enumerate(seq):
            encoded_aa = [1.0 if seq_aa == standard_aa else 0.0 for standard_aa in standard_AA]
            encoded_seq.append(encoded_aa)
        return np.array(encoded_seq).flatten()

    def trasfer_num_alphabeta(self, numbers):
        """
        Converts a sequence from numerical to alphabet representation.

        Args:
            numbers (numpy.ndarray): Numerical representation of the sequence.

        Returns:
            str: Alphabet representation of the sequence.
        """
        standard_AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        reshped_num = numbers.reshape(-1, len(standard_AA))
        seq = []
        for i in range(reshped_num.shape[0]):
            reshped_num[i, :] = self.softmax(reshped_num[i, :])
            idx = np.argmax(reshped_num[i, :])
            seq.append(standard_AA[idx])
        return "".join(seq)

    def softmax(self, x):
        """
        Applies the softmax function to an array, commonly used in machine learning for classification.

        Args:
            x (numpy.ndarray): The input array.

        Returns:
            numpy.ndarray: The array after applying the softmax function.
        """
        a = np.exp(x)
        b = np.sum(np.exp(x))
        return a / b

    def optimize(self, constraintName=None, constraint_lb_tol_percent=0.9, constraint_ub_tol_percent=1.1, maxiter=1000):
        """
        Performs global optimization on the sequence with optional constraints.

        Args:
            constraintName (str, optional): The name of the constraint to be applied during optimization.
            constraint_lb_tol_percent (float, optional): Lower bound tolerance percentage for the constraint.
            constraint_ub_tol_percent (float, optional): Upper bound tolerance percentage for the constraint.
            maxiter (int, optional): Maximum number of iterations for the optimization process.

        Returns:
            tuple: A tuple containing the optimized sequence and a dictionary of results including optimized sequence, antifungal index, and any applied constraint.

        Raises:
            ValueError: If an invalid constraint is specified.
        """
        if constraintName is None:
            # nlc = None
            nlc = ()

        elif constraintName in ['MW', 'Charge', 'ChargeDensity', 'pI', 'InstabilityInd', 'Aromaticity', 'AliphaticInd',
                                'BomanInd', 'HydrophRatio']:
            lb = self.constraint(self.transfer_alphabeta_num(self.sequence), constraintName) * constraint_lb_tol_percent
            ub = self.constraint(self.transfer_alphabeta_num(self.sequence), constraintName) * constraint_ub_tol_percent
            nlc = NonlinearConstraint(lambda seq_num: self.constraint(seq_num, constraintName), lb, ub)

        else:
            raise ValueError(
                "The constraint is not valid, please choose from ['MW', 'Charge', 'ChargeDensity', 'pI', 'InstabilityInd', 'Aromaticity', 'AliphaticInd', 'BomanInd', 'HydrophRatio']")

        # minimize
        minimization = minimize(self.cost,
                                x0=self.transfer_alphabeta_num(self.sequence),
                                method='COBYLA',
                                constraints=nlc,
                                options={'maxiter': maxiter, 'disp': False})
        optimized_seq = self.trasfer_num_alphabeta(minimization.x)

        if constraintName is None:
            dict_result = {'optimized_seq': optimized_seq,
                           'antifungal_idx': self.cost(minimization.x)}
        elif constraintName in ['MW', 'Charge', 'ChargeDensity', 'pI', 'InstabilityInd', 'Aromaticity', 'AliphaticInd',
                                'BomanInd', 'HydrophRatio']:
            dict_result = {'optimized_seq': optimized_seq,
                           'antifungal_idx': self.cost(minimization.x),
                           'constraintName': constraintName,
                           'constraint': self.constraint(minimization.x, constraintName)}

        return optimized_seq, dict_result
