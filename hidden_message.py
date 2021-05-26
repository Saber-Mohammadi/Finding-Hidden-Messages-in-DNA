import numpy as np
import random

class Genome:

    def __init__(self, sequence):
        self.sequence = sequence
        self.kmer_map = {}
        self.kmer_list = []
        self.CheckValidSeq()

    def GetKmer(self, kmer):
        for i in range(0, len(self.sequence) - kmer + 1):
            pattern = self.sequence[i: i + kmer]
            try:
                if self.kmer_map[pattern]:
                    self.kmer_map[pattern] += 1
            except KeyError:
                self.kmer_map.update({pattern: 1})
        return self.kmer_map

    def DisplaySequence(self):
        return self.sequence

    def CheckValidSeq(self):
        nucleotide = ["A", "T", "C", "G"]
        for i in self.sequence:
            if i in nucleotide:
                continue
            else:
                raise Exception("Not valid sequence")

    def GetKmerList(self, kmer):
        for i in range(0, len(self.sequence) - kmer + 1):
            pattern = self.sequence[i: i + kmer]
            self.kmer_list.append(pattern)
        return list(dict.fromkeys(self.kmer_list))

    # this method calculates the hamming distance between two sequences
    # input: two sequences
    # output: hamming distance
    @staticmethod
    def Hamming_Distance(first_seq, second_seq):
        distance = 0
        if len(first_seq) != len(second_seq):
            raise ValueError("Not Equal Lengths")
        for i in range(len(first_seq)):
            if first_seq[i] != second_seq[i]:
                distance += 1
        return distance

    # Skew function calculates #G-C for finding location of Ori. It takes sequence as input and gives you a numpy.array
    # which indicates the value of G-C in each position. it always starts with 0.
    @staticmethod
    def Skew(genome):
        skew = 0
        output_list = [0]
        for i in genome:
            if i == "C":
                skew -= 1
                output_list.append(skew)
            elif i == "G":
                skew += 1
                output_list.append(skew)
            else:
                output_list.append(skew)
                continue
        return np.array(output_list)

    # this method calculates all the mutants of a pattern with
    # d mismatches
    @staticmethod
    def GetNeighbors(Pattern, d):
        if d == 0:
            return [Pattern]
        if len(Pattern) == 1:
            return ["A", "C", "G", "T"]
        neighborhood = []
        suffix_neighbor = Genome.GetNeighbors(Pattern[1:], d)
        a = Pattern[1:]
        for i in suffix_neighbor:
            if Genome.Hamming_Distance(Pattern[1:], i) < d:
                for j in ["A", "T", "G", "C"]:
                    neighborhood.append(j + i)
            else:
                neighborhood.append(Pattern[0] + i)
        return neighborhood

    @staticmethod
    # this method checks whether a pattern with d mismatch is present in
    # a list of strings
    # input: PatternCheckInListOfString(*["ATTA", "ATCC"], pattern="AGT", d=1)
    # output: [True, False]
    def PatternCheckInListOfStrings(*args, pattern, d):
        pattern_match = []
        for string in args:
            flag = False
            for i in range(len(string) - len(pattern) + 1):
                sub_str = string[i: i + len(pattern)]
                if Genome.Hamming_Distance(sub_str, pattern) <= d:
                    flag = True
                    break
                else:
                    continue
            pattern_match.append(flag)
        return pattern_match

    # Given a collection of strings Dna and an integer d,
    # a k-mer is a (k,d)-motif if it appears in every string
    # from Dna with at most d mismatches
    # this function returns (k,d)-motif
    @staticmethod
    def MotifEnumeration(*args, k, d):
        kmers_map = {}
        neighbors_list = []
        motif_list = []
        for dna in args:
            dna_object = Genome(dna)
            kmer_for_dna = dna_object.GetKmer(k)
            for key in kmer_for_dna.keys():
                if key in kmers_map:
                    kmers_map[key] += 1
                else:
                    kmers_map.update({key: kmer_for_dna[key]})
        for kmer in kmers_map.keys():
            temp_neighbor = Genome.GetNeighbors(kmer, d)
            neighbors_list.extend(temp_neighbor)
        neighbors_list = list(set(neighbors_list))
        for i in neighbors_list:
            status = Genome.PatternCheckInListOfStrings(*args, pattern=i, d=d)
            sum_of_elements = 0
            for s in status:
                sum_of_elements += s
            if sum_of_elements == len(args):
                motif_list.append(i)
        return list(set(motif_list))

    # This method finds a kmer in a larger text that is minimum
    @staticmethod
    def MinHammingDistance(pattern, text):
        hamming_map = {}
        final_map = {}
        for i in range(len(text) - len(pattern) + 1):
            sub_text = text[i: i + len(pattern)]
            dist = Genome.Hamming_Distance(pattern, sub_text)
            hamming_map.update({sub_text: dist})
        min_map = min(hamming_map.values())
        return min_map
        # for key in hamming_map.keys():
        #    if hamming_map[key] == min_map:
        #        final_map.update({key: hamming_map[key]})
        # return min(hamming_map.items(), key=lambda x: x[1])

    @staticmethod
    def MedianString(*args, k):
        kmers_list = []
        for dna in args:
            dna_object = Genome(dna)
            kmers = Genome.GetKmer(dna_object, k)
            kmers_list.extend(list(kmers.keys()))
        kmers_set = list(set(kmers_list))
        dna_min_maps = {}
        for i in kmers_set:
            counter = 0
            for dna in args:
                min_distance = Genome.MinHammingDistance(pattern=i, text=dna)
                counter += min_distance
            dna_min_maps.update({i: counter})
        return min(dna_min_maps.items(), key=lambda x: x[1])

    # This method takes a list of sequences with equal length and returns the profile matrix
    @staticmethod
    def MakeProfile(*args):
        args = list(args)
        adenosine_list = []
        cytosine_list = []
        guanine_list = []
        thymine_list = []
        final_profile = []
        sequence_length = len(args[0])
        for i in range(sequence_length):
            frequency_map = {"A": 0, "T": 0, "C": 0, "G": 0}
            for j in range(len(args)):
                if args[j][i] == "A":
                    frequency_map["A"] += 1
                elif args[j][i] == "C":
                    frequency_map["C"] += 1
                elif args[j][i] == "G":
                    frequency_map["G"] += 1
                elif args[j][i] == "T":
                    frequency_map["T"] += 1
            adenosine_list.append(frequency_map["A"]/len(args))
            cytosine_list.append(frequency_map["C"]/len(args))
            guanine_list.append(frequency_map["G"]/len(args))
            thymine_list.append(frequency_map["T"]/len(args))
        final_profile.append(adenosine_list)
        final_profile.append(cytosine_list)
        final_profile.append(guanine_list)
        final_profile.append(thymine_list)
        return np.reshape(final_profile, newshape=(4, sequence_length)) + 1

    # This method accepts a list of list as profile matrix and returns the probability of text
    # argument occurrence based on profile matrix
    @staticmethod
    def GetProbability(*args, text):
        probability = 1
        array_profile = np.array(args)
        array_profile = array_profile.reshape(4, len(text))
        for i in range(len(text)):
            if text[i] == "A":
                probability *= float(array_profile[0][i])
                t = float(array_profile[0][i])
            elif text[i] == "C":
                probability *= float(array_profile[1][i])
                z = float(array_profile[1][i])
            elif text[i] == "G":
                probability *= float(array_profile[2][i])
                u = float(array_profile[2][i])
            elif text[i] == "T":
                probability *= float(array_profile[3][i])
                q = float(array_profile[3][i])
        return probability

    # This method uses a profile(probability of each nucleotide in each position) to
    # return most probable kmer in a sequence. k=kmer, *args should be given as list of list
    @staticmethod
    def ProfileMostProbableKmer(*args, sequence, k):
        max_probability = 0
        most_probable_kmer = []
        dna_object = Genome(sequence=sequence)
        kmers_list = dna_object.GetKmerList(kmer=k)
        for i in kmers_list:
            probability = Genome.GetProbability(*args, text=i)
            if probability > max_probability:
                max_probability = probability
                most_probable_kmer = i
        return most_probable_kmer

    # This method takes a list of sequences (1*n) and returns consensus sequence and
    # sequences profile
    @staticmethod
    def GetConsensusSequence(*args):
        data = args
        profile = Genome.MakeProfile(*data)
        consensus_sequence = ""
        profile_shape = profile.shape
        for i in range(profile_shape[1]):
            max_column = np.argmax(profile[:, i])
            if max_column == 0:
                consensus_sequence += "A"
            if max_column == 1:
                consensus_sequence += "C"
            if max_column == 2:
                consensus_sequence += "G"
            if max_column == 3:
                consensus_sequence += "T"
        return consensus_sequence

    # This method takes a list of strings as a motif matrix and consensus sequence of this
    # motif matrix. It calculates the hamming distance between each motif and consensus sequence.
    # It returns the sum of these hamming distances.
    # input: (["ACGTA", "ACCTA", "TTGTA"], "ACGTA")
    # output: int(sum of the hamming distances)
    @staticmethod
    def GetScore(*args):
        args = np.array(args)
        args = [str(item[0]) for item in args.tolist()]
        consensus = Genome.GetConsensusSequence(*args)
        score = 0
        for i in args:
            score += Genome.Hamming_Distance(i, consensus)
        return score

    # This method accepts a list of DNA sequences, kmer, the number of DNAs and gives
    # a list of best motifs
    # Hidden Messages in DNA course, week 3
    @staticmethod
    def GreedyMotifSearch(*Dna, KMER, number):
        Dna = list(Dna)
        dna_object_list = []
        consensus_dna = Genome.GetConsensusSequence(*Dna)
        for i in Dna:
            dna_object = Genome(sequence=i)
            dna_object_list.append(dna_object)
        best_motifs = []
        for k_mer in Dna:
            sub_dna = k_mer[:KMER]
            best_motifs.append(sub_dna)
        best_motifs = np.reshape(np.array(best_motifs), newshape=(number, 1))
        first_string_kmers = dna_object_list[0].GetKmerList(kmer=KMER)
        for k in range(len(first_string_kmers)):
            motif_1 = first_string_kmers[k]
            motifs_list_for_profile = [motif_1]
            motifs = [motif_1]
            for i in range(1, number):
                profile_matrix = Genome.MakeProfile(*motifs_list_for_profile)
                motif_i = Genome.ProfileMostProbableKmer(*profile_matrix, sequence=Dna[i], k=KMER)
                motifs_list_for_profile.append(motif_i)
                motifs.append(motif_i)
            motifs = np.reshape(np.array(motifs), newshape=(number, 1))
            score_motifs = Genome.GetScore(*motifs)
            score_best_motifs = Genome.GetScore(*best_motifs)
            # flattened_motif = [str(item[0]) for item in motifs.tolist()]
            # flattened_best_motif = [str(item[0]) for item in best_motifs.tolist()]
            if score_motifs < score_best_motifs:
                best_motifs = motifs
        return best_motifs

    # This method calculates the minimum hamming distance between a pattern
    # and a set of DNAs (should be given as a list to method)
    @staticmethod
    def DistanceBetweenPatternAndString(pattern, *args):
        distance = 0
        for dna in args:
            hamming_distance = 10000
            dna_object = Genome(sequence=dna)
            dna_object_kmers = dna_object.GetKmerList(kmer=len(pattern))
            for kmer in dna_object_kmers:
                kmer_pattern_hamming_dist = Genome.Hamming_Distance(pattern, kmer)
                if kmer_pattern_hamming_dist < hamming_distance:
                    hamming_distance = kmer_pattern_hamming_dist
            distance += hamming_distance
        return distance

    @staticmethod
    def RandomizedMotifSearch(*Dna, KMER):
        best_motifs = []
        for dna in Dna:
            dna_object = Genome(dna)
            dna_object_kmers = dna_object.GetKmerList(kmer=KMER)
            random_kmer_index = random.randint(0, len(dna_object_kmers) - 1)
            random_kmer = dna_object_kmers[random_kmer_index]
            best_motifs.append(random_kmer)
        while True:
            motifs = []
            profile = Genome.MakeProfile(*best_motifs)
            for dna in Dna:
                probable_kmer = Genome.ProfileMostProbableKmer(*profile, sequence=dna, k=KMER)
                motifs.append(probable_kmer)
            best_motifs_score = Genome.GetScore(*best_motifs)
            motifs_score = Genome.GetScore(*motifs)
            if motifs_score < best_motifs_score:
                best_motifs = motifs
            else:
                return best_motifs

