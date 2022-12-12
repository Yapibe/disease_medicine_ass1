import sys


# if __name__ == '__main__':
#
#     srr_list = find_srr(sys.argv[1])
#     if srr_list:
#         for srr in srr_list[:-1]:
#             print(f'{srr[0]},{srr[1]}', end=";")
#         print(f'{srr_list[-1][0]},{srr_list[-1][1]}')
#     else:
#         print('No simple repeats in DNA sequence')
#
#     reversed_transcribed = reverse_transcribe(sys.argv[2])
#     print("DNA sequence:", reversed_transcribed)
#
#     protein = translate(sys.argv[3], int(sys.argv[4]))
#     if protein:
#         print("Translation: ", end="")
#         for char in protein[:-1]:
#             print(char, end=";")
#         print(protein[-1:])
#     else:
#         print("Non-coding RNA")

# Create cell class
class Cell:
    """
    Cell class represents a cell in a cell division process
    :param name: name of the cell (str)
    :param genome: genome of the cell (list(str))
    :param reading_frame: reading frame of the cell (list(int))
    """

    def __init__(self, name, genome, reading_frame):
        self.name = name
        self.genome = genome
        self.reading_frame = reading_frame

    def __str__(self):
        return f'<{self.name}, {self.genome}>'

    def index_loop(self, index):
        """
        loop index to the beginning of the genome
        :param index: index we want to loop (int)
        :return: looped index (int)
        """
        return index % len(self.genome)

    def find_srr(self, genome_index):
        """
        find simple repeats in the genome of the cell
        :param genome_index: index of the gene we want to search (int)
        :return: list of tuple of simple repeats (list(tuple(str, int)))
        """

        dna_sequence = self.genome[self.index_loop(genome_index)]
        # receive list of all sub elements in the string
        list_of_subs = self.find_all_subs(dna_sequence)
        # create list of tuples containing substring and amount of repeat
        list_of_repeat = []
        # send each substring to the max_repeating function
        for sub in list_of_subs:
            # value is the max value of repeated element
            value = self.max_repeating(dna_sequence, sub)
            # only take value above 3
            if value >= 3:
                list_of_repeat.append((sub, value))
        lst = len(list_of_repeat)
        for i in range(0, lst):
            for j in range(0, lst - i - 1):
                if list_of_repeat[j][1] > list_of_repeat[j + 1][1]:
                    temp = list_of_repeat[j]
                    list_of_repeat[j] = list_of_repeat[j + 1]
                    list_of_repeat[j + 1] = temp
        # if list is empty, meaning no repeat elements
        if not list_of_repeat:
            return None
        # return list of tuple
        return list_of_repeat

    def __find_all_subs(self, dna_sequence):
        """
        find all substrings in the dna sequence
        :param dna_sequence: dna sequence (str)
        :return: list of substrings (list(str))
        """
        lookup = []
        n = len(dna_sequence)
        for i in range(1, 6):
            for j in range(i, n + 1):
                # if this sub element isn't in lookup-add it
                if dna_sequence[j - i:j] not in lookup:
                    lookup.append(dna_sequence[j - i:j])
        return lookup

    def __max_repeating(self, dna_sequence, sub_string):
        """
        returns the max value of each repeating element
        :param dna_sequence: given sequence (str) we wish to find repeats in
        :param sub_string: current sub string (str) we are working on
        :return: amount of consecutive repeats (int) this substring has
        """
        count_occurrences = dna_sequence.count(sub_string)
        # Concatenate sub_string count_occurrences times
        concat = sub_string * count_occurrences
        # Iterate over the string dna_sequences
        # while concat is not present in dna_sequence
        while concat not in dna_sequence:
            # Update count_occurrences
            count_occurrences -= 1
            # Update concat
            concat = sub_string * count_occurrences
        return count_occurrences

    # transcribe function that receives a dna sequence and returns the rna sequence
    def transcribe(self, genome_index):
        """
        :param genome_index:
        :return:
        """
        dna_sequence = self.genome[self.index_loop(genome_index)]
        dic = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        for key in dic:
            dna_sequence = dna_sequence.replace(key, dic[key])
        dna_sequence.reverse()
        return dna_sequence.upper()

    # translate function that receives a genome_index and returns the protein sequence
    def translate(self, genome_index):
        """
        :param genome_index:
        :return:
        """
        index = self.index_loop(genome_index)
        rna_sequence = self.transcribe(index)[self.reading_frame[index] - 1:]
        if len(rna_sequence) % 3 == 1:
            reading = rna_sequence[:-1]
        elif len(rna_sequence) % 3 == 2:
            reading = rna_sequence[:-2]
        # table of amino acids
        table = {
            'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
            'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
            'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
            'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
            'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
            'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
            'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
            'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
            'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
            'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
            'UAA': '_', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
            'UAG': '_', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
            'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
            'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
            'UGA': '_', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
            'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G',
        }
        protein_final = ""
        protein = ""
        # methionine flag, should start read only after M
        met_flag = False
        # read 3 bases every time, 3 bases = codon
        for i in range(0, len(reading), 3):
            codon = table[reading[i:i + 3]]
            # flag change
            if codon == 'M':
                met_flag = True
            if met_flag:
                # address option of stop codon
                if codon == '_':
                    # we want to return longest read
                    if len(protein_final) < len(protein):
                        protein_final = protein
                        protein = ""
                        met_flag = False
                else:
                    protein += codon
        if len(protein_final) < len(protein):
            protein_final = protein
        if not protein_final:
            return None
        else:
            return protein_final

    def repertoire(self):
        # return list of tuples of the form [(find_ssr(genome_index), transcribe(genome_index)]
        return [(self.find_ssr(i), self.transcribe(i)) for i in range(len(self.genome))]


# class NerveCell(Cell):
#     """
#     NerveCell class
#     """
#
#     def __init__(self, genome, reading_frame, StemCell):
#         """
#         :param genome: list of strings
#         :param reading_frame: list of ints
#         """
#         StemCell.__init__(genome, reading_frame)
#         self.coefficient
#         self.StemCell = StemCell
#
#     def __str__(self):
#         return "NerveCell"
#
#     def receive(self, strength):
#         pass
#
#     def send(self, strength):
#         pass


# class StemCell(Cell):
#     """
#     StemCell class
#     """
#
#     def __init__(self, genome, reading_frame):
#         """
#         :param genome: list of strings
#         :param reading_frame: list of ints
#         """
#         super().__init__(genome, reading_frame)
#         self.name = "StemCell"
#
#
#     def __str__(self):
#         return self.name
#
#     def differentiate(self):
#         """
#         :return: NerveCell object
#         """
#         return NerveCell(self.genome, self.reading_frame, self)
