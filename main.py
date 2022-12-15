import sys
import csv


# Create cell class
class Cell:
    """
    Cell class, represents a cell in the body
    :param name: name of the cell (str)
    :param genome: genome of the cell (list(str))
    :param reading_frame: reading frame of the cell (list(int))
    """

    def __init__(self, name, genome, reading_frame):

        self.name = name
        self.genome = genome
        self.reading_frame = reading_frame

    def __str__(self):
        """
        :print: string representation of the cell and its genome
        """
        print(f'<{self.name}, {self.genome}>')

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
        list_of_subs = self.__find_all_subs(dna_sequence)
        # create list of tuples containing substring and amount of repeat
        list_of_repeat = []
        # send each substring to the max_repeating function
        for sub in list_of_subs:
            # value is the max value of repeated element
            value = self.__max_repeating(dna_sequence, sub)
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
        return [(self.find_srr(i), self.translate(i)) for i in range(len(self.genome))]


class StemCell(Cell):
    """
    StemCell class
    """

    def __init__(self, name, genome, reading_frame):
        """
        :param genome: list of strings
        :param reading_frame: list of ints
        """
        super().__init__(name, genome, reading_frame)

    # multiply function that given an int, returns a list of deepcopy of self without external libraries
    def __mul__(self, n):
        """
        :param n: int
        :return: list of deepcopy of self
        """
        # create list that will contain copies of self where the first element is self
        copies = [self]
        # iterate n-1 times
        for i in range(n - 1):
            # append copy of self to copies
            copies.append(type(self)(**vars(self)))
        return copies

    # mitosis function that returns a list of two stem cells one that is self and the other is a deepcopy of self
    def mitosis(self):
        """
        :return: list of two stem cells one that is self and the other is a deepcopy of self
        """
        return [self, type(self)(**vars(self))]

    # differentiate function that returns a new cell of type cell_type and parameters args
    def differentiate(self, cell_type, args):
        if cell_type == "NC":
            return NerveCell(self, args)
        elif cell_type == "MC":
            return MuscleCell(self, *args)


class NerveCell(Cell):
    """
    NerveCell class
    """

    def __init__(self, origin_cell: StemCell, args):
        super().__init__("NerveCell", origin_cell.genome, origin_cell.reading_frame)
        self.signal = None
        self.coef = args

    def receive(self, strength):
        """
        :param strength: int
        :return: int
        """
        self.signal = strength * self.coef

    def send(self):
        """
        :return: int
        """
        return self.signal


class MuscleCell(Cell):
    """
    MuscleCell class
    """

    def __init__(self, origin_cell: StemCell, args):
        super().__init__("MuscleCell", origin_cell.genome, origin_cell.reading_frame)
        list_args = args.split(',')
        self.file = list_args[0]
        self.threshold = list_args[1]

    def receive(self, strength):
        if strength > self.threshold:
            # write to file
            with open(self.file, "w") as f:
                f.write(str(strength) + "I like to move it\n")


class NerveNetwork:

    def __init__(self, nerve_list, musclecell):
        self.nerve_list = nerve_list
        self.muscle = musclecell

    def send_signal(self, strength):
        # iterate over nerve_list and send strength to each nerve cell
        for nerve in self.nerve_list:
            nerve.receive(strength)
            strength = nerve.send()

        self.muscle.receive(strength)


if __name__ == '__main__':
    # receive argv from command line where argv[1] is the file name and argv[2] is list of strength
    argv = sys.argv
    list_signal = argv[2].split(",")
    list_of_nerve = []
    muscle = None

    # Open the TSB file
    with open(argv[1], 'r') as tsb_file:
        # Create a DictReader object
        tsb_reader = csv.DictReader(tsb_file, delimiter='\t')

        # Iterate over the rows in the TSB file
        for row in tsb_reader:
            # Create a new cell
            stemy = StemCell("StemCell", row['DNA'].split(), [int(i) for i in row['reading_frames'].split(',')])
            if row['type'] == 'NC':
                nerve1, nerve2 = stemy.mitosis()
                nerve1.differentiate("NC", row['parameter'])
                nerve2.differentiate("NC", row['parameter'])
                list_of_nerve.append(nerve1)
                list_of_nerve.append(nerve2)
            else:
                muscle = stemy.differentiate("MC", row['parameter'])

    # create a nerve network
    nerve_network = NerveNetwork(list_of_nerve, muscle)
    for signal in list_signal:
        nerve_network.send_signal(int(signal))
        print(muscle.repertoire())




