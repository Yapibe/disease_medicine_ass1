import sys


def find_srr(dna_sequence):
    """
    function receives a dna sequence (str), and prints a list of tuple
    containing the repeat sequence (str) and amount of repeats (int)
    :param dna_sequence:
    :return:
    """
    # receive list of all sub elements in the string
    list_of_subs = find_all_subs(dna_sequence)
    # create list of tuples containing substring and amount of repeat
    list_of_repeat = []
    # send each substring to the max_repeating function
    for sub in list_of_subs:
        value = max_repeating(dna_sequence, sub)
        if value >= 3:
            list_of_repeat.append((sub, value))
    lst = len(list_of_repeat)
    for i in range(0, lst):
        for j in range(0, lst - i - 1):
            if list_of_repeat[j][1] > list_of_repeat[j + 1][1]:
                temp = list_of_repeat[j]
                list_of_repeat[j] = list_of_repeat[j + 1]
                list_of_repeat[j + 1] = temp
    if not list_of_repeat:
        return None
    return list_of_repeat


def find_all_subs(dna_sequence):
    lookup = []
    n = len(dna_sequence)
    for i in range(1, 6):
        for j in range(i, n + 1):
            if dna_sequence[j - i:j] not in lookup:
                lookup.append(dna_sequence[j - i:j])
    return lookup


def max_repeating(sequence, sub_string):
    # Stores the count of consecutive
    # occurrences of str2 in str1
    count_occurrences = sequence.count(sub_string)
    # Concatenate str2 cntOcc times
    concat = sub_string * count_occurrences
    # Iterate over the string str1
    # while Contstr is not present in str1
    while concat not in sequence:
        # Update count_occurrences
        count_occurrences -= 1
        # Update concat
        concat = sub_string * count_occurrences
    return count_occurrences


def reverse_transcribe(rna_seq):
    dna_seq = []
    for base in rna_seq.upper():
        if base == 'A':
            dna_seq.append('T')
        elif base == 'U':
            dna_seq.append('A')
        elif base == 'G':
            dna_seq.append('C')
        elif base == 'C':
            dna_seq.append('G')
    dna_seq.reverse()
    dna_seq_comp = ""
    for char in dna_seq:
        dna_seq_comp += char
    return dna_seq_comp


def translate(rna_seq, reading_frame):
    dna_seq = []
    for base in rna_seq.upper():
        if base == 'U':
            dna_seq.append('T')
        else:
            dna_seq.append(base)
    rna_seq_t = ""
    for base in dna_seq:
        rna_seq_t += base
    reading = rna_seq_t[reading_frame - 1:]
    if len(reading) % 3 == 1:
        reading = reading[:-1]
    elif len(reading) % 3 == 2:
        reading = reading[:-2]
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein_final = ""
    protein = ""
    met_flag = False

    for i in range(0, len(reading), 3):
        codon = table[reading[i:i + 3]]
        if codon == 'M':
            met_flag = True
        if met_flag:
            if codon == '_':
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


if __name__ == '__main__':

    srr_list = find_srr(sys.argv[1])
    if srr_list:
        for srr in srr_list[:-1]:
            print(f'{srr[0]},{srr[1]}', end=";")
        print(f'{srr_list[-1][0]},{srr_list[-1][1]}')
    else:
        print('No simple repeats in DNA sequence')

    reversed_transcribed = reverse_transcribe(sys.argv[2])
    print("DNA sequence:", reversed_transcribed)

    protein = translate(sys.argv[3], int(sys.argv[4]))
    if protein:
        print("Translation: ", end="")
        for char in protein[:-1]:
            print(char, end=";")
        print(protein[-1:])
    else:
        print("Non-coding RNA")



