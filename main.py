import sys


def find_srr(sequence):
    list_of_subs = find_all_subs(sequence)
    list_of_repeat = []
    for sub in list_of_subs:
        value = max_repeating(sequence, sub)
        if value >= 3:
            list_of_repeat.append((sub, value))

    print(list_of_repeat)


def find_all_subs(sequence):
    lookup = []
    n = len(sequence)
    for i in range(1, 6):
        for j in range(i, n + 1):
            if sequence[j - i:j] not in lookup:
                lookup.append(sequence[j - i:j])
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
    translate(rna_seq, 2)
    dna_seq = []
    for char in rna_seq.upper():
        if char == 'A':
            dna_seq.append('T')
        elif char == 'U':
            dna_seq.append('A')
        elif char == 'G':
            dna_seq.append('C')
        elif char == 'C':
            dna_seq.append('G')
    dna_seq.reverse()
    dna_seq_comp = ""
    for char in dna_seq:
        dna_seq_comp += char
    print(dna_seq_comp)


def translate(rna_seq, reading_frame):
    dna_seq = []
    for char in rna_seq.upper():
        if char == 'U':
            dna_seq.append('T')
        else:
            dna_seq.append(char)
    rna_seq_t = ""
    for char in dna_seq:
        rna_seq_t += char
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
    print(protein_final)


reverse_transcribe(sys.argv[1])
