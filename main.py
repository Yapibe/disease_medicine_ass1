import sys


def find_srr(sequence):
    lookup = []
    n = len(sequence)
    for i in range(1, 6):
        for j in range(i, n+1):
            if sequence[j-i:j] not in lookup:
                lookup.append(sequence[j-i:j])
    print(lookup)
