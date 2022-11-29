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

    # Driver Code


find_srr(sys.argv[1])
