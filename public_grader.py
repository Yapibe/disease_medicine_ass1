import sys
import os
import subprocess


def run_on_input(stud_sol, input):
    """ will return the student output as list(output split by  """
    input = input.strip().replace("  ", " ").split(" ")
    command = ["python3", stud_sol]
    command.extend(input)
    p = subprocess.run(command, capture_output=True, text=True)
    return p.stdout.strip().replace("  ", " ").replace("\n\n", "\n").split("\n"), p.stderr


def check_srr(deserved_srrs, out):
    srrs_d = deserved_srrs.split(";")
    srrs_out = out.split(";")
    mistakes = 0
    if not all([s in srrs_out for s in srrs_d]) or not all([s in srrs_d for s in srrs_out]):
        print("WRONG srr list !!! should be: {}, but got: {}".format(srrs_d, out))
        mistakes = +1
    srr_numbers = [int(s.split(",")[-1]) for s in srrs_out]
    if not all([srr_numbers[i] <= srr_numbers[i + 1] for i in range(len(srr_numbers) - 1)]):
        print("WRONG srr order !!!")
        mistakes = +1
    return mistakes


def check_output(outs):
    mistakes = 0
    deserved_outs = [
        "Translation: M;P;P;L;S",
        "No simple repeats in DNA sequence",
        "Non-coding RNA",
        "ATG,3",
        "Translation: M;H;H;H",
    ]
    deserved_srrs = "A,3;GA,3;ATCAA,3;AG,4;G,5"
    mistakes += check_srr(deserved_srrs, outs[0])
    for des, out in zip(deserved_outs, outs[1:]):
        if des != out:
            mistakes += 1
            print("WRONG OUTPUT!!! should be: {}, but got: {}".format(des, out))
    return mistakes


def check_file(file="output.txt"):
    if not os.path.exists(file):
        print("ERROR!!! no output file")
        exit()
    file_deservesed = ["50.80319999999999, I like to move it\n", "76.20479999999999, I like to move it\n"]
    mistakes = 0
    with open(file, "r") as fd:
        file_lines = fd.readlines()
        if len(file_deservesed) != len(file_lines):
            print(
                "ERROR!!! file should have {} lines, but have {} lines".format(
                    len(file_deservesed), len(file_lines)
                )
            )
            exit()
        for des, out in zip(file_deservesed, file_lines):
            if des != out:
                mistakes += 1
                print("WRONG OUTPUT!!! should be {}, but got: {}".format(des, out))
    return mistakes


if __name__ == "__main__":
    solution = sys.argv[1]
    if not os.path.exists(solution):
        print("ERROR!!! pythom solution file {} , not found".format(solution))
        exit(0)
    input = "input.txt 30,200,300"
    print("input is: ", input)
    outs, err = run_on_input(solution, input)
    if err:
        print("ERROR!!! your program crasing with err{}".format(err))
        exit(0)
    print("your output is: ", outs)
    mistakes = check_output(outs)
    mistakes += check_file()
    print("finished  checking with {} mistakes".format(mistakes))
