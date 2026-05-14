import subprocess
import shutil
import os
from pathlib import Path
import re

MARS_JAR = "Mars4_5.jar"
ASM_FILE = "wiener_filter.asm"

TESTCASE_DIR = "Testcases"


def run_mips():

    cmd = [
        "java",
        "-jar",
        MARS_JAR,
        "nc",
        ASM_FILE
    ]

    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )


def read_file(filename):

    with open(filename, "r") as f:
        return f.read().strip()


def cleanup_temp_files():

    temp_files = [
        "input.txt",
        "desired.txt",
        "output.txt"
    ]

    for file in temp_files:
        if os.path.exists(file):
            os.remove(file)


def run_testcase(test_num):

    input_src = f"{TESTCASE_DIR}/input_{test_num}.txt"
    expected_src = f"{TESTCASE_DIR}/expected_{test_num}.txt"
    desired_src = f"{TESTCASE_DIR}/desired.txt"

    saved_output = f"{TESTCASE_DIR}/output_{test_num}.txt"

    shutil.copy(input_src, "input.txt")
    shutil.copy(desired_src, "desired.txt")

    run_mips()

    shutil.copy("output.txt", saved_output)

    actual_output = read_file(saved_output)
    expected_output = read_file(expected_src)

    # Handle error messages
    if expected_output.startswith("Error"):
        passed = actual_output.strip() == expected_output.strip()

    else:
        # Extract all float/int numbers
        actual_nums = list(map(float, re.findall(r'-?\d+\.?\d*', actual_output)))
        expected_nums = list(map(float, re.findall(r'-?\d+\.?\d*', expected_output)))

        if len(actual_nums) != len(expected_nums):
            passed = False
        else:
            passed = all(
                abs(a - e) <= 0.21
                for a, e in zip(actual_nums, expected_nums)
            )

    print(f"\nTest {test_num}: ", end="")

    print(f"\n===== TEST {test_num} =====")

    if passed:
        print("PASSED")
    else:
        print("FAILED")

    print("\nEXPECTED OUTPUT:")
    print(expected_output)

    print("\nACTUAL OUTPUT:")
    print(actual_output)

    cleanup_temp_files()

    return passed


def main():

    total = 0
    passed = 0

    testcase_files = sorted(
        Path(TESTCASE_DIR).glob("input_*.txt")
    )

    for file in testcase_files:

        test_num = file.stem.split("_")[1]

        total += 1

        if run_testcase(test_num):
            passed += 1

    print(f"\nPassed {passed}/{total} tests")


if __name__ == "__main__":
    main()