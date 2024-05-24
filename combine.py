import sys
import os


total = len(os.listdir(path=sys.argv[1]))


with open(sys.argv[2], "w") as fp:
    itr = 0
    for f in os.listdir(path=sys.argv[1]):
        print(f"{itr}/{total} | Reading {f}", end="\r")
        with open(f"{sys.argv[1]}/{f}", "r") as src:
            for line in src.readlines():
                fp.write(line)
        fp.write("\n")
        itr += 1
    print()
