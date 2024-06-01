#! /usr/bin/env python
import glob
log_file = glob.glob('*_log.txt')
print(log_file)
file1 = open("affinity.out","w")
for i in log_file:
    file2 = open(i,'r')
    data = file2.readlines()
    file2.close()
    for j in data:
        if len(j.split()) >= 2 and j.split( )[0] == "1":
            affinity = i + "   " + j.split()[1]
            print(affinity)
            file1.writelines(affinity + '\n')
file1.close()
