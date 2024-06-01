def my_sort(line):
    line_fields = line.strip().split('   ')
    energy = float(line_fields[1])
    return energy
    
# opening file affinity.out
# and getting contents into a list
fp = open('affinity.out')
contents = fp.readlines()

file1 = open("sorted_final_energy.txt","w")

# sorting using our custom logic
contents.sort(key=my_sort)
# printing the sorting contents to stdout
for line in contents:
    print(line)
    file1.writelines(line)

fp.close()
file1.close()
