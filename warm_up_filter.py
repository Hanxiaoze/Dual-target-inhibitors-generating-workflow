input_file = './E_M_gen_affinity.txt'
output_file = './E_M_gen_affinity_warm_up.txt'

with open(input_file, 'r') as inf:
    with open(output_file, 'w') as outf:
        for infl in inf.readlines():
            id = int(infl.split('_')[1])
            if id > 100:
                outf.writelines(infl)