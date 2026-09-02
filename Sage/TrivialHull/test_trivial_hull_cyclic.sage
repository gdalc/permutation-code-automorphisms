load('../utils.sage')
import time
import csv
import re

#proof.all(False)

q = 3 #finite field size
num_test = 100 #number of tests for each value of k
do_proof = False #set to True if you want to test each of the found automorphisms

min_k = 3
max_k = 100

max_attempts_per_k = 20000 #safety cap on sampling attempts for each k

test_type = "OUR" #test only our algorithm
#test_type = "BUILTIN" #test only builtin
#test_type = "BOTH" #test both algorithms

print("Test for q = ", q)
print("num_test = ", num_test)
print("do_proof = ",do_proof)

print("k from ",min_k, " to ", max_k-1)

print("max_attempts_per_k = ",max_attempts_per_k)
print("test_type = ",test_type)
print("-------------------\n")

#Prepare files for output
output_dir = "Results"
os.makedirs(output_dir, exist_ok=True)

if (test_type == "OUR")|(test_type == "BOTH"):
    file_name_OUR = "Results/q"+str(q)+"_h0_cyclic_OUR.csv"
    #Prepare csv file
    with open(file_name_OUR, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(["k", "Our"])

if (test_type == "BUILTIN")|(test_type == "BOTH"):
    file_name_BUILTIN = "Results/q"+str(q)+"_h0_cyclic_BUILTIN.csv"
    #Prepare csv file
    with open(file_name_BUILTIN, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(["k", "BuiltIn"])


Fq = GF(q)

#For each value of k, sample code with automorphism group of size at least k.
#Then, check if it has trivial hull and, when this happens, compare Sage-built in function
#with our reduction to the Automorphisms Graph Problem

for k in range(min_k, max_k):

    n = 2*k

    our = 0 #average time for our solver
    sage_built_in = 0 #average time for built-in function

    num_codes = 0 #number of tested codes

    data = []

    while num_codes < num_test:

        #Generate code
        cycles, P, G = gen_code(q, n, k)
        h = k-rank(G*G.transpose()) #hull dimension

        errors = 0 #set to 1 if something bad happens

        if h == 0:
            num_codes += 1

            if (test_type == "OUR")|(test_type == "BOTH"):
                #Our reduction
                start = time.time()

                auto_group = find_code_aut_trivial_hull(G)

                end = time.time()

                our += end-start

                #Check validity
                if do_proof:
                    initial_auto_group = PermutationGroup([decrement_cycles(cycles)], domain = [i for i in range(n)])

                    #Check if each of the initial automorphisms is included in the found group
                    for A in initial_auto_group:
                        if A not in auto_group:
                            errors = 1

                    #Check if each automorphism is indeed an automorphism of the code
                    H = codes.LinearCode(G).parity_check_matrix()
                    for A in auto_group:
                        A_matrix = A.matrix()
                        s = G*P*H.transpose()
                        for i in range(k):
                            for j in range(k):
                                if s[i,j]:
                                    errors = 1

                    #Check if some of found automorphisms were not ok
                    if errors:
                        print("############## SOMETHING WENT WRONG")

            if (test_type == "BUILTIN")|(test_type == "BOTH"):
                #Now, Sage built-in function
                C = codes.LinearCode(G)

                start = time.time()
    #            gens = C.automorphism_group_gens()
                gens = C.permutation_automorphism_group()
                end = time.time()

                sage_built_in += end-start

    #Print data
    print("Avg times: k = ",k)
    if (test_type == "OUR")|(test_type == "BOTH"):
        print("- Our = ",N(our/num_test))

        data = [k, N(our/num_test)]
        with open(file_name_OUR, 'a', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(data)

    if (test_type == "BUILTIN")|(test_type == "BOTH"):
        print("- Built-in = ", N(sage_built_in/num_test))
        data = [k, N(sage_built_in/num_test)]
        with open(file_name_BUILTIN, 'a', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(data)
    print("-------------------")
