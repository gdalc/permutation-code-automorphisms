load('../utils.sage')
import time
import csv
import re

#proof.all(False)

q = 5 #finite field size
num_test = 100 #number of tests for each value of k
do_proof = False #set to True if you want to test each of the found automorphisms

min_k = 3
max_k = 100

max_attempts_per_k = 20000 #safety cap on sampling attempts for each k
progress_every = 500 #print a progress line every this many attempts

#test_type = "OUR" #test only our algorithm
#test_type = "BUILTIN" #test only builtin
test_type = "BOTH" #test both algorithms

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
    file_name_OUR = "Results/q"+str(q)+"_h1_random_OUR.csv"
    with open(file_name_OUR, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(["k", "Our"])

if (test_type == "BUILTIN")|(test_type == "BOTH"):
    file_name_BUILTIN = "Results/q"+str(q)+"_h1_random_BUILTIN.csv"
    with open(file_name_BUILTIN, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(["k", "BuiltIn"])


Fq = GF(q)



#For each value of k, sample a fully random [n,k]_q code (n = 2k, matching the
#length scale used elsewhere). Keep only codes with dim(Hull(C)) = 1 and
#(x,1) != 0 (i.e. GI-reducible per Theorem 4.5).

for k in range(min_k, max_k):

    n = 2*k

    our = 0 #average time for our solver (includes finding a,b)
    sage_built_in = 0 #average time for built-in function

    num_codes = 0 #number of tested codes
    attempts = 0 #number of sampling attempts (for progress reporting / safety cap)

    data = []

    while (num_codes < num_test) and (attempts < max_attempts_per_k):

        attempts += 1

        #if attempts % progress_every == 0:
        #    #print("  ...k =",k,": attempt",attempts,"/",max_attempts_per_k,
        #          "- found",num_codes,"/",num_test,"valid hull-1 codes so far")

        #Generate a fully random code (no imposed automorphism structure)
        G = gen_random_code(q, n, k)
        GGt = G*G.transpose()
        h = k-rank(GGt) #hull dimension (standard inner product)

        errors = 0 #set to 1 if something bad happens

        if h == 1:

            if is_reducible(G, GGt):
                continue

            num_codes += 1

            if (test_type == "OUR")|(test_type == "BOTH"):
                #Our reduction - timing now includes finding a,b (the search
                #for the structure matrix M), not just the graph solve
                start = time.time()

                auto_group = find_code_aut_onedim_hull(G, q, n)
                
                end = time.time()

                our += end-start

                #Check validity
                if do_proof:
                    #No known automorphism to validate against here (the
                    #code has no imposed structure), so we only check that
                    #each found automorphism is genuinely an automorphism
                    #of the code.
                    H = codes.LinearCode(G).parity_check_matrix()
                    for Aut in auto_group:
                        Aut_matrix = Aut.matrix()
                        s = G*Aut_matrix*H.transpose()
                        for i in range(k):
                            for j in range(H.nrows()):
                                if s[i,j]:
                                    errors = 1

                    if errors:
                        print("############## SOMETHING WENT WRONG")

            if (test_type == "BUILTIN")|(test_type == "BOTH"):
                #Now, Sage built-in function
                C = codes.LinearCode(G)

                start = time.time()
                gens = C.permutation_automorphism_group()
                end = time.time()

                sage_built_in += end-start

    #Print data
    print("Avg times: k = ",k," (found",num_codes,"/",num_test,"in",attempts,"attempts)")

    if num_codes < num_test:
        print("  !! hit max_attempts_per_k before reaching num_test - skipping this k")
        print("-------------------")
        continue

    if (test_type == "OUR")|(test_type == "BOTH"):
        print("- Our = ",N(our/num_codes))

        data = [k, N(our/num_codes)]
        with open(file_name_OUR, 'a', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(data)

    if (test_type == "BUILTIN")|(test_type == "BOTH"):
        print("- Built-in = ", N(sage_built_in/num_codes))
        data = [k, N(sage_built_in/num_codes)]
        with open(file_name_BUILTIN, 'a', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(data)
    print("-------------------")
