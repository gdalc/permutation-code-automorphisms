load('utils.sage')
import time
import csv
import re

#proof.all(False)

q = 3 #finite field size
num_test = 100 #number of tests for each value of k
do_proof = False

file_name = "Results/q"+str(q)+"_h0_cyclic.csv"

#Prepare csv file
with open(file_name, 'w', newline='') as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow(["k", "BuiltIn", "Our"])


Fq = GF(q)

#For each value of k, sample code with automorphism group of size at least k.
#Then, check if it has trivial hull and, when this happens, compare Sage-built in function
#with our reduction to the Automorphisms Graph Problem


for k in range(3,30):

    n = 2*k

    our = 0 #average time for our solver
    sage_built_in = 0 #average time for built-in function

    num = 0 #number of tested codes

    data = []

    while num < num_test:

        #Generate code
        cycles, P, G = gen_code(Fq, n, k)
        h = k-rank(G*G.transpose()) #hull dimension

        errors = 0 #set to 1 if something bad happens

        if h == 0:
            num += 1

            #Our reduction
            start = time.time()

            A = G.transpose()*(G*G.transpose())^-1*G #adjacency matrix
            graph_A = Graph(A.change_ring(ZZ)) #graph
            auto_group = graph_A.automorphism_group() #solver for automorphisms of graph
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

            if errors:
                print("############## SOMETHING WENT WRONG")


            #Now, Sage built-in function
            C = codes.LinearCode(G)

            start = time.time()
#            gens = C.automorphism_group_gens()
            gens = C.permutation_automorphism_group()
            end = time.time()

            sage_built_in += end-start

    print("Avg times: k = ",k,": Our = ",N(our/num_test),", Built-in = ", N(sage_built_in/num_test))

    #Write in csv file
    data = [k, N(sage_built_in/num_test), N(our/num_test)]
    with open(file_name, 'a', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(data)
