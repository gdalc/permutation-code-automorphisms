load('../utils.sage')
import time
import csv
import re

#proof.all(False)

q = 3 #finite field size
num_test = 100 #number of tests for each value of k
do_proof = False #set to True if you want to test each of the found automorphisms

max_attempts_per_k = 20000 #safety cap on sampling attempts for each k
progress_every = 500 #print a progress line every this many attempts

#test_type = "OUR" #test only our algorithm
#test_type = "BUILTIN" #test only builtin
test_type = "BOTH" #test both algorithms

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



def find_b(Fq, n):
    '''
    Closed-form (non-bruteforce) choice of b in Fq^* such that
    M = I + b*J is nondegenerate, i.e. det(M) = 1 + n*b != 0 (Remark 3.6, a=1).
    Only the single degenerate value of b is excluded - no rank/determinant
    search over candidate b's.
    '''
    n_Fq = Fq(n)

    if n_Fq == 0:
        return Fq(1)

    forbidden = -Fq(1)/n_Fq

    for b in Fq:
        if b != 0 and b != forbidden:
            return b

    return None


#For each value of k, sample a fully random [n,k]_q code (n = 2k, matching the
#length scale used elsewhere). Keep only codes with dim(Hull(C)) = 1 and
#(x,1) != 0 (i.e. GI-reducible per Theorem 4.5).

for k in range(3,100):

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

            #Recover the hull vector x = c*G, where c spans the left kernel of G*G^T
            ker = GGt.left_kernel()
            c = ker.basis()[0]
            x = c*G

            #Theorem 4.5: C is GI-reducible iff (x,1) != 0. This eligibility
            #check (like the h==1 check above) is treated as instance
            #selection - common "test harness" overhead needed to pick a
            #comparable set of codes for OUR vs BUILTIN - and is therefore
            #NOT included in the timed reduction below, consistent with how
            #the h==0 filter is untimed in test_trivial_hull.sage.
            if sum(x) == 0:
                continue

            num_codes += 1

            if (test_type == "OUR")|(test_type == "BOTH"):
                #Our reduction - timing now starts BEFORE finding a,b, since
                #for a random code this search is a genuine part of running
                #the algorithm (Algorithm 1, Steps 1-2), not something we can
                #precompute from known structure as in test_hull_one.sage.
                start = time.time()

                b = find_b(Fq, n) #closed-form choice of b (a = 1 fixed)
                J = matrix(Fq, n, n, lambda i,j: 1)
                M = identity_matrix(Fq, n) + b*J

                GMGt = G*M*G.transpose() #guaranteed nonsingular by Theorem 4.5
                A = M*G.transpose()*(GMGt)^-1*G #generalized orthogonal projector Pi_{C,M}

                #Pi_{C,M} is generally not symmetric (unless M = I), so it is
                #the weighted adjacency matrix of a directed graph (Section 6)
                graph_A = DiGraph(A.change_ring(ZZ), weighted=True)
                auto_group = graph_A.automorphism_group(edge_labels=True) #solver for automorphisms of digraph
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
