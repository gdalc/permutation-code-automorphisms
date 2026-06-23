def sample_perm_with_order(k):
    '''
    Samples permutation if length n = 2*k and with order k. This is done by sampling
    two disjoint cycles of length k. The cycles are joined to create a string, then the string is used to create a permutation
    of length n = 2*k
    '''
    cycles = '('
    P = Permutations(2*k).random_element()

    for i in range(k-1):
        cycles+=str(P[i])+','
    cycles += str(P[k-1])+')('

    for i in range(k-1):
            cycles+=str(P[k+i])+','
    cycles +=str(P[2*k-1])+')'

    return cycles, Permutation(cycles).to_matrix()


def gen_code(Fq, n, k):
    '''
    Generate random code whose automorphism group contains P (and its powers).
    P is sampled as a permutation of length n and order k
    This is done by picking a random codeword c of length n, and then we compute
    c*P^i, for each i = 0,1, ..., k-1
    This way, we obtain the k rows of the generator matrix
    We also check if the code rate is ok
    '''
    ok = 0
    while ok == 0:

        cycles, P = sample_perm_with_order(k)

        G = matrix(Fq, k, n)
        a = random_matrix(Fq, 1, n)
        G[0,:] = a
        for i in range(1,k):
            a = a*P
            G[i,:] = a
        if rank(G) == k:
              ok = 1

    #change of basis
    S = random_matrix(Fq,k,k)
    while rank(S)<k:
        S = random_matrix(Fq,k,k)

    return cycles, P, S*G


def decrement_cycles(s):
    '''
    Decrements by 1 indices in input cycles decomposition
    '''
    return re.sub(r'\d+', lambda m: str(int(m.group()) - 1), s)
