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

def gen_random_code(q, n, k):
    '''
    Purely random [n,k]_q code: no imposed automorphism structure at all
    (unlike utils.sage's gen_code, which forces an order-k automorphism).
    '''
    Fq = GF(q)
    G = random_matrix(Fq, k, n)
    while rank(G) < k:
        G = random_matrix(Fq, k, n)
    return G

def gen_code(q, n, k):
    '''
    Generate random code whose automorphism group contains P (and its powers).
    P is sampled as a permutation of length n and order k
    This is done by picking a random codeword c of length n, and then we compute
    c*P^i, for each i = 0,1, ..., k-1
    This way, we obtain the k rows of the generator matrix
    We also check if the code rate is ok
    '''
    Fq = GF(q)
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

def find_b(Fq, n):
    '''
    Closed-form (non-bruteforce) choice of b in Fq^* such that
    M = I + b*J is nondegenerate, i.e. det(M) = 1 + n*b != 0 (Remark 3.6, a=1).

    By Theorem 4.5 / the proof's Schur-complement computation, for a code with
    dim(Hull(C)) = 1 and hull vector x satisfying (x,1) != 0, *every* choice of
    a != 0, b != 0 with a+n*b != 0 makes GMG^T nonsingular. So we do not search
    over candidate b's by testing rank(G*M*G^T) (that would cost O(q) rank
    computations); we only need to avoid the single value of b that makes M
    itself degenerate, which is a single field-element check.
    '''
    n_Fq = Fq(n)

    if n_Fq == 0:
        #characteristic divides n: 1 + n*b = 1 for every b, so any nonzero b works
        return Fq(1)

    forbidden = -Fq(1)/n_Fq #the unique b with 1 + n*b = 0

    for b in Fq:
        if b != 0 and b != forbidden:
            return b

    return None #only happens if q <= 2 and forbidden coincides with the only nonzero element

def is_reducible(G, GGt):
    '''
    Theorem 4.5: C is GI-reducible iff (x,1) != 0.
    '''
    #Recover the hull vector x = c*G, where c spans the left kernel of G*G^T
    #(Hull(C) = {c*G : c*G*G^T = 0}, see proof of Lemma 3.1)
    ker = GGt.left_kernel()
    c = ker.basis()[0]
    x = c*G
    if sum(x) == 0:
        return true
    return False

def find_code_aut_onedim_hull(G, q, n):
    #Reduction
    Fq = GF(q)
    b = find_b(Fq, n)
    J = matrix(Fq, n, n, lambda i,j: 1)
    M = identity_matrix(Fq, n) + b*J

    GMGt = G*M*G.transpose() #guaranteed nonsingular by Theorem 4.5
    A = M*G.transpose()*(GMGt)^-1*G #generalized orthogonal projector Pi_{C,M}

    #Pi_{C,M} is generally not symmetric (unless M = I), so it is
    #the weighted adjacency matrix of a directed graph (Section 6)
    graph_A = DiGraph(A.change_ring(ZZ), weighted=True)
    auto_group = graph_A.automorphism_group(edge_labels=True) #solver for automorphisms of digraph
    return auto_group

def find_code_aut_trivial_hull(G):
    A = G.transpose()*(G*G.transpose())^-1*G #adjacency matrix
    graph_A = Graph(A.change_ring(ZZ)) #graph
    auto_group = graph_A.automorphism_group() #solver for automorphisms of graph
    return auto_group