# 1 -> 2nd elt, 2 -> 3rd, 4 -> 4th, etc
def subset_binary(I, ss):
    first_set = []
    second_set = []
    mask = 1
    first_pass = True
    for x in I:
        if first_pass:
            first_set.append(x)
            first_pass = False
        else:
            if mask&ss!=0:
                second_set.append(x)
            else:
                first_set.append(x)
            mask *= 2
    return first_set, second_set

# enumerate the partitions of a set into non-empty subsets
# partitions are not ordered
def set_partitions(init_set):
    n = len(init_set)
    if n==0:
        return [[]]
    res = [[init_set]]
    for ss in xrange(1, 2^(n-1)):
        first_set, second_set = subset_binary(init_set, ss)
        for A in set_partitions(second_set):
            res.append([first_set]+A)
    return res

def char_Lpacket_SO_odd(n):
    res = []
    a = int((n+1)/2)
    b = int(n/2)
    for i in xrange(b+1):
        for set1 in [x[0] for x in set_partitions_2(a, i, range(a))]:
            for set2 in [y[0] for y in set_partitions_2(b, i, range(b))]:
                A = [Mod(0, 2)]*n
                for z in set1:
                    A[2*z] = Mod(1, 2)
                for z in set2:
                    A[2*z+1] = Mod(1, 2)
                res.append(A)
    return res

def char_Lpacket_SO_even(n):
    res = []
    a = Integer(n/2)
    for i in xrange(a+1):
        for set1 in [x[0] for x in set_partitions_2(a, i, range(a))]:
            for set2 in [y[0] for y in set_partitions_2(a, i, range(a))]:
                A = [Mod(0, 2)]*n
                for z in set1:
                    A[2*z] = Mod(1, 2)
                for z in set2:
                    A[2*z+1] = Mod(1, 2)
                res.append(A)
    return res

def char_Lpacket_Sp(n):
    if n==0:
        return [[]]
    res = []
    for x in char_Lpacket_Sp(n-1):
        res.append([Mod(0,2)]+x)
        res.append([Mod(1,2)]+x)
    return res

def formal_archi_Arthur_param_BC(n):
    if n==0:
        return [[[],0]]
    res = [[[], n]]
    for d in xrange(1, n+1):
        for x in formal_archi_Arthur_param_BC(n-d):
            res.append([[d]+x[0], x[1]])
    return res

def mask_of_formal_archi_Arthur_param_BC(P):
    res = 0
    cur_m = 1
    for i in xrange(P[1]):
        res |= cur_m
        cur_m *= 2
    for d in P[0][::-1]:
        cur_m *= 2
        for i in xrange(d-1):
            res |= cur_m
            cur_m *= 2
    return res

# Input: formal archimedean Arthur parameter P, discrete series given as its tempered character temp_char
# Output: sign in Johnson's resolution and character on non-tempered S_Psi
def char_and_sign_in_formal_archi_Arthur_param_SO_odd(P, temp_char):
    cur_pos = 0
    non_temp_char = []
    sign_Johnson = Mod(0, 2)
    for d in P[0]:
        a = len([i for i in xrange(d) if Mod(temp_char[cur_pos+i]+i, 2)==Mod(0, 2)])
        sign_Johnson += a*(d-a)
        non_temp_char.append(sum(temp_char[cur_pos:cur_pos+d]))
        cur_pos += d
    d = P[1]
    a = 2*len([i for i in xrange(d) if Mod(temp_char[cur_pos+i]+i, 2)==Mod(0, 2)])
    if Mod(d, 2)==Mod(0, 2):
        a += 1
    sign_Johnson += Mod((a*(2*d+1-a))/2, 2)
    if d>0:
        non_temp_char.append(sum(temp_char[cur_pos:cur_pos+d]))
    return tuple(non_temp_char), sign_Johnson

def char_and_sign_in_formal_archi_Arthur_param_SO_even(P, temp_char):
    cur_pos = 0
    non_temp_char = []
    sign_Johnson = Mod(0, 2)
    for d in P[0]:
        a = len([i for i in xrange(d) if Mod(temp_char[cur_pos+i]+i, 2)==Mod(0, 2)])
        sign_Johnson += a*(d-a)
        non_temp_char.append(sum(temp_char[cur_pos:cur_pos+d]))
        cur_pos += d
    d = P[1]
    if d>0:
        non_temp_char.append(sum(temp_char[cur_pos:cur_pos+d]))
    return tuple(non_temp_char), sign_Johnson

def char_and_sign_in_formal_archi_Arthur_param_Sp(P, temp_char):
    cur_pos = 0
    non_temp_char = []
    sign_Johnson = Mod(P[1]*(P[1]+1)/2, 2)
    for d in P[0]:
        a = len([i for i in xrange(d) if Mod(temp_char[cur_pos+i]+i, 2)==Mod(0, 2)])
        sign_Johnson += a*(d-a)
        non_temp_char.append(sum(temp_char[cur_pos:cur_pos+d]))
        cur_pos += d
    return tuple(non_temp_char), sign_Johnson

def formal_archi_Arthur_param_with_char_SO_odd(n):
    discrete_series = char_Lpacket_SO_odd(n)
    res = []
    for P in formal_archi_Arthur_param_BC(n):
        char_with_mult = {}
        for ds in discrete_series:
            non_temp_char, sign_Johnson = char_and_sign_in_formal_archi_Arthur_param_SO_odd(P, ds)
            if not (non_temp_char in char_with_mult):
                char_with_mult[non_temp_char] = (-1)^sign_Johnson
            else:
                char_with_mult[non_temp_char] += (-1)^sign_Johnson
        res.append([P, mask_of_formal_archi_Arthur_param_BC(P), char_with_mult])
    return res

def formal_archi_Arthur_param_with_char_SO_odd_cpct(n):
    discrete_series = [[Mod(i+n, 2) for i in xrange(n)]]
    res = []
    for P in formal_archi_Arthur_param_BC(n):
        char_with_mult = {}
        for ds in discrete_series:
            non_temp_char, sign_Johnson = char_and_sign_in_formal_archi_Arthur_param_SO_odd(P, ds)
            if not (non_temp_char in char_with_mult):
                char_with_mult[non_temp_char] = (-1)^sign_Johnson
            else:
                char_with_mult[non_temp_char] += (-1)^sign_Johnson
        res.append([P, mask_of_formal_archi_Arthur_param_BC(P), char_with_mult])
    return res

def formal_archi_Arthur_param_with_char_SO_even(n, split=True):
    if split:
        discrete_series = char_Lpacket_SO_even(n)
    else: # holomorphic: signature (2n-2,2) (n odd) or (2n,0) (n even)
        discrete_series = [[Mod(i,2) for i in xrange(n)]]
    res = []
    for P in formal_archi_Arthur_param_BC(n):
        char_with_mult = {}
        for ds in discrete_series:
            non_temp_char, sign_Johnson = char_and_sign_in_formal_archi_Arthur_param_SO_even(P, ds)
            if not (non_temp_char in char_with_mult):
                char_with_mult[non_temp_char] = (-1)^sign_Johnson
            else:
                char_with_mult[non_temp_char] += (-1)^sign_Johnson
        res.append([P, mask_of_formal_archi_Arthur_param_BC(P), char_with_mult])
    return res

def formal_archi_Arthur_param_with_char_Sp(n):
    discrete_series = char_Lpacket_Sp(n)
    res = []
    for P in formal_archi_Arthur_param_BC(n):
        char_with_mult = {}
        for ds in discrete_series:
            non_temp_char, sign_Johnson = char_and_sign_in_formal_archi_Arthur_param_Sp(P, ds)
            if not (non_temp_char in char_with_mult):
                char_with_mult[non_temp_char] = (-1)^sign_Johnson
            else:
                char_with_mult[non_temp_char] += (-1)^sign_Johnson
        res.append([P, mask_of_formal_archi_Arthur_param_BC(P), char_with_mult])
    return res

def formal_archi_Arthur_param_with_char_Siegel(n):
    discrete_series = [[Mod(i,2) for i in xrange(n)]] # hol DS
    res = []
    for P in formal_archi_Arthur_param_BC(n):
        if P[1]==0: # holomorphic discrete series is part of the A-packet
            char_with_mult = {}
            for ds in discrete_series:
                non_temp_char, sign_Johnson = char_and_sign_in_formal_archi_Arthur_param_Sp(P, ds)
                if not (non_temp_char in char_with_mult):
                    char_with_mult[non_temp_char] = (-1)^sign_Johnson
                else:
                    char_with_mult[non_temp_char] += (-1)^sign_Johnson
            res.append([P, mask_of_formal_archi_Arthur_param_BC(P), char_with_mult])
    return res

def formal_global_Arthur_param_with_char_BC(n, group_type):
    if group_type=='B': # SO_odd
        archi_param = formal_archi_Arthur_param_with_char_SO_odd(n)
        parity = 0
    elif group_type=='B_cpct':
        archi_param = formal_archi_Arthur_param_with_char_SO_odd_cpct(n)
        parity = 0
    elif group_type=='C': # Sp
        archi_param = formal_archi_Arthur_param_with_char_Sp(n)
        parity = 1
    elif group_type=='Siegel':
        archi_param = formal_archi_Arthur_param_with_char_Siegel(n)
        parity = 1
    else:
        print "Error in formal_global_Arthur_param_with_char_BC: group_type must be B, C, or Siegel."
        return None
    res = []
    for param, mask, archi_char_dict in archi_param:
        positions = [0]
        cur_pos = 0
        for d in param[0]:
            cur_pos += d
            positions.append(cur_pos)
        d_list = list(set(param[0]))
        num_d = len(d_list)
        d_indices = [[i for i in xrange(len(param[0])) if param[0][i]==d] for d in d_list]
        if not (2*param[1]+parity in d_list):
            for partis in product_lists_gen([set_partitions(A) for A in d_indices]):
                if all([all([Mod(len(x), 2)==Mod(0,2) for x in partis[i]]) for i in xrange(num_d) if Mod(d_list[i]+parity, 2)==Mod(0, 2)]):
                    global_char_dict = {}
                    for archi_char in archi_char_dict:
                        global_char_as_list = []
                        for i in xrange(num_d):
                            for x in partis[i]:
                                global_char_as_list.append(sum([archi_char[y] for y in x]))
                        if parity==0 and param[1]>0:
                            global_char_as_list.append(archi_char[len(param[0])])
                        global_char = tuple(global_char_as_list)
                        if global_char in global_char_dict:
                            global_char_dict[global_char] += archi_char_dict[archi_char]
                        else:
                            global_char_dict[global_char] = archi_char_dict[archi_char]
                    global_param0 = []
                    for i in xrange(num_d):
                        for x in partis[i]:
                            global_param0.append([d_list[i], [positions[y] for y in x]])
                    if parity==1 or param[1]>0 or d_list!=[1] or len(partis[0])>1: # we don't want the tempered non-endoscopic parameter
                        res.append([[global_param0, param[1], []], mask, global_char_dict])
        else:
            spe_d_index = 0
            while d_list[spe_d_index]!=2*param[1]+parity:
                spe_d_index += 1
            for partis in product_lists_gen([set_partitions(A) for A in d_indices]):
                if all([all([Mod(len(x), 2)==Mod(0,2) for x in partis[i]]) for i in xrange(num_d) if i!=spe_d_index and Mod(d_list[i]+parity, 2)==Mod(0, 2)]):
                    must_be_paired = [i for i in xrange(len(partis[spe_d_index])) if Mod(len(partis[spe_d_index][i]),2)==Mod(1,2)]
                    if len(must_be_paired)==1:
                        m = must_be_paired[0]
                        global_char_dict = {}
                        for archi_char in archi_char_dict:
                            global_char_as_list = []
                            for i in xrange(num_d):
                                if i!=spe_d_index:
                                    for x in partis[i]:
                                        global_char_as_list.append(sum([archi_char[y] for y in x]))
                                else:
                                    for j in xrange(len(partis[spe_d_index])):
                                        if j!=m:
                                            global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][j]]))
                            if parity==0:
                                global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][m]])+archi_char[len(param[0])])
                            global_char = tuple(global_char_as_list)
                            if global_char in global_char_dict:
                                global_char_dict[global_char] += archi_char_dict[archi_char]
                            else:
                                global_char_dict[global_char] = archi_char_dict[archi_char]
                        global_param0 = []
                        for i in xrange(num_d):
                            if i!=spe_d_index:
                                for x in partis[i]:
                                    global_param0.append([d_list[i], [positions[y] for y in x]])
                            else:
                                for j in xrange(len(partis[spe_d_index])):
                                    if j!=m:
                                        global_param0.append([d_list[spe_d_index], [positions[y] for y in partis[spe_d_index][j]]])
                        if parity==0 or param[1]>0 or num_d>1 or len(partis[0])>1:
                            res.append([[global_param0, param[1], [positions[y] for y in partis[spe_d_index][m]]], mask, global_char_dict])
                    elif len(must_be_paired)==0:
                        # first we do not pair the central block with any other one
                        global_char_dict = {}
                        for archi_char in archi_char_dict:
                            global_char_as_list = []
                            for i in xrange(num_d):
                                for x in partis[i]:
                                    global_char_as_list.append(sum([archi_char[y] for y in x]))
                            if parity==0 and param[1]>0:
                                global_char_as_list.append(archi_char[len(param[0])])
                            global_char = tuple(global_char_as_list)
                            if global_char in global_char_dict:
                                global_char_dict[global_char] += archi_char_dict[archi_char]
                            else:
                                global_char_dict[global_char] = archi_char_dict[archi_char]
                        global_param0 = []
                        for i in xrange(num_d):
                            for x in partis[i]:
                                global_param0.append([d_list[i], [positions[y] for y in x]])
                        res.append([[global_param0, param[1], []], mask, global_char_dict])
                        # next we pair it
                        if parity==0 or param[1]>0 or num_d>1 or len(partis[0])>1:
                            for m in xrange(len(partis[spe_d_index])):
                                global_char_dict = {}
                                for archi_char in archi_char_dict:
                                    global_char_as_list = []
                                    for i in xrange(num_d):
                                        if i!=spe_d_index:
                                            for x in partis[i]:
                                                global_char_as_list.append(sum([archi_char[y] for y in x]))
                                        else:
                                            for j in xrange(len(partis[spe_d_index])):
                                                if j!=m:
                                                    global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][j]]))
                                    if parity==0:
                                        global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][m]])+archi_char[len(param[0])])
                                    global_char = tuple(global_char_as_list)
                                    if global_char in global_char_dict:
                                        global_char_dict[global_char] += archi_char_dict[archi_char]
                                    else:
                                        global_char_dict[global_char] = archi_char_dict[archi_char]
                                global_param0 = []
                                for i in xrange(num_d):
                                    if i!=spe_d_index:
                                        for x in partis[i]:
                                            global_param0.append([d_list[i], [positions[y] for y in x]])
                                    else:
                                        for j in xrange(len(partis[spe_d_index])):
                                            if j!=m:
                                                global_param0.append([d_list[spe_d_index], [positions[y] for y in partis[spe_d_index][j]]])
                                res.append([[global_param0, param[1], [positions[y] for y in partis[spe_d_index][m]]], mask, global_char_dict])
    return [[P for P in res if P[1]&i==P[1]] for i in xrange(2^n)]


def formal_global_Arthur_param_with_char_SO_even(n, split=True):
    archi_param = formal_archi_Arthur_param_with_char_SO_even(n, split)
    res = []
    for param, mask, archi_char_dict in archi_param:
        positions = [0]
        cur_pos = 0
        for d in param[0]:
            cur_pos += d
            positions.append(cur_pos)
        d_list = list(set(param[0]))
        num_d = len(d_list)
        d_indices = [[i for i in xrange(len(param[0])) if param[0][i]==d] for d in d_list]
        if param[1]==0: # no odd orth part
            for partis in product_lists_gen([set_partitions(A) for A in d_indices]):
                if all([all([Mod(len(x), 2)==Mod(0,2) for x in partis[i]]) for i in xrange(num_d) if Mod(d_list[i], 2)==Mod(1, 2)]):
                    global_char_dict = {}
                    for archi_char in archi_char_dict:
                        global_char_as_list = []
                        for i in xrange(num_d):
                            for x in partis[i]:
                                global_char_as_list.append(sum([archi_char[y] for y in x]))
                        if param[1]>0:
                            global_char_as_list.append(archi_char[len(param[0])])
                        global_char = tuple(global_char_as_list)
                        if global_char in global_char_dict:
                            global_char_dict[global_char] += archi_char_dict[archi_char]
                        else:
                            global_char_dict[global_char] = archi_char_dict[archi_char]
                    global_param0 = []
                    for i in xrange(num_d):
                        for x in partis[i]:
                            global_param0.append([d_list[i], [positions[y] for y in x]])
                    if d_list!=[1] or len(partis[0])>1: # we don't want the tempered non-endoscopic parameter
                        res.append([[global_param0, param[1]], mask, global_char_dict])
        else: # two odd orth
            spe_d = 2*param[1]-1 # odd >=1
            for partis in product_lists_gen([set_partitions(A) for A in d_indices]):
                if all([all([Mod(len(x), 2)==Mod(0,2) for x in partis[i]]) for i, d in enumerate(d_list) if (Mod(d, 2)==Mod(1, 2) and d!=1 and d!=spe_d)]):
                    if not (1 in d_list):
                        first_pairings = [(-1,0)]
                    else:
                        one_index = 0
                        while d_list[one_index]!=1:
                            one_index += 1
                        must_be_paired = [m for m in xrange(len(partis[one_index])) if Mod(len(partis[one_index][m]), 2)==Mod(1,2)]
                        if len(must_be_paired)==1:
                            m = must_be_paired[0]
                            first_pairings = [(one_index, m)]
                        elif len(must_be_paired)==0 and spe_d>1:
                            first_pairings = [(-1,0)]
                            for j, x in enumerate(partis[one_index]):
                                first_pairings.append((one_index, j))
                        else:
                            first_pairings = []
                    if first_pairings: # no point of going on otherwise
                        if spe_d==1:
                            second_pairings = [(-1,0)]
                            if not (1 in d_list):
                                print "Error!"
                            for j, x in enumerate(partis[one_index]):
                                if Mod(len(x), 2)==Mod(0, 2):
                                    second_pairings.append((one_index, j))
                        else:
                            if not (spe_d in d_list):
                                second_pairings = [(-1,0)]
                            else:
                                spe_d_index = 0
                                while d_list[spe_d_index]!=spe_d:
                                    spe_d_index += 1
                                must_be_paired = [m for m in xrange(len(partis[spe_d_index])) if Mod(len(partis[spe_d_index][m]), 2)==Mod(1,2)]
                                if len(must_be_paired)==1:
                                    m = must_be_paired[0]
                                    second_pairings = [(spe_d_index, m)]
                                elif len(must_be_paired)==0:
                                    second_pairings = [(-1,0)]
                                    for j, x in enumerate(partis[spe_d_index]):
                                        second_pairings.append((spe_d_index, j))
                                else:
                                    second_pairings = []
                        for pair1 in first_pairings:
                            for pair2 in second_pairings:
                                global_char_dict = {}
                                for archi_char in archi_char_dict:
                                    global_char_as_list = []
                                    for i in xrange(num_d):
                                        for j, x in enumerate(partis[i]):
                                            if (i,j)!=pair1 and (i,j)!=pair2:
                                                global_char_as_list.append(sum([archi_char[y] for y in x]))
                                #    global_char_odd_part = archi_char[len(param[0])]
                                #    if pair1[0]>=0:
                                #        i, j = pair1
                                #        global_char_odd_part += sum([archi_char[k] for k in partis[i][j]])
                                #    if pair2[0]>=0:
                                #        i, j = pair2
                                #        global_char_odd_part += sum([archi_char[k] for k in partis[i][j]])
                                #    if sum(global_char_as_list)!=global_char_odd_part:
                                #        print "Error"
                                    global_char = tuple(global_char_as_list)
                                    if global_char in global_char_dict:
                                        global_char_dict[global_char] += archi_char_dict[archi_char]
                                    else:
                                        global_char_dict[global_char] = archi_char_dict[archi_char]
                                global_param0 = []
                                for i in xrange(num_d):
                                    for j, x in enumerate(partis[i]):
                                        if (i,j)!=pair1 and (i,j)!=pair2:
                                            global_param0.append([d_list[i], [positions[y] for y in x]])
                                if pair1[0]<0:
                                    pair1_pos = []
                                else:
                                    i, j = pair1
                                    pair1_pos = [positions[y] for y in partis[i][j]]
                                if pair2[0]<0:
                                    pair2_pos = []
                                else:
                                    i, j = pair2
                                    pair2_pos = [positions[y] for y in partis[i][j]]
                                res.append([[global_param0, param[1], [pair1_pos, pair2_pos]], mask, global_char_dict])
    return [[P for P in res if P[1]&i==P[1]] for i in xrange(2^n)]


def epsilon_archi_product(H1, H2):
    res = Mod(0, 2)
    for w1 in H1:
        for w2 in H2:
            res += Mod(1+2*max(w1, w2),2)
    return res

def endo_contrib_SO_odd(n, dom_weight, formal_param, new_sympl_para_list, new_orth_odd_para_list, new_orth_even_para_list):
    res = 0
    for P, Pmask, char_dict in formal_param:
        sub_Hodge_weights = []
        epsilon_Psi_as_list = []
        for d, pos_list in P[0]:
            sub_H_wt = [dom_weight[pos]+n-pos-d/2 for pos in pos_list]
            sub_Hodge_weights.append(sub_H_wt)
        if P[1]==0: # no (odd orth)[even] in the parameter
            for x in xrange(len(P[0])):
                d = P[0][x][0]
                if Mod(d,2)==Mod(1,2):
                    epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(0, 2) and P[0][y][0]>d]), 2))
                else:
                    epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(1, 2) and P[0][y][0]<d]), 2))
        else: # one component of type (odd orth)[even] in the parameter
            orth_odd_Hodge_weights = [dom_weight[pos]+n-pos-P[1] for pos in P[2]]
            for x in xrange(len(P[0])):
                d = P[0][x][0]
                if Mod(d,2)==Mod(1,2):
                    if d<2*P[1]:
                        epsilon_Psi_as_list.append(epsilon_archi_product(sub_Hodge_weights[x], orth_odd_Hodge_weights)+Mod(sum([w+1/2 for w in sub_Hodge_weights[x]]), 2)+Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(0, 2) and P[0][y][0]>d]), 2))
                    else:
                        epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(0, 2) and P[0][y][0]>d]), 2))
                else:
                    epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(1, 2) and P[0][y][0]<d]), 2))
            epsilon_Psi_as_list.append(Mod(sum(epsilon_Psi_as_list), 2))
        epsilon_Psi = tuple(epsilon_Psi_as_list)
        if epsilon_Psi in char_dict:
            P_contrib = char_dict[epsilon_Psi]
            for x in xrange(len(P[0])):
                d = P[0][x][0]
                rk_x = len(P[0][x][1])
                if Mod(d,2)==Mod(1,2):
                    sub_dom_weight = [sub_Hodge_weights[x][l]-rk_x+1/2+l for l in xrange(rk_x)]
                    P_contrib *= new_sympl_para_list[rk_x-1][index_from_dom_weight(rk_x, sub_dom_weight)]
                else:
                    sub_dom_weight = [sub_Hodge_weights[x][l]-rk_x+l+1 for l in xrange(rk_x)]
                    P_contrib *= new_orth_even_para_list[rk_x/2-1][index_from_dom_weight(rk_x, sub_dom_weight)]
            if P[1]>0:
                rk_x = len(P[2])
                sub_dom_weight = [orth_odd_Hodge_weights[l]-rk_x+l for l in xrange(rk_x)]
                P_contrib *= new_orth_odd_para_list[rk_x][index_from_dom_weight(rk_x, sub_dom_weight)]
            res += P_contrib
    return res

def endo_contrib_SO_even(n, dom_weight, formal_param, new_sympl_para_list, new_orth_odd_para_list, new_orth_even_para_list, comp_later = [{},{},{}]):
    res = 0
    count = 0
    for P, Pmask, char_dict in formal_param:
        count += 1
        if Mod(count,100000)==0:
            print "count=", count, ", found", res, ", unknown:", [len(di) for di in comp_later]
        sub_Hodge_weights = []
        epsilon_Psi_as_list = []
        for d, pos_list in P[0]:
            sub_H_wt = [dom_weight[pos]+n-1-pos-(d-1)/2 for pos in pos_list]
            sub_Hodge_weights.append(sub_H_wt)
        if P[1]>0:
            first_orth_odd_Hodge_weights = [dom_weight[pos]+n-1-pos for pos in P[2][0]]
            second_orth_odd_Hodge_weights = [dom_weight[pos]+n-1-pos-(P[1]-1) for pos in P[2][1]]
        for x in xrange(len(P[0])):
            d = P[0][x][0]
            if Mod(d,2)==Mod(0,2): # sympl[even]
                eps = Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(1, 2) and P[0][y][0]<d]), 2)
                if P[1]>0:
                    eps += epsilon_archi_product(sub_Hodge_weights[x], first_orth_odd_Hodge_weights)
                    if d>2*P[1]-1:
                        eps += epsilon_archi_product(sub_Hodge_weights[x], second_orth_odd_Hodge_weights)
                    else:
                        eps += Mod(sum([w+1/2 for w in sub_Hodge_weights[x]]), 2)
                epsilon_Psi_as_list.append(Mod(eps,2))
            else: # (even orth)[odd]
                epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(0, 2) and P[0][y][0]>d]), 2))
        epsilon_Psi = tuple(epsilon_Psi_as_list)
        if epsilon_Psi in char_dict:
            unknown_term = -1
            P_contrib = char_dict[epsilon_Psi]
            for x in xrange(len(P[0])):
                d = P[0][x][0]
                rk_x = len(P[0][x][1])
                if Mod(d,2)==Mod(0,2):
                    sub_dom_weight = [sub_Hodge_weights[x][l]-rk_x+1/2+l for l in xrange(rk_x)]
                    if rk_x<=max_computed_rk_SO_odd:
                        P_contrib *= new_sympl_para_list[rk_x-1][index_from_dom_weight(rk_x, sub_dom_weight)]
                    else:
                        unknown_term = 0
                        unknown_dom_weight = tuple(sub_dom_weight)
                else:
                    sub_dom_weight = [sub_Hodge_weights[x][l]-rk_x+l+1 for l in xrange(rk_x)]
                    if rk_x<=max_computed_rk_SO_even:
                        P_contrib *= new_orth_even_para_list[rk_x/2-1][index_from_dom_weight(rk_x, sub_dom_weight)]
                    else:
                        unknown_term = 2
                        unknown_dom_weight = tuple(sub_dom_weight)
            if P[1]>0:
                rk_x = len(P[2][0])
                if rk_x>0:
                    sub_dom_weight = [first_orth_odd_Hodge_weights[l]-rk_x+l for l in xrange(rk_x)]
                    if rk_x<=max_computed_rk_Sp:
                        P_contrib *= new_orth_odd_para_list[rk_x][index_from_dom_weight(rk_x, sub_dom_weight)]
                    else:
                        unknown_term = 1
                        unknown_dom_weight = tuple(sub_dom_weight)
                rk_x = len(P[2][1])
                if rk_x>0:
                    sub_dom_weight = [second_orth_odd_Hodge_weights[l]-rk_x+l for l in xrange(rk_x)]
                    if rk_x<=max_computed_rk_Sp:
                        P_contrib *= new_orth_odd_para_list[rk_x][index_from_dom_weight(rk_x, sub_dom_weight)]
                    else:
                        unknown_term = 1
                        unknown_dom_weight = tuple(sub_dom_weight)
            if dom_weight[n-1]==0 and P[1]==0:
                P_contrib *= 2
            if P_contrib!=0:
                if unknown_term==-1:
                    res += P_contrib
                else:
                    if unknown_dom_weight in comp_later[unknown_term]:
                        comp_later[unknown_term][unknown_dom_weight] += P_contrib
                    else:
                        comp_later[unknown_term][unknown_dom_weight] = P_contrib
    return res

def endo_to_str_Sp(P, sub_Hodge_weights, orth_odd_Hodge_weights):
    res = ""
    for x in xrange(len(P[0])):
        d = P[0][x][0]
        if Mod(d,2)==Mod(0,2):
            res += "S"+str(tuple(sub_Hodge_weights[x]))+"["+str(d)+"] + "
        else:
            res += "Oe"+str(tuple(sub_Hodge_weights[x]))
            if d>1:
                res += "["+str(d)+"]"
            res += " + "
    if len(P[2])>0:
        res += "Oo"+str(tuple(orth_odd_Hodge_weights))
    if P[1]>0:
        res += "["+str(2*P[1]+1)+"]"
    elif len(P[2])==0:
        res += "1"
    return res

def endo_contrib_Sp(n, dom_weight, formal_param, new_sympl_para_list, new_orth_odd_para_list, new_orth_even_para_list, verbose=False):
    res = 0
    for P, Pmask, char_dict in formal_param:
        sub_Hodge_weights = []
        epsilon_Psi_as_list = []
        for d, pos_list in P[0]:
            sub_H_wt = [dom_weight[pos]+n-pos-(d-1)/2 for pos in pos_list]
            sub_Hodge_weights.append(sub_H_wt)
        orth_odd_Hodge_weights = [dom_weight[pos]+n-pos-P[1] for pos in P[2]]
        for x in xrange(len(P[0])):
            d = P[0][x][0]
            if Mod(d,2)==Mod(0,2): # sympl[even]
                if d>2*P[1]+1:
                    epsilon_Psi_as_list.append(epsilon_archi_product(sub_Hodge_weights[x], orth_odd_Hodge_weights)+Mod(sum([w+1/2 for w in sub_Hodge_weights[x]]), 2)+Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(1, 2) and P[0][y][0]<d]), 2))
                else:
                    epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(1, 2) and P[0][y][0]<d]), 2))
            else: # (even orth)[odd]
                epsilon_Psi_as_list.append(Mod(sum([epsilon_archi_product(sub_Hodge_weights[x], sub_Hodge_weights[y]) for y in xrange(len(P[0])) if Mod(P[0][y][0],2)==Mod(0, 2) and P[0][y][0]>d]), 2))
        epsilon_Psi = tuple(epsilon_Psi_as_list)
        if epsilon_Psi in char_dict:
            P_contrib = char_dict[epsilon_Psi]
            for x in xrange(len(P[0])):
                d = P[0][x][0]
                rk_x = len(P[0][x][1])
                if Mod(d,2)==Mod(0,2):
                    sub_dom_weight = [sub_Hodge_weights[x][l]-rk_x+1/2+l for l in xrange(rk_x)]
                    P_contrib *= new_sympl_para_list[rk_x-1][index_from_dom_weight(rk_x, sub_dom_weight)]
                else:
                    sub_dom_weight = [sub_Hodge_weights[x][l]-rk_x+l+1 for l in xrange(rk_x)]
                    P_contrib *= new_orth_even_para_list[rk_x/2-1][index_from_dom_weight(rk_x, sub_dom_weight)]
            rk_x = len(P[2])
            if rk_x>0:
                sub_dom_weight = [orth_odd_Hodge_weights[l]-rk_x+l for l in xrange(rk_x)]
                P_contrib *= new_orth_odd_para_list[rk_x][index_from_dom_weight(rk_x, sub_dom_weight)]
            if verbose and P_contrib!=0:#all([ddd==dom_weight[0] for ddd in dom_weight]):
                print "  "+endo_to_str_Sp(P, sub_Hodge_weights, orth_odd_Hodge_weights)+": "+str(P_contrib)
            res += P_contrib
    return res



# TODO: remove
#def formal_global_Arthur_param_with_char_SO_odd(n):
#    archi_param = formal_archi_Arthur_param_with_char_SO_odd(n)
#    res = []
#    for param, mask, archi_char_dict in archi_param:
#        positions = [0]
#        cur_pos = 0
#        for d in param[0]:
#            cur_pos += d
#            positions.append(cur_pos)
#        d_list = list(set(param[0]))
#        num_d = len(d_list)
#        d_indices = [[i for i in xrange(len(param[0])) if param[0][i]==d] for d in d_list]
#        if not (2*param[1] in d_list):
#            for partis in product_lists([set_partitions(A) for A in d_indices]):
#                if all([all([Mod(len(x), 2)==Mod(0,2) for x in partis[i]]) for i in xrange(num_d) if Mod(d_list[i], 2)==Mod(0, 2)]):
#                    global_char_dict = {}
#                    for archi_char in archi_char_dict:
#                        global_char_as_list = []
#                        for i in xrange(num_d):
#                            for x in partis[i]:
#                                global_char_as_list.append(sum([archi_char[y] for y in x]))
#                        if param[1]>0:
#                            global_char_as_list.append(archi_char[len(param[0])])
#                        global_char = tuple(global_char_as_list)
#                        if global_char in global_char_dict:
#                            global_char_dict[global_char] += archi_char_dict[archi_char]
#                        else:
#                            global_char_dict[global_char] = archi_char_dict[archi_char]
#                    global_param0 = []
#                    for i in xrange(num_d):
#                        for x in partis[i]:
#                            global_param0.append([d_list[i], [positions[y] for y in x]])
#                    if param[1]>0 or d_list!=[1] or len(partis[0])>1:
#                        res.append([[global_param0, param[1], []], mask, global_char_dict])
#        else:
#            spe_d_index = 0
#            while d_list[spe_d_index]!=2*param[1]:
#                spe_d_index += 1
#            for partis in product_lists([set_partitions(A) for A in d_indices]):
#                if all([all([Mod(len(x), 2)==Mod(0,2) for x in partis[i]]) for i in xrange(num_d) if i!=spe_d_index and Mod(d_list[i], 2)==Mod(0, 2)]):
#                    must_be_paired = [i for i in xrange(len(partis[spe_d_index])) if Mod(len(partis[spe_d_index][i]),2)==Mod(1,2)]
#                    if len(must_be_paired)==1:
#                        m = must_be_paired[0]
#                        global_char_dict = {}
#                        for archi_char in archi_char_dict:
#                            global_char_as_list = []
#                            for i in xrange(num_d):
#                                if i!=spe_d_index:
#                                    for x in partis[i]:
#                                        global_char_as_list.append(sum([archi_char[y] for y in x]))
#                                else:
#                                    for j in xrange(len(partis[spe_d_index])):
#                                        if j!=m:
#                                            global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][j]]))
#                            global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][m]])+archi_char[len(param[0])])
#                            global_char = tuple(global_char_as_list)
#                            if global_char in global_char_dict:
#                                global_char_dict[global_char] += archi_char_dict[archi_char]
#                            else:
#                                global_char_dict[global_char] = archi_char_dict[archi_char]
#                        global_param0 = []
#                        for i in xrange(num_d):
#                            if i!=spe_d_index:
#                                for x in partis[i]:
#                                    global_param0.append([d_list[i], [positions[y] for y in x]])
#                            else:
#                                for j in xrange(len(partis[spe_d_index])):
#                                    if j!=m:
#                                        global_param0.append([d_list[spe_d_index], [positions[y] for y in partis[spe_d_index][j]]])
#                        res.append([[global_param0, param[1], [positions[y] for y in partis[spe_d_index][m]]], mask, global_char_dict])
#                    elif len(must_be_paired)==0:
#                        # first we do not pair the central block with any other one
#                        global_char_dict = {}
#                        for archi_char in archi_char_dict:
#                            global_char_as_list = []
#                            for i in xrange(num_d):
#                                for x in partis[i]:
#                                    global_char_as_list.append(sum([archi_char[y] for y in x]))
#                            if param[1]>0:
#                                global_char_as_list.append(archi_char[len(param[0])])
#                            global_char = tuple(global_char_as_list)
#                            if global_char in global_char_dict:
#                                global_char_dict[global_char] += archi_char_dict[archi_char]
#                            else:
#                                global_char_dict[global_char] = archi_char_dict[archi_char]
#                        global_param0 = []
#                        for i in xrange(num_d):
#                            for x in partis[i]:
#                                global_param0.append([d_list[i], [positions[y] for y in x]])
#                        res.append([[global_param0, param[1], []], mask, global_char_dict])
#                        # next we pair it
#                        for m in xrange(len(partis[spe_d_index])):
#                            global_char_dict = {}
#                            for archi_char in archi_char_dict:
#                                global_char_as_list = []
#                                for i in xrange(num_d):
#                                    if i!=spe_d_index:
#                                        for x in partis[i]:
#                                            global_char_as_list.append(sum([archi_char[y] for y in x]))
#                                    else:
#                                        for j in xrange(len(partis[spe_d_index])):
#                                            if j!=m:
#                                                global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][j]]))
#                                global_char_as_list.append(sum([archi_char[y] for y in partis[spe_d_index][m]])+archi_char[len(param[0])])
#                                global_char = tuple(global_char_as_list)
#                                if global_char in global_char_dict:
#                                    global_char_dict[global_char] += archi_char_dict[archi_char]
#                                else:
#                                    global_char_dict[global_char] = archi_char_dict[archi_char]
#                            global_param0 = []
#                            for i in xrange(num_d):
#                                if i!=spe_d_index:
#                                    for x in partis[i]:
#                                        global_param0.append([d_list[i], [positions[y] for y in x]])
#                                else:
#                                    for j in xrange(len(partis[spe_d_index])):
#                                        if j!=m:
#                                            global_param0.append([d_list[spe_d_index], [positions[y] for y in partis[spe_d_index][j]]])
#                            res.append([[global_param0, param[1], [positions[y] for y in partis[spe_d_index][m]]], mask, global_char_dict])
#    return [[P for P in res if P[1]&i==P[1]] for i in xrange(2^n)]

# TODO: remove
#def formal_endo_param_SO_odd(n):
#    res = []
#    char_list = char_Lpacket_SO_odd(n)
#    for f_endo in set_partitions(range(n))[1:]: # [1:] to remove the non-endoscopic parameter
#        num_char = 0
#        for c in char_list:
#            if all([sum([c[i] for i in pa])==Mod(0, 2) for pa in f_endo]):
#                num_char += 1
#        res.append([f_endo, num_char])
#    return res

