# issquare
# bit_digits, issquare
def bit_digits(n):
    tmp = n
    cnt = 0
    while tmp:
        tmp >>= 1
        cnt += 1
    return cnt
def issquare(n: int):
    l, r = 0, n + 1
    d = (l + r) // 2
    cnt = bit_digits(n+1)
    for _ in range(cnt + 1):
        if d * d <= n:
            l = d
            d = (l + r) // 2
        else:
            r = d
            d = (l + r) // 2
    return d

# integer operation
# gcd, factorization
def gcd(a, b):
    k, l = a, b
    while l:
        k, l = l, k % l
    return k
def factorization(n):
    D = {}
    p = 2
    a = n
    while a != 1:
        cnt = 0
        while a % p == 0:
            cnt += 1
            a //= p
        if cnt:
            D[p] = cnt
        p += 1
        if p * p > n and a != 1:
            D[a] = 1
            break
    return D

# fraction operation
# frac_sum, frac_prod, fracted
def frac_sum(f1, f2, sign=1): # f1[0] / f1[1] + sign * f2[0] / f2[1]
    f = [f1[0]*f2[1]+sign*f1[1]*f2[0], f1[1]*f2[1]]
    g = gcd(f[0], f[1])
    f[0] //= g
    f[1] //= g
    return f
def frac_prod(f1, f2): # (f1[0] / f1[1]) * (f2[0] / f2[1])
    f = [f1[0]*f2[0], f1[1]*f2[1]]
    g = gcd(f[0], f[1])
    f[0] //= g
    f[1] //= g
    return f
def fracted(f):
    g = gcd(f[0], f[1])
    a = f[0] // g
    b = f[1] // g
    return [a, b]

# sqrt operation
# sqrt_sum, parse_sqrt, sqrt_prod, sqrt_frac, parse_number
# ----- note: 0 is {}, so DO NOT use 0 as sqrt dictionary's key. -----
def sqrt_sum(d1, d2, sign = 1): #d = {s1: c1, s2: c2, ...} => c1*sqrt(s1) + c2*sqrt(s2) + ...
    dres = {i: d1[i][:] for i in d1}
    for i in d2:
        if i in d1:
            dres[i] = frac_sum(dres[i][:], frac_prod([sign, 1], d2[i]))
        else:
            dres[i] = frac_prod([sign, 1], d2[i])
    real_dres = {}
    for i in dres:
        if dres[i][0]:
            real_dres[i] = dres[i][:]
    return real_dres
def parse_sqrt(n: int): # n -> sqrt(n)
    sign = 1
    if n < 0:
        sign = -1
        n = -n
    D = factorization(n)
    coeff = 1
    res = 1
    for i in D:
        coeff *= i ** (D[i] // 2)
        res *= i * (D[i] % 2) + 1 * (1 - D[i] % 2)
    return {sign * res: [coeff, 1]}
def sqrt_prod(d1, d2): # d1 * d2
    dres = {}
    for i in d1:
        for j in d2:
            coeff = frac_prod(d1[i], d2[j])
            res = i * j
            if i < 0 and j < 0:
                coeff[0] *= -1
            light_d = parse_sqrt(res)
            light_L = [[i, light_d[i]] for i in light_d]
            coeff = frac_prod(coeff, light_L[0][1])
            res = light_L[0][0]
            if res in dres:
                dres[res] = frac_sum(dres[res][:], coeff)
            else:
                dres[res] = coeff[:]
    dres2 = {}
    for i in dres:
        if dres[i][0]: dres2[i] = dres[i][:]
    return dres2
def sqrt_frac(d1, d2): # d1 / d2
    upper = {}
    downer = {}
    for i in d1:
        if d1[i][0]: upper[i] = d1[i][:]
    for i in d2:
        if d2[i][0]: downer[i] = d2[i][:]
    if downer == {}:
        raise ZeroDivisionError("division by zero")
    while len(downer) != 1 or 1 not in downer:
        conj_d2 = {}
        left = {}
        right = {}
        k = 0
        for j in downer:
            if j != 1 and k == 0:
                k = j
        for j in downer:
            if j % k == 0:
                right[j] = downer[j][:]
                conj_d2[j] = frac_prod([-1, 1], downer[j][:])
            else:
                left[j] = downer[j][:]
                conj_d2[j] = downer[j][:]
        left = sqrt_prod(left, left)
        right = sqrt_prod(right, right)
        downer = sqrt_sum(left, right, -1)
        upper = sqrt_prod(upper, conj_d2)
    res = {}
    for i in upper:
        res[i] = frac_prod(upper[i][:], [downer[1][1], downer[1][0]])
    return res
def parse_number(d1: dict):
    res = 0
    for i in d1:
        res += d1[i][0] * i**(0.5) / d1[i][1]
    return res

# matrix operation
# kronecker_product, mat_sum, mat_prod, identity_matrix
def kronecker_product(mat1, mat2):
    n1, m1 = len(mat1), len(mat1[0])
    n2, m2 = len(mat2), len(mat2[0])
    res = [[{} for _ in range(m1*m2)] for __ in range(n1*n2)]
    for i in range(n2*n1):
        for j in range(m1*m2):
            res[i][j] = sqrt_prod(mat1[i//n2][j//m2], mat2[i%n2][j%m2])
    return res
def mat_sum(mat1, mat2, sign=1):
    n1, m1 = len(mat1), len(mat1[0])
    n2, m2 = len(mat2), len(mat2[0])
    if n1 != n2 or m1 != m2:
        raise ValueError("matrix size is invalid")
    res = [[{} for _ in range(m1)] for __ in range(n1)]
    for i in range(n1):
        for j in range(m1):
            res[i][j] = sqrt_sum(mat1[i][j], mat2[i][j], sign)
    return res
def mat_prod(mat1, mat2):
    n1, m1 = len(mat1), len(mat1[0])
    n2, m2 = len(mat2), len(mat2[0])
    if m1 != n2:
        raise ValueError("matrix size is invalid")
    res = [[{} for _ in range(m2)] for __ in range(n1)]
    for i in range(n1):
        for j in range(m2):
            S = {}
            for k in range(m1):
                S = sqrt_sum(S, sqrt_prod(mat1[i][k], mat2[k][j]))
            res[i][j] = S
    return res
def identity_matrix(matrix_size):
    res = [[{} for _ in range(matrix_size)] for __ in range(matrix_size)]
    for i in range(matrix_size):
        res[i][i] = {1: [1, 1]}
    return res

# vector operator
def inner_product(vec1, vec2):
    n = len(vec1)
    m = len(vec2)
    res = {}
    if n != m:
        raise ValueError("vector length is invalid")
    for i in range(n):
        res = sqrt_sum(res, sqrt_prod(vec1[i], vec2[i]))
    return res
def narrow_normalization(vec):
    n = len(vec)
    res = [vec[i] for i in range(n)]
    norm = inner_product(vec, vec)
    if norm == {}:
        raise ValueError("zero-vector cannot be normalize")
    else:
        norm = sqrt_frac(parse_sqrt(norm[1][0]), parse_sqrt(norm[1][1]))
        for i in range(n):
            res[i] = sqrt_frac(res[i], norm)
        return res
def scalar_multiplication(c, vec):
    n = len(vec)
    res = [sqrt_prod(vec[i], c) for i in range(n)]
    return res
def vec_sum(vec1, vec2, sign=1):
    n = len(vec1)
    m = len(vec2)
    if n != m:
        raise ValueError("invalid vectors")
    res = [sqrt_sum(vec1[i], vec2[i], sign) for i in range(n)]
    return res
def normalization(vec):
    n = len(vec)
    norm = 0
    for i in range(n):
        norm += vec[i]*vec[i]
    if norm:
        norm = norm**0.5
        return [vec[i] / norm for i in range(n)]
    else:
        raise ValueError("zero-vector cannot be normalize.")
def gramm_schmidt(vecs):
    vec_num = len(vecs)
    dim_flg = True
    for i in range(vec_num-1):
        dim_flg *= (len(vecs[i]) == len(vecs[i+1]))
    if not dim_flg:
        raise ValueError("invalid vectors")
    vec_dim = len(vecs[0])
    res = [[vecs[i][j] for j in range(vec_dim)] for i in range(vec_num)]
    for i in range(vec_num):
        for j in range(i):
            res[i] = vec_sum(res[i], scalar_multiplication(inner_product(vecs[i], res[j]), res[j]), -1)
        res[i] = narrow_normalization(res[i])
    return res

# equation solver
# quadratic_equation_int, quadratic_equation_float
def quadratic_equation_int(a: int, b: int, c: int):
    if a == 0:
        if b == 0:
            if c == 0:
                return "any number"
            return "no solution"
        return [{1: frac_sum([0, 1], [-c, b])}]
    x = parse_sqrt(b*b - 4*a*c)
    y = {}
    if b:
        y[1] = [-b, 1]
    x1 = sqrt_sum(y, x)
    x2 = sqrt_sum(y, x, -1)
    for i in x1:
        x1[i][1] *= 2*a
        x1[i] = frac_prod([1, 1], x1[i])
    for i in x2:
        x2[i][1] *= 2*a
        x2[i] = frac_prod([1, 1], x2[i])
    return [x1, x2]
def quadratic_equation_float(a, b, c):
    if a == 0:
        if b == 0:
            if c == 0:
                return "any number"
            return "no solution"
        return [-c / b]
    x = (b*b - 4*a*c) ** 0.5
    y = -b
    return [(y + x) / (2*a), (y - x) / (2*a)]
def linear_equation(mat, vec): # type == dict; key: "solution", "subsolution"
    n, m = len(mat), len(mat[0])
    nv = len(vec)
    ext_mat = [[{} for _ in range(m+1)] for __ in range(n)]
    if n != nv:
        raise ValueError("invalid equation, and there is no solution")
    for i in range(n):
        for j in range(m):
            ext_mat[i][j] = mat[i][j]
    for i in range(n):
        ext_mat[i][m] = vec[i]
    pnt = 0
    for i in range(m+1):
        row_pnt = pnt
        flg = False
        while row_pnt < n:
            if ext_mat[row_pnt][i] == {}:
                row_pnt += 1
            else:
                flg = True
                break
        if flg:
            for j in range(m+1):
                ext_mat[pnt][j], ext_mat[row_pnt][j] = ext_mat[row_pnt][j], ext_mat[pnt][j]
            keep = ext_mat[pnt][i]
            for j in range(m+1):
                ext_mat[pnt][j] = sqrt_frac(ext_mat[pnt][j], keep)
            for j in range(pnt+1, n):
                raw_keep = ext_mat[j][i]
                for k in range(m+1):
                    ext_mat[j][k] = sqrt_sum(ext_mat[j][k], sqrt_prod(raw_keep, ext_mat[pnt][k]), -1)
            pnt += 1
    flg = True
    for i in range(n):
        fflg = True
        for j in range(m):
            if ext_mat[i][j] != {}:
                fflg = False
        if fflg and ext_mat[i][m] != {}:
            flg = False
    if flg:
        for i in range(n, 0, -1):
            pnt = 0
            fflg = False
            while pnt <= m:
                if ext_mat[i-1][pnt] == {}:
                    pnt += 1
                else:
                    fflg = True
                    break
            if fflg:
                keep = ext_mat[i-1][pnt]
                for j in range(i-1):
                    raw_keep = ext_mat[j][pnt]
                    for k in range(m+1):
                        ext_mat[j][k] = sqrt_sum(ext_mat[j][k], sqrt_prod(raw_keep, ext_mat[i-1][k]), -1)
        chk = [-1 for _ in range(m)]
        D = {}
        cnt = 0
        for i in range(n):
            for j in range(m):
                if ext_mat[i][j] != {}:
                    chk[j] = i
                    D[i] = j
                    cnt += 1
                    break
        solution = {"solution": [{} for _ in range(m)], "subspace": [[{} for _ in range(m)] for __ in range(n-cnt)]}
        pos = 0
        for i in range(m):
            if chk[i] == -1:
                solution["subspace"][pos][i] = {1: [1, 1]}
                for k in D:
                    solution["subspace"][pos][D[k]] = sqrt_prod({1: [-1, 1]}, ext_mat[k][i])
                pos += 1
            else:
                solution["solution"][i] = ext_mat[chk[i]][m]
        for i in range(n-cnt):
            keep = {}
            for j in range(m):
                keep = sqrt_sum(sqrt_prod(solution["subspace"][i][j], solution["subspace"][i][j]), keep)
            keep = sqrt_frac(parse_sqrt(keep[1][0]), parse_sqrt(keep[1][1]))
            for j in range(m):
                solution["subspace"][i][j] = sqrt_frac(solution["subspace"][i][j], keep)
        return solution
    return "no solution"

# angular momentum algebra
# updown_operator, momentum_x, momentum_y, momentum_z, momentum_synthesis
# J_x = (J_+ + J_-) / 2
# J_y = (J_+ - J_-) / 2i
# J^2 = (J_+J_- + J_+J_-) / 2 + J_z^2
def updown_operator(n, op):
    res = [[{} for _ in range(n)] for __ in range(n)]
    if op == "upper":
        for i in range(n-1):
            res[i][i+1] = parse_sqrt((i+1)*(n-i-1))
        return res
    elif op == "lower":
        for i in range(n-1):
            res[i+1][i] = parse_sqrt((i+1)*(n-i-1))
        return res
    else:
        raise ValueError("operator valuable (2nd valuable) is not correct. it must be 'upper' or 'lower'.")
def momentum_x(n):
    upper = updown_operator(n, "upper")
    lower = updown_operator(n, "lower")
    res = mat_sum(upper, lower)
    for i in range(n):
        for j in range(n):
            res[i][j] = sqrt_prod(res[i][j], {1: [1, 2]})
    return res
def momentum_y(n):
    upper = updown_operator(n, "upper")
    lower = updown_operator(n, "lower")
    res = mat_sum(upper, lower, -1)
    for i in range(n):
        for j in range(n):
            res[i][j] = sqrt_prod(res[i][j], {-1: [-1, 2]})
    return res
def momentum_z(n):
    val = fracted([n-1, 2])
    res = [[{} for _ in range(n)] for __ in range(n)]
    for i in range(n):
        if val[0]: res[i][i] = {1: val[:]}
        val = frac_sum(val, [1, 1], -1)
    return res
def momentum_synthesis(n, m): #J_1 + J_2, dim(J_1) = n, dim(J_2) = m
    left_x = kronecker_product(momentum_x(n), identity_matrix(m))
    left_y = kronecker_product(momentum_y(n), identity_matrix(m))
    left_z = kronecker_product(momentum_z(n), identity_matrix(m))
    right_x = kronecker_product(identity_matrix(n), momentum_x(m))
    right_y = kronecker_product(identity_matrix(n), momentum_y(m))
    right_z = kronecker_product(identity_matrix(n), momentum_z(m))
    synth_x = mat_sum(left_x, right_x)
    synth_y = mat_sum(left_y, right_y)
    synth_z = mat_sum(left_z, right_z)
    synth_x = mat_prod(synth_x, synth_x)
    synth_y = mat_prod(synth_y, synth_y)
    synth_z = mat_prod(synth_z, synth_z)
    res = [[{} for __ in range(n*m)] for _ in range(n*m)]
    res = mat_sum(res, synth_x)
    res = mat_sum(res, synth_y)
    res = mat_sum(res, synth_z)
    return res

# test input / output
N = int(input()) # <- number of spins
A = list(map(int, input().split())) # <- sizes of spins' matrix (0 -> 1, 1/2 -> 2, ...)
dim = 1
ms = [[[[{1: [1, 1]}]] for __ in range(3)] for _ in range(N)]
for i in range(N):
    dim *= A[i]
    for j in range(N):
        if i == j:
            ms[j][0] = kronecker_product(ms[j][0], momentum_x(A[i]))
            ms[j][1] = kronecker_product(ms[j][1], momentum_y(A[i]))
            ms[j][2] = kronecker_product(ms[j][2], momentum_z(A[i]))
        else:
            for k in range(3):
                ms[j][k] = kronecker_product(ms[j][k], identity_matrix(A[i]))
msx = [[{} for _ in range(dim)] for __ in range(dim)]
msy = [[{} for _ in range(dim)] for __ in range(dim)]
msz = [[{} for _ in range(dim)] for __ in range(dim)]
for i in range(N):
    msx = mat_sum(msx, ms[i][0])
    msy = mat_sum(msy, ms[i][1])
    msz = mat_sum(msz, ms[i][2])
mssz = [[{} for _ in range(dim)] for __ in range(dim)]
for i in range(dim):
    for j in range(dim):
        mssz[i][j] = msz[i][j]
msx = mat_prod(msx, msx)
msy = mat_prod(msy, msy)
msz = mat_prod(msz, msz)
mss = mat_sum(msx, msy)
mss = mat_sum(mss, msz)
spin_size = sum(A) - N
zeros = [{} for _ in range(dim)]
for i in range(spin_size, -1, -2):
    SS = {}
    if i:
        SS = {1: fracted([i*(i+2), 4])}
    tmps = [[{} for __ in range(dim)] for _ in range(dim)]
    for j in range(dim):
        for k in range(dim):
            tmps[j][k] = mss[j][k]
            if j == k:
                tmps[j][k] = sqrt_sum(tmps[j][k], SS, -1)
    sol = linear_equation(tmps, zeros)
    sstr = ""
    if i % 2:
        sstr = str(i) + "/" + str(2)
    else:
        sstr = str(i//2)
    if sol["subspace"] != {}:
        dictM = {}
        for j in sol["subspace"]:
            vec = [j[k] for k in range(dim)]
            nve = [{} for _ in range(dim)]
            for k in range(dim):
                SSS = {}
                for l in range(dim):
                    SSS = sqrt_sum(SSS, sqrt_prod(vec[l], mssz[k][l]))
                nve[k] = SSS
            MM = {}
            for k in range(dim):
                MM = sqrt_sum(MM, sqrt_prod(vec[k], nve[k]))
            mstr = ""
            if MM == {}:
                mstr = "0"
            else:
                if MM[1][1] == 1:
                    mstr = str(MM[1][0])
                else:
                    mstr = str(MM[1][0]) + "/" + str(MM[1][1])
            if mstr in dictM:
                dictM[mstr][0] += 1
                dictM[mstr].append(j)
            else:
                dictM[mstr] = [1, j[:]]
        for j in dictM:
            z = gramm_schmidt(dictM[j][1:])
            for k in range(dictM[j][0]):
                print("|S = {}, M = {}, cnt = {}> =".format(sstr, j, k+1), end = ' ')
                ccnt = 0
                for l in range(dim):
                    if z[k][l] != {}:
                        tmpk = l
                        tmpt = [0 for _ in range(N)]
                        for m in range(N, 0, -1):
                            tmpt[m-1] = tmpk % A[m-1]
                            tmpk //= A[m-1]
                        if ccnt:
                            print(" +", end = '')
                        print(z[k][l], end = '')
                        for m in range(N):
                            print("|{}>".format(tmpt[m]), end = '')
                        ccnt += 1
                print()
