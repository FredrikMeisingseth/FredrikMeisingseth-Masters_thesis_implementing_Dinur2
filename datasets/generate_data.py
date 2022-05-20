
import numpy as np
from sage.all import *
import sys

def BruteForceSystem(E, n):
    Sols = []

    for x_hat in range(2**n):
        x_hat = list(map(int, list(bin(x_hat)[2:].zfill(n))))
        #print(f"xhat: {x_hat} ")
        
        try:
            if np.all([P_j(*x_hat) == 0 for P_j in E]):
                Sols.append(x_hat)
        except TypeError:
            print(f"# TYPEERROR! #")
            print(f"x_hat: {x_hat}")

    return Sols



# TODO: have a look at whether this truly gives uniformly random polynomials
def GenerateSystem(m, d):
    E = np.array([F2_pols_x.random_element(degree = d, terms = m**4) for j in range(m)])
    return E


def coeffs_to_dict(monomials, coefficients):
    monomials = list(map(tuple, monomials))
    return dict(zip(monomials, coefficients))


def test_solution(E, sol):
    #print(f"# Testing candidate x_hat = {sol}")
    for j,P_j in enumerate(E):
        #print(f"# P_{j} = {P_j(*sol)}")
        if P_j(*sol) == 1:
            return False
    
    return True


def write_to_file(m,n,d, matrix, filename, sols = []):

    with open(filename+".txt", 'w') as file:
        file.write(str(n)+"\n")
        file.write(str(m)+"\n")
        file.write(str(d)+"\n")
        file.write(str(len(matrix[0]))+"\n")
        for j, line in enumerate(matrix):
            line = str(line).replace(",","").replace("[", "").replace("]","")
            file.write(line+"\n")
            # if(j < m-1):
            # else: file.write(line)

    if(len(sols) > 0):
        with open(filename + "_solutions.txt", 'w') as file:
            file.write(str(sols)+"\n")
        
        
    



if __name__ == "__main__":

    if(len(sys.argv)<5):
        raise ValueError("You must give 4 arguments!")

    n = int(sys.argv[1])
    m = int(sys.argv[2])
    d = int(sys.argv[3])
    if(sys.argv[4] == "True"):
        with_sols = True
    elif(sys.argv[4] == "False"):
        with_sols = False
    else:
        raise ValueError("Invalid fourth argument, must be True or False")

    F2_pols_x = BooleanPolynomialRing(n, 'x', order = 'lex')
    x = F2_pols_x.gens()
    E = GenerateSystem(m, d)
    monoms = monomials(x,[2 for i in range(n)])
    monoms.reverse()
    #print(f"monoms: {monoms}")
    full_matrix = [[0 for i in monoms] for j in range(m)]

    E_seq = Sequence(list(E), F2_pols_x)
    coeffs, used_monoms = E_seq.coefficient_matrix(sparse = False)

    #print(f"used monomials:\n {used_monoms}")

    for j in range(m):
        for idx, monomial in enumerate(used_monoms):
            full_matrix[j][monoms.index(monomial[0])] = coeffs[j][idx]
    

    #print(coeffs)
    #print(f"full matrix:\n {full_matrix}")
    #print(str(coeffs[0]).replace(",","").replace("(", "").replace(")",""))

    filename = f"system_{n}_{m}_{d}"
    if(with_sols):
        sols = BruteForceSystem(E, n)
        write_to_file(m,n,d,full_matrix, filename, sols)
    else:
        write_to_file(m,n,d,full_matrix, filename)

    #print(f"sols:\n {sols}")














