import numpy as np
import math
from fractions import Fraction as fr

def gomory_cut_algo():
    input_file = open(r"test_case_6.txt", "r")
    data = input_file.readlines()

    boo = 1 if (data[1].strip().lower() == "maximize") else 0
    data[4] = data[4].strip().split(",")
    no_of_vars = len(data[4])
    for j in range(no_of_vars):
        data[4][j] = fr(data[4][j].strip())
    a = np.array(data[4], dtype = fr).reshape((1, -1))

    i = 5
    while (data[i] != '\n'):
        data[i] = data[i].strip().split(",")
        for j in range(no_of_vars):
            data[i][j] = fr(data[i][j].strip())
        temp = np.array(data[i], dtype = fr).reshape((1, -1))
        a = np.concatenate((a, temp), axis = 0, dtype = fr)
        i += 1
    m, n = a.shape
    i += 2

    b = np.full((m, 1), fr("0"), dtype = fr)
    for j in range(m):
        b[j, 0] = fr(data[i].strip())
        i += 1
    i += 2

    aug_lst = []
    cnt = 0
    for j in range(m):
        if data[i].strip() == ">=":
            for new_i in range(no_of_vars):
                a[j, new_i] = a[j, new_i] * fr("-1")
            b[j, 0] *= fr("-1")
            aug_lst.append(cnt)
            cnt += 1
        elif data[i].strip() == "<=":
            aug_lst.append(cnt)
            cnt += 1
        else:
            aug_lst.append(-1)
        i += 1
    aug_mat = np.full((m, cnt), fr("0"), dtype = fr)
    for j in range(m):
        if aug_lst[j] != -1:
            aug_mat[j, aug_lst[j]] = fr("1")
    a = np.concatenate((a, aug_mat), axis = 1, dtype = fr)
    i += 2

    for j in range(m):
        if b[j, 0] < fr("0"):
            for new_i in range(no_of_vars+cnt):
                a[j, new_i] = a[j, new_i] * fr("-1")
            b[j, 0] *= fr("-1")

    data[i] = data[i].strip().split(",")
    for j in range(no_of_vars):
        data[i][j] = fr(data[i][j].strip())
    c = np.array(data[i], dtype = fr).reshape((-1, 1))
    if boo:
        for new_i in range(no_of_vars):
            c[new_i, 0] = c[new_i, 0] * fr("-1")
    c = np.concatenate((c, np.full((cnt, 1), fr("0"), dtype = fr)), axis = 0, dtype = fr)
    input_file.close()

    answer_dict = {"number_of_cuts" : 0, "initial_solution" : [0]*no_of_vars, "solution_status" : None, "final_solution" : [0]*no_of_vars, "optimal_value" : 0}

    def simplex(tabl, basic_var):
        m, n=  tabl.shape
        m -= 1
        n -= 1
        ans = {"solution_status": None, "optimal_tableau": tabl, "basic_variables": basic_var}
        while True:
            bool_var = 0
            for i in range(1, n+1):
                if tabl[0, i] < fr("0"):
                    pivot_col = i
                    bool_var = 1
                    break
            if bool_var == 0:
                break
            pivot_row = None
            theta = None
            for i in range(1, m+1):
                if tabl[i, pivot_col] > fr("0"):
                    if (theta is None) or (tabl[i, 0]/tabl[i, pivot_col] < theta):
                        theta = tabl[i, 0]/tabl[i, pivot_col]
                        pivot_row = i
            if pivot_row is None:
                ans["solution_status"] = "unbounded"
                return ans
            tabl[pivot_row] /= tabl[pivot_row, pivot_col]
            for i in range(m+1):
                if i == pivot_row:
                    continue
                tabl[i, :] -= (tabl[pivot_row, :]*tabl[i, pivot_col])
            ans["basic_variables"][pivot_row - 1] = pivot_col - 1
        ans["solution_status"] = "optimal"
        ans["optimal_tableau"] = tabl
        return ans
    
    def dual_simplex(tabl, basic_var):
        m, n=  tabl.shape
        m -= 1
        n -= 1
        ans = {"solution_status": None, "optimal_tableau": tabl, "basic_variables": basic_var}
        while True:
            bool_var = 0
            for i in range(1, m+1):
                if tabl[i, 0] < fr("0"):
                    pivot_row = i
                    bool_var = 1
                    break
            if bool_var == 0:
                break
            pivot_col = None
            theta = None
            for i in range(1, n+1):
                if tabl[pivot_row, i] < fr("0"):
                    if theta is None or (fr("-1") * tabl[0, i]/tabl[pivot_row, i]) < theta:
                        theta = (fr("-1") * tabl[0, i]/tabl[pivot_row, i])
                        pivot_col = i
            if pivot_col is None:
                ans["solution_status"] = "unbounded"
                return ans
            print(theta)
            tabl[pivot_row] /= tabl[pivot_row, pivot_col]
            for i in range(m+1):
                if i == pivot_row:
                    continue
                tabl[i, :] -= (tabl[pivot_row, :]*tabl[i, pivot_col])
            ans["basic_variables"][pivot_row - 1] = pivot_col - 1
        ans["solution_status"] = "optimal"
        ans["optimal_tableau"] = tabl
        return ans
    
    ## create an auxiliary problem for initial bfs
    tableau = np.full((m+1, no_of_vars+1+cnt+m), fr("0"), dtype = fr)
    for i in range(1, m+1):
        for j in range(1, no_of_vars+cnt+1):
            tableau[i, j] = a[i-1, j-1]
        tableau[i, no_of_vars+cnt+i] = fr("1")
        tableau[i, 0] = b[i-1, 0]
    c_aux = np.concatenate((np.full((no_of_vars+cnt, 1), fr("0"), dtype = fr), np.full((m, 1), fr("1"), dtype = fr)), axis = 0, dtype = fr)
    init_aux_sol = np.full((no_of_vars+cnt+m, 1), fr("0"), dtype = fr)
    for j in range(m):
        init_aux_sol[n+cnt+j, 0] = b[j, 0]
    tableau[0, 0] = fr("-1") * (c_aux.T @ init_aux_sol)[0,0]
    tableau[0, 1:] = c_aux.T - np.full((1,m), fr("1"), dtype = fr) @ tableau[1:, 1:]
    initial_basic_vars = list(range(no_of_vars+cnt, no_of_vars+cnt+m, 1))
    ans = simplex(tableau, initial_basic_vars)
    if (ans["solution_status"] == "unbounded") or (ans["solution_status"] == "optimal" and ans["optimal_tableau"][0, 0] < fr("0")):
        answer_dict["solution_status"] = "infeasible"
        return
    redundant = []
    for i in range(len(ans["basic_variables"])):
        if ans["basic_variables"][i] >= no_of_vars+cnt:
            basis_change = None
            for j in range(1, no_of_vars+cnt+1):
                if ans["optimal_tableau"][i+1, j] == fr("0"):
                    continue
                basis_change = j
                break
            if basis_change is None:
                redundant.append(i)
            else:
                for j in range(no_of_vars+cnt+m+1):
                    ans["optimal_tableau"][i+1, j] = ans["optimal_tableau"][i+1, j]/ans["optimal_tableau"][i+1, basis_change]
                for j in range(m+1):
                    if j == i+1:
                        continue
                    for k in range(no_of_vars+cnt+m+1):
                        ans["optimal_tableau"][j, k] = ans["optimal_tableau"][j, k] - (ans["optimal_tableau"][i+1, k]*ans["optimal_tableau"][j, basis_change])
                ans["basic_variables"][i] = j - 1
    for i in range(len(redundant) - 1, -1, -1):
        ans["basic_variables"].pop(redundant[i])
        ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], i+1, axis = 0)
    ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], slice(no_of_vars+cnt+1, no_of_vars+cnt+1+m), axis = 1)

    # recalculate first row of the tableau
    ans["optimal_tableau"][0, 0] = fr("0")
    m, n = ans["optimal_tableau"].shape
    m -= 1
    n -= 1
    for i in range(1, no_of_vars+cnt+1):
        ans["optimal_tableau"][0, i] = c[i-1, 0]
    for i in range(no_of_vars+cnt+1):
        temp = fr("0")
        for j in range(1, m+1):
            temp = temp - (c[ans["basic_variables"][j-1], 0] * ans["optimal_tableau"][j,i])
        ans["optimal_tableau"][0, i] = ans["optimal_tableau"][0, i] + temp
    
    ans = simplex(ans["optimal_tableau"], ans["basic_variables"])
    if (ans["solution_status"] == "unbounded"):
        answer_dict["solution_status"] = "unbounded"
        return answer_dict
    for i in range(len(ans["basic_variables"])):
        if ans["basic_variables"][i] < no_of_vars:
            answer_dict["initial_solution"][ans["basic_variables"][i]] = float(ans["optimal_tableau"][i+1, 0])
    
    ## start gomory cut
    while True:
        temp_lst = []
        for i in range(len(ans["basic_variables"])):
            if ans["basic_variables"][i] >= no_of_vars+cnt:
                temp_lst.append(i)
        for i in range(len(temp_lst)-1 , -1, -1):
            temp_var = ans["basic_variables"][temp_lst[i]]
            ans["basic_variables"].pop(temp_lst[i])
            ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], temp_lst[i]+1, axis = 0)
            temp_lst[i] = temp_var
        temp_lst.sort()
        for i in range(len(temp_lst)-1 , -1, -1):
            ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], temp_lst[i]+1, axis = 1)
        
        m, n = ans["optimal_tableau"].shape
        m -= 1
        n -= 1
        bool_var = 1
        idx = None
        for i in range(1, m+1):
            if ans["optimal_tableau"][i, 0].denominator == 1:
                continue
            else:
                bool_var = 0
                idx = i
                break
        if bool_var == 1:
            break
        ans["optimal_tableau"] = np.concatenate((ans["optimal_tableau"], np.full((m+1, 1), fr("0"), dtype = fr)), axis = 1, dtype = fr)
        n += 1
        ans["optimal_tableau"] = np.concatenate((ans["optimal_tableau"], np.full((1, n+1), fr("0"), dtype = fr)), axis = 0, dtype = fr)
        m += 1
        ans["optimal_tableau"][m, n] = fr("1")
        for i in range(n):
            ans["optimal_tableau"][m, i] = fr(math.floor(ans["optimal_tableau"][idx, i])) - ans["optimal_tableau"][idx, i]
        ans["basic_variables"].append(n)
        answer_dict["number_of_cuts"] += 1
        print(ans["optimal_tableau"][0,0])
        ans = dual_simplex(ans["optimal_tableau"], ans["basic_variables"])
        if ans["solution_status"] == "unbounded":
            answer_dict["solution_status"] = "infeasible"
            return answer_dict
        
    answer_dict["solution_status"] = "optimal"
    for i in range(len(ans["basic_variables"])):
        if ans["basic_variables"][i] < no_of_vars:
            answer_dict["final_solution"][ans["basic_variables"][i]] = float(ans["optimal_tableau"][i+1, 0])
    answer_dict["optimal_value"] = float(ans["optimal_tableau"][0, 0]) if boo == 1 else -float(ans["optimal_tableau"][0, 0])
    return answer_dict

fn_call = gomory_cut_algo()
print("initial_solution: ", end = "")
print(*fn_call["initial_solution"], sep = ", ")
print("final_solution: ", end = "")
print(*fn_call["final_solution"], sep = ", ")
print("solution_status:", fn_call["solution_status"])
print("number_of_cuts:", fn_call["number_of_cuts"])
print("optimal_value:", (fn_call["optimal_value"]))