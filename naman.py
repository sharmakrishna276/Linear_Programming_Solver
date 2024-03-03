import numpy as np

def simplex_algo():
    input_file = open("input.txt", "r")
    data = input_file.readlines()
    variables = {}

    boo = 1 if (data[1].strip().lower() == "maximize") else 0
    data[4] = data[4][:-1].split()
    no_of_vars = len(data[4])
    for j in range(no_of_vars):
        if data[4][j][-1] == ",":
            data[4][j] = data[4][j][:-1]
        data[4][j] = float(data[4][j])
        data[4].append(-1*data[4][j])
        variables[j+1] = (j, len(data[4])-1)
    a = np.array(data[4]).reshape((1, -1))

    i = 5
    while (data[i] != '\n'):
        data[i] = data[i][:-1].split()
        for j in range(no_of_vars):
            if data[i][j][-1] == ",":
                data[i][j] = data[i][j][:-1]
            data[i][j] = float(data[i][j])
            data[i].append(-1*data[i][j])
        temp = np.array(data[i]).reshape((1, -1))
        a = np.concatenate((a, temp), axis = 0)
        i += 1
    m, n = a.shape
    i += 2

    b = np.zeros(shape=(m, 1))
    for j in range(m):
        b[j, 0] = float(data[i])
        i += 1
    i += 2

    aug_lst = []
    cnt = 0
    for j in range(m):
        if data[i].strip() == ">=":
            a[j] *= -1
            b[j] *= -1
            aug_lst.append(cnt)
            cnt += 1
        elif data[i].strip() == "<=":
            aug_lst.append(cnt)
            cnt += 1
        else:
            aug_lst.append(-1)
        i += 1
    aug_mat = np.zeros(shape=(m, cnt))
    for j in range(m):
        if aug_lst[j] != -1:
            aug_mat[j, aug_lst[j]] = 1
            try:
                variables["slack"].append(j+2*no_of_vars)
            except:
                variables["slack"] = [j+2*no_of_vars]
    a = np.concatenate((a, aug_mat), axis = 1)
    i += 2

    for j in range(m):
        if b[j, 0] < 0:
            a[j] *= -1
            b[j] *= -1

    data[i] = data[i][:-1].split()
    for j in range(no_of_vars):
        if data[i][j][-1] == ",":
            data[i][j] = data[i][j][:-1]
        data[i][j] = float(data[i][j])
        data[i].append(-1*data[i][j])
    c = np.array(data[i]).reshape((-1, 1))
    if boo:
        c *= -1
    c = np.concatenate((c, np.zeros(shape=(cnt, 1))), axis = 0)
    input_file.close()

    for j in range(m):
        try:
            variables["artificial"].append(n+cnt+j)
        except:
            variables["artificial"] = [n+cnt+j]
    a = np.concatenate((a, np.eye(m)), axis = 1)

    answer_dict = {"initial_tableau" : None, "final_tableau" : None, "solution_status" : None, "optimal_solution" : None, "optimal_value" : None}

    c_aux = np.zeros(shape = (n+cnt+m, 1))
    n_, m_ = n+cnt+m, m 
    tableau = np.zeros(shape=(m_+1, n_+1))
    init_aux_sol = np.zeros(shape=(n_, 1))
    for j in range(m_):
        c_aux[n+cnt+j, 0] = 1
        init_aux_sol[n+cnt+j] = b[j]

    tableau[0, 0] = -1 * (c_aux.T @ init_aux_sol)
    tableau[0, 1:n_+1] = c_aux.T - np.ones(shape=(1,m_))@a
    tableau[1:m_+1, 0] = np.reshape(b, (m,))
    tableau[1:m_+1, 1:n_+1] = a

    def simplex(tableau, basic_vars):
        m, n = tableau.shape
        m -= 1
        n -= 1
        answer = {"solution_status" : "", "optimal_tableau" : None, "basic_variables" : None, "not_none" : 0}
        while True:
            for i in range(m+1):
                for j in range(n+1):
                    if abs(tableau[i, j]) < 1e-6:
                        tableau[i, j] = 0
            bool_var = 0
            for j in range(1, n+1):
                if (tableau[0, j]) < 0:
                    pivot_col = j
                    break
            else:
                bool_var = 1
            if bool_var == 1:
                break
            theta = None
            pivot_row = None
            for j in range(1, m+1):
                if tableau[j, pivot_col] > 0:
                    if theta == None or tableau[j, 0]/tableau[j, pivot_col] < theta:
                        theta = tableau[j, 0]/tableau[j, pivot_col]
                        pivot_row = j
            if pivot_row == None:
                answer["solution_status"] = "unbounded"
                answer["optimal_tableau"] = tableau
                answer["not_none"] = 1
                return answer
            tableau[pivot_row] /= tableau[pivot_row, pivot_col]
            for j in range(m+1):
                if j == pivot_row:
                    continue
                tableau[j] -= tableau[pivot_row]*(tableau[j, pivot_col])
            basic_vars[pivot_row-1] = pivot_col-1
        answer["solution_status"] = "optimal"
        answer["basic_variables"] = basic_vars
        answer["optimal_tableau"] = tableau
        answer["not_none"] = 1
        return answer

    basic_vars = list(range(n+cnt, n+cnt+m, 1))
    answer_dict["initial_tableau"] = tableau[1:, :]
    ans = simplex(tableau, basic_vars)
    if ans["solution_status"] == "unbounded" or (ans["solution_status"] == "optimal" and ans["optimal_tableau"][0, 0] + 1e-6 < 0):
        answer_dict["solution_status"] = "infeasible"
    else:
        redundant_constraints = []
        for j in range(m):
            if ans["basic_variables"][j] >= n+cnt:
                for k in range(1, n+cnt+1):
                    if ans["optimal_tableau"][j+1, k] != 0:
                        ans["optimal_tableau"][j+1] /= ans["optimal_tableau"][j+1, k]
                        for z in range(m+1):
                            if z == j+1:
                                continue
                            ans["optimal_tableau"][z] -= ans["optimal_tableau"][j+1]*ans["optimal_tableau"][z, k]
                        ans["basic_variables"][j] = k-1
                        break
                else:
                    redundant_constraints.append(j)
        z = len(redundant_constraints)
        for j in range(z-1, -1, -1):
            ans["basic_variables"].pop(redundant_constraints[j])
            ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], redundant_constraints[j]+1, 0)
        ans["optimal_tableau"] = np.delete(ans["optimal_tableau"], slice(n+cnt,n+cnt+m), axis = 1)

        ans["optimal_tableau"][0] *= 0
        for j in range(ans["optimal_tableau"].shape[0] - 1):
            ans["optimal_tableau"][0] += ans["optimal_tableau"][j+1] * c[ans["basic_variables"][j], 0]
        ans["optimal_tableau"][0] = np.concatenate((np.zeros(shape=(1,1)), c.T), axis = 1) - ans["optimal_tableau"][0]
        
        # answer_dict["initial_tableau"] = ans["optimal_tableau"][1:, :]
        ans = simplex(ans["optimal_tableau"], ans["basic_variables"])

        answer_dict["solution_status"] = ans["solution_status"]
        answer_dict["final_tableau"] = ans["optimal_tableau"][1:, :] if ans["not_none"] == 1 else None
        # answer_dict["final_tableau"] = ans["optimal_tableau"]
        if ans["solution_status"] == "optimal":
            answer_dict["optimal_value"] = ans["optimal_tableau"][0, 0] if boo == 1 else -1 * ans["optimal_tableau"][0, 0]
            answer_dict["optimal_solution"] = [0]*no_of_vars
            ### correct ordering of optimal vector solution
            for j in range(ans["optimal_tableau"].shape[0] - 1):
                if ans["basic_variables"][j] < n:
                    if ans["basic_variables"][j] >= no_of_vars:
                        answer_dict["optimal_solution"][ans["basic_variables"][j] - no_of_vars] -= ans["optimal_tableau"][j+1, 0]
                    else:
                        answer_dict["optimal_solution"][ans["basic_variables"][j]] += ans["optimal_tableau"][j+1, 0]

        f = open("my_output.txt", "w")
        f.write(answer_dict["solution_status"])
        if (answer_dict["solution_status"]=="optimal"):
            f.write(answer_dict["optimal_value"])
        f.close()
        
    return answer_dict

print(simplex_algo())