import json
import numpy as np

def simplex_algo(filename):

    ##### improve the part where you check if it is required to add artificial variables or not

    dct = {'[objective]':[],'[A]':[],'[b]':[],'[constraint_types]':[],'[c]':[]}
    ans = {'initial_tableau': None, 'final_tableau': None, 'solution_status': '', 'optimal_solution': None, 'optimal_value': None}

    with open(filename, 'r') as file:
        # Read all lines into a list
        lines = [line.strip() for line in file.readlines()]

        motive=''

        for line in lines:
            if line=='':
                continue
            elif line[0]=='[':
                motive = line
            else:
                if motive=='[A]' or motive=='[b]' or motive=='[c]':
                    lst = json.loads('[' + line + ']')
                    dct[motive].append(lst)
                else:
                    dct[motive].append(line)
        
    var = 'x'
    dct['vars']=[var+str(i) for i in range(len(dct['[A]'][0]))]
    n_initial = len(dct['vars'])

    sgn = [0 for _ in range(len(dct['[A]'][0]))]

    # checking number of sign constraints
    _ = 0
    popy = []
    for lst in dct['[A]']:
        if (lst.count(0)==len(lst)-1 and dct['[constraint_types]'][_]=='>='):
            for i in range(len(lst)):
                if (lst[i]!=0):
                    ind =i
                    break
            popy.append(_)  # constraint to be removed 
            if (lst[ind]<0):
                for i in range(len(lst)):
                    dct['[A]'][i][ind]=(-1)*dct['[A]'][i][ind] #changing <= constraint to >= constraint
            sgn[ind]=1
        _+=1

    for idx in range(len(popy)-1,-1,-1):
        dct['[A]'].pop(popy[idx])
        dct['[b]'].pop(popy[idx])
        dct['[constraint_types]'].pop(popy[idx])

    n_vars = len(dct['[A]'][0])

    # getting rid of free variables
    for i in range(len(sgn)):
        if (sgn[i]==0):
            for j in range(len(dct['[A]'])):
                dct['[A]'][j].append(-1*dct['[A]'][j][i])
            dct['[c]'][0].append(-1*dct['[c]'][0][i])
            sgn.append(1)
            dct['vars'].append(dct['vars'][i] + '-')
            dct['vars'][i] = dct['vars'][i] + '+'
            sgn[i]=1
            
    #adding slack variables            
    _ = 0
    for constraint in dct['[constraint_types]']:
        if constraint=='<=':
            dct['vars'].append(var + str(n_vars))
            n_vars+=1
            dct['[constraint_types]'][_] = '='
            for i in range(len(dct['[A]'])):
                if i==_:
                    dct['[A]'][i].append(1)
                else:
                    dct['[A]'][i].append(0)
            dct['[c]'][0].append(0)
            sgn.append(1)

        elif constraint=='>=':
            dct['vars'].append(var + str(n_vars))
            n_vars+=1
            dct['[constraint_types]'][_] = '='
            for i in range(len(dct['[A]'])):
                if i==_:
                    dct['[A]'][i].append(-1)
                else:
                    dct['[A]'][i].append(0)
            dct['[c]'][0].append(0)
            sgn.append(1)

        else:
            _+=1
            continue

        _+=1

    # changing maximize to minimize

    if (dct['[objective]'][0]=='maximize'):
        dct['[c]'][0]=[(-1)*dct['[c]'][0][i] for i in range(len(dct['[c]'][0]))]
        
    # conversion to numpy arrays

    for key in dct:
        if (key!='vars'):
            dct[key] = np.array(dct[key])

    # 2 phase simplex begins

    # making b>=0
    n_original = len(dct['[A]'][0])
    unbdd = False

    for i in range(len(dct['[b]'])):
        if dct['[b]'][i][0]<0:
            dct['[A]'][i] = (-1)*dct['[A]'][i]
            dct['[b]'][i] = (-1)*dct['[b]'][i]

    #basis dictionary 
    bss = {'vars':[] , 'idx': []}

    m=len(dct['[A]']) #number of constraints
    # checking if it is required to add aritificial variables

    nb = 0 # number of basis vectors obtained
    ids = []

    for i in range(len(dct['[A]'][0])):
        # iterating over all the variables
        if (nb==m):
            break

        column_vector = dct['[A]'][:, i]
        
        cond1 = (np.sum(column_vector==0) == (len(column_vector) - 1))
        is_condition_met = False
        if (cond1):
            for j in range(len(column_vector)):
                if (column_vector[j]>0):
                    ids.append(j)
                    is_condition_met = True
                    break

        if (is_condition_met):
            nb+=1
            bss['idx'].append(i)
            bss['vars'].append(dct['vars'][i])

    flag = True  #whether artificial simplex is required or not
    if (nb==m):
        flag = False

    dct.pop('[constraint_types]')

    if (flag):
        #perform tableau method for artificial simplex and drive artificial variables out of basis

        #new cost vector for phase 1
        c_new = np.zeros_like(dct['[c]'])

        #adding artificial varibles
        var = 'y'

        art = 0

        for i in range(m):
            if i not in ids:
                bss['vars'].append(var+str(art))
                bss['idx'].insert(i,len(dct['[A]'][0]))
                dct['vars'].append(var+str(art))
                col = np.zeros((m,1))
                col[i][0]=1
                dct['[A]'] = np.hstack((dct['[A]'],col))
                c_new = np.append(c_new, [[1]], axis=1) 
                art+=1       

        m = len(dct['[A]'])
        n = len(dct['[A]'][0])

        #initializing tableau
        dct['tableau']= np.zeros((m+1,n+1))
        bss['mtx'] = dct['[A]'][:,bss['idx']]
        B = bss['mtx']
        B_1 = np.copy(bss['mtx'])
        np.fill_diagonal(B_1, 1 / np.diag(B))

        # B_1b
        xB = B_1@dct['[b]']
        dct['tableau'][1:,0]= (B_1@dct['[b]']).reshape(m,)

        # B-1 A
        h = B_1@dct['[A]']
        dct['tableau'][1:,1:] = h

        # basic costs
        cB = c_new[:,bss['idx']]
        # reduced costs
        dct['tableau'][0,1:] = (c_new - cB@h).reshape(n,)

        # negative of current cost
        dct['tableau'][0,0] = (-1*cB@xB)[0][0]

        # applying tableau algorithm

        dct['[c_new]'] = c_new
        ans['initial_tableau'] = dct['tableau'][1:]

        while True:
            found = False
            j = 0
            for idx in range(1,n+1):
                if (dct['tableau'][0,idx]<(-1e-6)):
                    found = True
                    j = idx  #entering column
                    break

            if not found:  # all reduced costs are non-negative
                break  #hence we have achieved optimal cost

            u = dct['tableau'][1:,j]

            if (np.all(u<=0)):
                unbdd = True  # optimal cost is -inf
                break

            ratios = []

            for i in range(len(u)):
                if (u[i]>1e-6):
                    ratios.append(dct['tableau'][i+1,0]/u[i])

            mini = min(ratios)

            l = 0  # minimizing index, exiting column
            for i in range(len(u)):
                if u[i]>1e-6 and dct['tableau'][i+1,0]/u[i] == mini:
                    l = i+1
                    break

            dct['tableau'][l] = dct['tableau'][l]/dct['tableau'][l][j]  # making pivot element = 1

            for k in range(m+1):  # making all other entries of pivot column = 0
                if k!=l:
                    dct['tableau'][k] = dct['tableau'][k] - dct['tableau'][k][j]*dct['tableau'][l]

            #updating basis indices and basis varibles vector
                    
            bss['vars'][l-1] = dct['vars'][j-1]
            bss['idx'][l-1] = j-1

        if (unbdd):
            ans['final_tableau'] = dct['tableau'][1:]
            ans['solution_status'] = 'unbounded'
            print(ans)
            return

        if (abs(dct['tableau'][0][0])>1e-6):
            ans['final_tableau'] = dct['tableau'][1:]
            ans['solution_status'] = 'infeasible'
            print(ans)
            return 

        drive_out = []
        for k in range(len(bss['vars'])):
            if bss['vars'][k][0]=='y':
                drive_out.append(k)

        popy = []
        while (len(drive_out)!=0):  # aritficial variables are there in the basis
            l = drive_out.pop()
            all_zero = True
            for j in range(n_original):
                if dct['tableau'][l+1][j+1]!=0:
                    all_zero=False

                    dct['tableau'][l+1] = dct['tableau'][l+1]/dct['tableau'][l+1][j+1]  # making pivot element = 1
                    for k in range(m+1):  # making all other entries of pivot column = 0
                        if k!=l+1:
                            dct['tableau'][k] = dct['tableau'][k] - dct['tableau'][k][j+1]*dct['tableau'][l+1]

                    bss['vars'][l] = dct['vars'][j]
                    bss['idx'][l] = j

            if (all_zero):
                popy.append(l+1)

        for id in popy:
            dct['[A]'] = np.delete(dct['[A]'] ,id-1, axis = 0)
            dct['[b]'] = np.delete(dct['[b]'] ,id-1, axis =0)
            dct['tableau'] = np.delete(dct['tableau'] ,id, axis=0)
            bss['vars'].pop(id-1)
            bss['idx'].pop(id-1)

        bss['mtx'] = dct['[A]'][:,bss['idx']]
        dct.pop('[c_new]')

        start_idx = 0 ### index where artificial variables begin

        for k in range(len(dct['vars'])):
            if (dct['vars'][k][0]=='y'):
                start_idx=k
                break

        if (k!=n+1):
            tableau = np.delete(dct['tableau'] , np.s_[start_idx+1:], axis=1)
            dct['tableau'] = tableau
            A = np.delete(dct['[A]'], np.s_[start_idx:], axis=1)
            dct['[A]'] = A
            del dct['vars'][start_idx:]

        n = len(dct['[A]'][0])
        h = dct['tableau'][1:,1:] 
        xB = dct['tableau'][1:,0]

        # basic costs
        cB = dct['[c]'][:,bss['idx']]
        # reduced costs
        dct['tableau'][0,1:] = (dct['[c]'] - cB@h).reshape(n,)

        # negative of current cost
        dct['tableau'][0,0] = (-1*cB@xB)[0]

    else: 

        # artificial simplex not required, initialize tableau for normal simplex
        m = len(dct['[A]'])
        n = len(dct['[A]'][0])

        #initializing tableau
        dct['tableau']= np.zeros((m+1,n+1))
        bss['mtx'] = dct['[A]'][:,bss['idx']]
        B = bss['mtx']
        B_1 = np.copy(bss['mtx'])
        np.fill_diagonal(B_1, 1 / np.diag(B))

        # B_1b
        xB = B_1@dct['[b]']
        dct['tableau'][1:,0]= (B_1@dct['[b]']).reshape(m,)

        # B-1 A
        h = B_1@dct['[A]']
        dct['tableau'][1:,1:] = h

        # basic costs
        cB = dct['[c]'][:,bss['idx']]

        # reduced costs
        dct['tableau'][0,1:] = (dct['[c]'] - cB@h).reshape(n,)

        # negative of current cost
        dct['tableau'][0,0] = (-1*cB@xB)[0][0]

        ans['initial_tableau'] = dct['tableau'][1:]

    ######### phase-2 begins ###############
        
    m = len(dct['[A]'])
    n = len(dct['[A]'][0])
    unbdd = False

    while True:
        found = False
        j = 0
        for idx in range(1,n+1):
            if (dct['tableau'][0,idx]<(-1e-6)):
                found = True
                j = idx  #entering column
                break

        if not found:  # all reduced costs are non-negative
            break  #hence we have achieved optimal cost

        u = dct['tableau'][1:,j]

        if (np.all(u<=0)):
            unbdd = True  # optimal cost is -inf
            break

        ratios = []

        for i in range(len(u)):
            if (u[i]>1e-6):
                ratios.append(dct['tableau'][i+1,0]/u[i])

        mini = min(ratios)
        l = 0  # minimizing index, exiting column
        for i in range(len(u)):
            if u[i]>1e-6 and dct['tableau'][i+1,0]/u[i] == mini:
                l = i+1
                break

        dct['tableau'][l] = dct['tableau'][l]/dct['tableau'][l][j]  # making pivot element = 1

        for k in range(m+1):  # making all other entries of pivot column = 0
            if k!=l:
                dct['tableau'][k] = dct['tableau'][k] - dct['tableau'][k][j]*dct['tableau'][l]

        #updating basis indices and basis varibles vector
                
        bss['vars'][l-1] = dct['vars'][j-1]
        bss['idx'][l-1] = j-1

    if (unbdd):
        ans['final_tableau'] = dct['tableau'][1:]
        ans['solution_status'] = 'unbounded'
        print(ans)
        return 
    
    #### optimal solutions reach this point #####
    ans['final_tableau'] = dct['tableau'][1:]
    ans['solution_status'] = 'optimal'
    if (dct['[objective]'][0]=='minimize'):
        ans['optimal_value'] = (-1)*dct['tableau'][0][0]
    else:
        ans['optimal_value'] = dct['tableau'][0][0]

    vec = {}
    var = 'x'
    for i in range(n_initial):
        vec[var + str(i) + '+'] = 0
        vec[var + str(i) + '-'] = 0

    for i in range(len(bss['vars'])):
        if (bss['vars'][i][-1] == '+' or bss['vars'][i][-1] == '-'):
            if bss['vars'][i] in vec:
                vec[bss['vars'][i]] = dct['tableau'][i+1][0]
        else:
            if (bss['vars'][i] + '+') in vec:
                vec[bss['vars'][i] + '+'] = dct['tableau'][i+1][0]

    ans['optimal_solution'] = []

    for i in range(n_initial):
        ans['optimal_solution'].append(vec[var + str(i) + '+'] - vec[var + str(i) + '-'])

    print(ans)
    # f = open("my_output.txt", "w")
    # f.write(ans["solution_status"])
    # if (ans["solution_status"]=="optimal"):
    #     f.write(str(ans["optimal_value"]))
    return

simplex_algo("testcase.txt")