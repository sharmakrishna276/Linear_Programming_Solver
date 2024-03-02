import pulp
import numpy as np
import random

pulp.LpSolverDefault.msg = 0

def input_generator(max_min, A, b, c, constraint):
    f = open("testcase.txt", "w")
    f.write("[objective]\n")

    if max_min:
        f.write("maximize\n\n[A]\n")
    else:
        f.write("minimize\n\n[A]\n")

    m = len(A)
    n = len(A[0])

    for i in range(m):
        f.write(str(A[i][0]))

        for j in range(1, n):
            f.write(", " + str(A[i][j]))
        f.write("\n")

    f.write("\n[b]\n")

    for i in range(m):
        f.write(str(b[i]) + "\n")

    f.write("\n[contraint_types]\n")

    for i in range(m):
        if constraint[i] == 0:
            f.write("=")
        elif constraint[i] == 1:
            f.write(">=")
        else:
            f.write("<=")

        f.write("\n")

    f.write("\n[c]\n")
    f.write(str(c[0]))

    for i in range(1, n):
        f.write(", " + str(c[i]))

    f.close()


def generate_lp_problem(n, m):
    max_min = random.randint(0, 1)
    c = np.round(np.random.uniform(-100, 100, n), 3)
    A = np.round(np.random.uniform(-100, 100, (m, n)), 3)
    b = np.round(np.random.uniform(-1000, 1000, m), 3)

    if max_min:
        prob = pulp.LpProblem("LP Problem", pulp.LpMaximize)
    else:
        prob = pulp.LpProblem("LP Problem", pulp.LpMinimize)

    x = [pulp.LpVariable(f'x{i}', lowBound=0) for i in range(n)]
    prob += pulp.lpDot(c, x)

    constraint = []

    for i in range(m):
        less_grt = random.randint(0, 2)
        if less_grt == 0:
            prob += pulp.lpDot(A[i], x) == b[i]
        elif less_grt == 1:
            prob += pulp.lpDot(A[i], x) >= b[i]
        else:
            prob += pulp.lpDot(A[i], x) <= b[i]
        constraint.append(less_grt)

    return max_min, prob, A, b, c, constraint


def test_lp_problem(prob):
    prob.solve()
    print("Status:", pulp.LpStatus[prob.status])

    for v in prob.variables():
        print(v.name, "=", v.varValue)

    print("Optimal value:", pulp.value(prob.objective))


n = random.randint(1, 10)
m = random.randint(1, 10)
max_min, prob, A, b, c, constraint = generate_lp_problem(n, m)
print("LP problem generated:")
print(prob)
print("Testing LP problem:")
test_lp_problem(prob)
input_generator(max_min, A, b, c, constraint)