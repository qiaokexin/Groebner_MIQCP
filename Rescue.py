from util4Groebner import *
from math import comb, log2


d = 3  # Number of S-boxes



def Rescue_round_model(R, m, r, t, var_in, con_in, deg_in, g_in, h_in, bse_in):
    # Modeling the rth round round function, t is the number of branches

    # Initializing the auxiliary variable dictionary
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # Variable Settings
    # Positions 0,1,2 are pending quaternions + base binary, directly from the current round of input states
    for i in range(t):
        var[i] = var_in[i]
        con[i] = con_in[i]
        deg[i] = deg_in[i]
        g[i] = g_in[i]
        h[i] = h_in[i]
        bse[i] = bse_in[i]

    # Positions t, t+1, ..., 5*t-1 are base binary + pending quaternion
    for i in range(t, 5*t):
        var[i], con[i] = gen_branching_vars_vc(m, r, i)
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r, i)
   

    # Add the equ variable to each arithmetic module
    # modules
    for i in range(t):
        equ[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ{}{}_{}r".format(i, i+t, r))
    for i in range(2*t, 3*t):
        equ[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ{}{}_{}r".format(i, i+t, r))

    equMDS1 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS1_{}r".format(r))
    equMDS2 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS2_{}r".format(r))
    m.update()

    # Adding constraints to each branch module
    for i in range(t, 3*t):
        pending_up(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r, i))
    
    for i in range(t):
        down_pending(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r, i))
    for i in range(3*t, 4*t):
        down_pending(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r, i))
    
    # Adding constraints to the algorithms

    for i in range(t):
        # i, i+t Univariate nonlinear modules
        NLUP_module(m, [var[i+t], var[i]], [con[i+t], con[i]], [deg[i+t], deg[i]], equ[i], "{}r_S{}".format(r, i+1), d)
        # i+2*t, i+3*t Univariate nonlinear modules
        NLUP_module(m, [var[i+2*t], var[i+3*t]], [con[i+2*t], con[i+3*t]], [deg[i+2*t], deg[i+3*t]], equ[i+2*t], "{}r_S{}".format(r, i+t+1), d)
    


    # MDS


        # 生成通用索引列表：t~2t-1 与 2t~3t-1
    indices = [t+k for k in range(t)] + [2*t + k for k in range(t)]

    # 第一次调用 MDS_module（对应原 t, 2t 组）
    MDS_module(m, t,
            [var[idx] for idx in indices],
            [con[idx] for idx in indices],
            [deg[idx] for idx in indices],
            [g[idx] for idx in indices],
            [h[idx] for idx in indices],
            [bse[idx] for idx in indices],
            equMDS1, "{}r_MDS1".format(r))

    # 第二次调用 MDS_module（对应原 3t, 4t 组）
    indices2 = [3*t + k for k in range(t)] + [4*t + k for k in range(t)]
    MDS_module(m, t,
            [var[idx] for idx in indices2],
            [con[idx] for idx in indices2],
            [deg[idx] for idx in indices2],
            [g[idx] for idx in indices2],
            [h[idx] for idx in indices2],
            [bse[idx] for idx in indices2],
            equMDS2, "{}r_MDS2".format(r))
    # MDS_module(m, t, [var[t], var[1+t], var[2+t], var[3+t], var[4+t], var[5+t], var[2*t], var[1+2*t], var[2+2*t], var[3+2*t], var[4+2*t], var[5+2*t]], \
    #             [con[t], con[1+t], con[2+t], con[3+t], con[4+t], con[5+t], con[2*t], con[1+2*t], con[2+2*t], con[3+2*t], con[4+2*t], con[5+2*t]], \
    #             [deg[t], deg[1+t], deg[2+t], deg[3+t], deg[4+t], deg[5+t], deg[2*t],deg[1+2*t], deg[2+2*t], deg[3+2*t], deg[4+2*t], deg[5+2*t]], \
    #             [g[t], g[1+t], g[2+t], g[3+t], g[4+t], g[5+t], g[2*t], g[1+2*t], g[2+2*t], g[3+2*t], g[4+2*t], g[5+2*t] ], \
    #             [h[t], h[1+t], h[2+t], h[3+t], h[4+t], h[5+t], h[2*t], h[1+2*t], h[2+2*t], h[3+2*t], h[4+2*t], h[5+2*t] ], \
    #             [bse[t], bse[1+t], bse[2+t], bse[3+t], bse[4+t], bse[5+t], bse[2*t], bse[1+2*t], bse[2+2*t], bse[3+2*t], bse[4+2*t], bse[5+2*t] ], \
    #             equMDS1, "{}r_MDS1".format(r))
    # MDS_module(m, t, [var[3*t], var[1+3*t], var[2+3*t], var[3+3*t], var[4+3*t], var[5+3*t], var[4*t], var[1+4*t], var[2+4*t], var[3+4*t], var[4+4*t], var[5+4*t]], \
    #             [con[3*t], con[1+3*t], con[2+3*t], con[3+3*t], con[4+3*t], con[5+3*t], con[4*t], con[1+4*t], con[2+4*t], con[3+4*t], con[4+4*t], con[5+4*t]], \
    #             [deg[3*t], deg[1+3*t], deg[2+3*t], deg[3+3*t], deg[4+3*t], deg[5+3*t], deg[4*t],deg[1+4*t], deg[2+4*t], deg[3+4*t], deg[4+4*t], deg[5+4*t]], \
    #             [g[3*t], g[1+3*t], g[2+3*t], g[3+3*t], g[4+3*t], g[5+3*t], g[4*t], g[1+4*t], g[2+4*t], g[3+4*t], g[4+4*t], g[5+4*t] ], \
    #             [h[3*t], h[1+3*t], h[2+3*t], h[3+3*t], h[4+3*t], h[5+3*t], h[4*t], h[1+4*t], h[2+4*t], h[3+4*t], h[4+4*t], h[5+4*t] ], \
    #             [bse[3*t], bse[1+3*t], bse[2+3*t], bse[3+3*t], bse[4+3*t], bse[5+3*t], bse[4*t], bse[1+4*t], bse[2+4*t], bse[3+4*t], bse[4+4*t], bse[5+4*t] ], \
    #             equMDS2, "{}r_MDS2".format(r))

    # Returns the output variable of the round
    var_out = []
    con_out = []
    deg_out = []
    g_out = []
    h_out = []
    bse_out = []
    for i in range(4*t, 5*t):
        var_out.append(var[i])
        con_out.append(con[i])
        deg_out.append(deg[i])  
        g_out.append(g[i])
        h_out.append(h[i])
        bse_out.append(bse[i])

    m.update()
    return var_out, con_out, deg_out, g_out, h_out, bse_out

        
def parsing_Rescue_solution(m, R, t):
    # each round
    for r in range(R):
        print('--------------r = {}-----------------'.format(r))
        
        # 输出 var
        print('\nvar_{}r'.format(r), end=" ")
        for j in range(t, 5 * t):  # 通用：t ~ 5t-1
            print(round(m.getVarByName("var_{}r_{}".format(r, j)).X), end=" ")
        
        # 输出 con
        print('\ncon_{}r'.format(r), end=" ")
        for j in range(t, 5 * t):
            print(round(m.getVarByName("con_{}r_{}".format(r, j)).X), end=" ")
        
        # 输出 deg
        print('\ndeg_{}r'.format(r), end=" ")
        for j in range(t, 5 * t):
            print(round(m.getVarByName("deg_{}r_{}".format(r, j)).X), end=" ")
        
        # 输出 bse
        print('\nbse_{}r'.format(r), end=" ")
        for j in range(t, 5 * t):
            print(round(m.getVarByName("bse_{}r_{}".format(r, j)).X), end=" ")
        
        # 输出 g ########################### 这里修复缩进！
        print('\ng_{}r'.format(r), end=" ")
        for j in range(t, 5 * t):
            print(round(m.getVarByName("g_{}r_{}".format(r, j)).X), end=" ")
        
        
        # 换行收尾
        print()


    for var in m.getVars():
        if "equ" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", var.X)


def gen_Rescue_model(R, t, numVar):
    # R: number of rounds; t: number of branches; numVar: number of introduced variables

    # Creating Models
    m = gp.Model("Rescue_CICO_MILP_{}r_{}t_{}v".format(R, t, numVar))
    # Create input variables: pending quaternion + input state binary
    var_in = {};
    con_in = {};
    deg_in = {};
    g_in = {};
    h_in = {};
    bse_in = {};
    for i in range(t):
        var_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_0r_{}".format(i))
        con_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_0r_{}".format(i))
        deg_in[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_0r_{}".format(i))
        g_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_0r_{}".format(i))
        h_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_0r_{}".format(i))
        bse_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_0r_{}".format(i))
    m.update()

    # Input Boundary Constraints
    # （1/3） Input constraints in CICO
    m.addConstr(con_in[0] == 1)
    # （2/3）Input boundaries add the constraint that var=0,con=0 then g=0
    varPcon = {}
    for i in range(t):
       varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_MDS_{}".format(-1, i))
    m.update()
    for i in range(t):
       m.addConstr(varPcon[i] == var_in[i] + con_in[i])
    m.update()
    # （3/3） constraints
    for i in range(t):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

   
    # Modeling per round
    for r in range(R):
        var_in, con_in, deg_in, g_in, h_in, bse_in = \
            Rescue_round_model(R, m, r, t, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # Constraints between output edges
    # （1/3）Output constraints in CICO
    m.addConstr(con_in[t - 1] == 1)
    # （2/3）Output bounds add the constraint that var=0,con=0 then g=0
    for i in range(t):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_{}".format(R, i))
    m.update()
    for i in range(t):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） communal restraint
    for i in range(t):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # Number of constants
    all_con_vars = [var for var in m.getVars() if "con_" in var.VarName[:4]]
    repeated_con_vars = []  # Take only one constant before and after the S-box and remove the duplicates
    for r in range(R):
        for j in range(t,2*t):
            repeated_con_vars.append(m.getVarByName("con_{}r_{}".format(r , j)))
        for j in range(3*t,4*t):
            repeated_con_vars.append(m.getVarByName("con_{}r_{}".format(r , j)))
    sum_cons = sum(all_con_vars) - sum(repeated_con_vars)
    m.addConstr(sum_cons == t)  # The number of constants is equal to t, then on average there is 1 solution

    # Given the number of variables
    all_vars = [var for var in m.getVars() if "var_" in var.VarName]
    m.addConstr(sum(all_vars) == numVar)

    # objective function
    # Sum of the number of equations
    all_equs = [var for var in m.getVars() if "equ" in var.VarName]
    obj = m.addVar(vtype=GRB.INTEGER, name="obj")
    m.update()
    m.addConstr(obj == sum(all_equs))
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()


    # Model setup
    m.setParam("TimeLimit", 600)  # Set the run to stop after 600 seconds and output the currently found optimal solution
    # m.setParam("BestObjStop", 100) #Set the objective function to stop when it falls below 100
    m.write("Rescue_CICO_MILP_{}rrrrrr_{}t_{}v.lp".format(R, t, numVar))
    m.optimize()
    m.write("Rescue_MILP_{}r_{}t_{}v.sol".format(R, t, numVar))
    parsing_Rescue_solution(m, R, t)
    print("Number of branches：{}，".format(t))
    print('round：{}，'.format(R))
    print('Number of variables：{}，'.format(numVar))
    print('Sum of the number of equations：{}'.format(obj.X))
    print('solution complexity：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))
    print('design Number of variables：{}，'.format(t * R + 1))
    print('design Sum of the number of equations：{}'.format((t * R + 1)* d))
    print('design complexity：2 ^ {}'.format(log2(comb((t * R + 1)* d+1, t * R +1) ** 2)))


if __name__ == '__main__':
    gen_Rescue_model(4, 12, 36)