from util4Groebner import *
from math import comb, log2

t = 8 #Number of branches
d = 5 #Number of S-boxes

def Griffin4_initM_model(m, t, var_in, con_in ,deg_in, g_in,h_in, bse_in):
    var_out = {}
    con_out = {}
    deg_out = {}
    g_out = {}
    h_out = {}
    bse_out = {}
    # Adding Output Variables
    for i in range(t):
        var_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_0r_{}".format(i))
        con_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_0r_{}".format(i))
        deg_out[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_0r_{}".format(i))
        g_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_0r_{}".format(i))
        h_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_0r_{}".format(i))
        bse_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_0r_{}".format(i))

    equ_M = m.addVar(lb=0, vtype = GRB.INTEGER, name="equMDS_-1r")
    m.update()

    # MDS constraints
    MDS_module(m, t, [var_in[0], var_in[1], var_in[2], var_in[3],var_in[4], var_in[5], var_in[6], var_in[7], var_out[0], var_out[1], var_out[2], var_out[3],var_out[4], var_out[5], var_out[6], var_out[7],], \
               [con_in[0], con_in[1], con_in[2], con_in[3],con_in[4], con_in[5], con_in[6], con_in[7], con_out[0], con_out[1], con_out[2], con_out[3],con_out[4], con_out[5], con_out[6], con_out[7]], \
               [deg_in[0], deg_in[1], deg_in[2], deg_in[3],deg_in[4], deg_in[5], deg_in[6], deg_in[7], deg_out[0], deg_out[1], deg_out[2], deg_out[3],deg_out[4], deg_out[5], deg_out[6], deg_out[7]], \
               [g_in[0], g_in[1], g_in[2], g_in[3],g_in[4], g_in[5], g_in[6], g_in[7], g_out[0], g_out[1], g_out[2], g_out[3],g_out[4], g_out[5], g_out[6], g_out[7]], \
               [h_in[0], h_in[1], h_in[2], h_in[3],h_in[4], h_in[5], h_in[6], h_in[7], h_out[0], h_out[1], h_out[2], h_out[3],h_out[4], h_out[5], h_out[6], h_out[7]],\
               [bse_in[0], bse_in[1], bse_in[2], bse_in[3],bse_in[4], bse_in[5], bse_in[6], bse_in[7], bse_out[0], bse_out[1], bse_out[2], bse_out[3],bse_out[4], bse_out[5], bse_out[6], bse_out[7]],\
               equ_M,  "{}r_MDS".format(-1))


    return var_out, con_out, deg_out, g_out, h_out, bse_out

def Griffin4_round_model(R, m, r, t, var_in, con_in ,deg_in, g_in,h_in, bse_in):
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
    # Positions 0,1,2,3,4,5,6,7 are pending quaternions + basis binary, directly from the current round of input states
    for i in range(t):
        var[i] = var_in[i]
        con[i] = con_in[i]
        deg[i] = deg_in[i]
        g[i] = g_in[i]
        h[i] = h_in[i]
        bse[i] = bse_in[i]
    # Position 8 is the vcde tuple
    var[8], con[8] = gen_branching_vars_vc(m, r, 8)
    deg[8], g[8], h[8], bse[8] = gen_pending_vars_dghb(m, r, 8)
    # Position 9 is the vcde tuple
    var[9], con[9], deg[9], equ[9]= gen_branching_vars_vcde(m, r, 9)
    # Position 10 is the pending quaternion
    deg[10], g[10], h[10], bse[10] = gen_pending_vars_dghb(m, r, 10)
    # Position 11 is the deg individual variable
    deg[11] = gen_vars_d(m, r, 11)
    # Position 12 is a triad
    var[12], con[12], deg[12] = gen_branching_vars_vcd(m, r, 12)
    # Position 13 is the vcde tuple
    var[13], con[13] = gen_branching_vars_vc(m, r, 13)
    deg[13], g[13], h[13], bse[13] = gen_pending_vars_dghb(m, r, 13)
    # Position 14 is deg individual variables
    deg[14] = gen_vars_d(m, r, 14)
    # Position 15 is a triad
    var[15], con[15], deg[15] = gen_branching_vars_vcd(m, r, 15)
    # Position 16 is the vcde tuple
    var[16], con[16] = gen_branching_vars_vc(m, r, 16)
    deg[16], g[16], h[16], bse[16] = gen_pending_vars_dghb(m, r, 16)
     # Position 17 is the deg individual variable
    deg[17] = gen_vars_d(m, r, 17)
    # Position 18 is a triad
    var[18], con[18], deg[18] = gen_branching_vars_vcd(m, r, 18)
    # Position 19 is the vcde tuple
    var[19], con[19] = gen_branching_vars_vc(m, r, 19)
    deg[19], g[19], h[19], bse[19] = gen_pending_vars_dghb(m, r, 19)
    # Position 20 is the deg individual variable
    deg[20] = gen_vars_d(m, r, 20)
    # Position 21 is a triad.
    var[21], con[21], deg[21] = gen_branching_vars_vcd(m, r, 21)
    # Position 22 is the vcde tuple
    var[22], con[22] = gen_branching_vars_vc(m, r, 22)
    deg[22], g[22], h[22], bse[22] = gen_pending_vars_dghb(m, r, 22)
    # Position 23 is deg individual variable
    deg[23] = gen_vars_d(m, r, 23)
    # Position 24 is a triad
    var[24], con[24], deg[24] = gen_branching_vars_vcd(m, r, 24)
    # Position 25 is the vcde tuple
    var[25], con[25] = gen_branching_vars_vc(m, r, 25)
    deg[25], g[25], h[25], bse[25] = gen_pending_vars_dghb(m, r, 25)
    # Position 26 is deg individual variable
    deg[26] = gen_vars_d(m, r, 26)
    # Position 27 is the triad
    var[27], con[27], deg[27] = gen_branching_vars_vcd(m, r, 27)
    # Position 28 is the vcde tuple
    var[28], con[28] = gen_branching_vars_vc(m, r, 28)
    deg[28], g[28], h[28], bse[28] = gen_pending_vars_dghb(m, r, 28)
    # Add output variable: next round 29,30,31,32,33,34,35,36 positions are pending quaternion + base binary
    for i in [29,30,31,32,33,34,35,36]:
        var[i], con[i] = gen_branching_vars_vc(m, r+1, i-29)
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r+1, i-29)


    # Add the equ variable to each arithmetic module
    # Individual computing modules
    equG2 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG2_{}r".format(r))
    equG1 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG1_{}r".format(r))
    equG3 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG3_{}r".format(r))
    equG4 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG4_{}r".format(r))
    equG5 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG5_{}r".format(r))
    equG6 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG6_{}r".format(r))
    equ69 = m.addVar(lb = 0, vtype=GRB.INTEGER, name="equ69_{}r".format(r))
    equ87 = m.addVar(lb = 0, vtype=GRB.INTEGER, name="equ87_{}r".format(r))
    equMDS = m.addVar(lb=0, vtype=GRB.INTEGER, name = "equMDS_{}r".format(r))
    equ51213 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ51213_{}r".format(r))
    equ41516 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ41516_{}r".format(r))
    equ31819 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ31819_{}r".format(r))
    equ22122 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ22122_{}r".format(r))
    equ12425 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ12425_{}r".format(r))
    equ02728 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ02728_{}r".format(r))
    m.update()

    # Adding constraints to each branch module
    # 0 pending-up; 1 pending-up(s),same as pending-up; 2 pending-up;
    # 5 pending-up(s),same as pending-up
    for i in [0,1,2,3,4,5,6,8]:
        pending_up(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r,i))
    # down-pending
    for i in [7,13,16,19,22,25,28]:
        down_pending(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r, i))
    # 7, 10 down-up
    for i in [12,15,18,21,24,27]:
        down_up(m, var[i], con[i], deg[i])
    # 9-10-11-14-17-20-23-26 down-pending-ups
    down_pending_ups(m, deg[9],deg[10], g[10],[deg[11], deg[14],deg[17],deg[20],deg[23],deg[26]], var[9], con[9],equ[9], "{}r_9".format(r))

    # Adding constraints to the algorithms
    #0,7 8 multiplication
    multiplication_module(m, [deg[0], deg[27]], var[28], con[28], deg[28], equ02728, "{}r_02728M".format(r))
    # 1,10,11 multiplication
    multiplication_module(m, [deg[1], deg[24]], var[25], con[25], deg[25], equ12425, "{}r_12425M".format(r))
    multiplication_module(m, [deg[2], deg[21]], var[22], con[22], deg[22], equ22122, "{}r_22122M".format(r))
    multiplication_module(m, [deg[3], deg[18]], var[19], con[19], deg[19], equ31819, "{}r_31819M".format(r))
    multiplication_module(m, [deg[4], deg[15]], var[16], con[16], deg[16], equ41516, "{}r_41516M".format(r))
    multiplication_module(m, [deg[5], deg[12]], var[13], con[13], deg[13], equ51213, "{}r_51213M".format(r))
    # 2,4 Univariate nonlinear modules
    NLUP_module(m, [var[6], var[9]],[con[6], con[9]], [deg[6], deg[9]], equ69, "{}r_S1".format(r), d)
    # 5,3 Univariate nonlinear modules
    NLUP_module(m, [var[8], var[7]], [con[8], con[7]], [deg[8], deg[7]], equ87, "{}r_S2".format(r), d)
    # G2
    NLMP_module(m, [deg[5], deg[8], deg[14]], var[15], con[15], deg[15], equG2, "{}r_G2".format(r), 2)
    # G1
    NLMP_module(m, [deg[8], deg[11]], var[12], con[12], deg[12], equG1, "{}r_G1".format(r), 2)
    NLMP_module(m, [deg[8], deg[17], deg[4]], var[18], con[18], deg[18], equG3, "{}r_G3".format(r), 2)
    NLMP_module(m, [deg[8], deg[20], deg[3]], var[21], con[21], deg[21], equG4, "{}r_G4".format(r), 2)
    NLMP_module(m, [deg[8], deg[23], deg[2]], var[24], con[24], deg[24], equG5, "{}r_G5".format(r), 2)
    NLMP_module(m, [deg[8], deg[26], deg[1]], var[27], con[27], deg[27], equG6, "{}r_G6".format(r), 2)
    # MDS

    MDS_module(m, 8, [var[28],var[25],var[22],var[19], var[16], var[13], var[9], var[8],var[29],var[30],var[31],var[32], var[33], var[34], var[35], var[36]], \
               [con[28], con[25], con[22], con[19], con[16], con[13], con[9], con[8],con[29], con[30], con[31], con[32], con[33], con[34], con[35], con[36]],\
               [deg[28],deg[25],deg[22], deg[19],deg[16],deg[13],deg[10],deg[8],deg[29],deg[30],deg[31], deg[32],deg[33],deg[34],deg[35],deg[36]],\
               [g[28],g[25],g[22], g[19],g[16],g[13],g[10],g[8],g[29],g[30],g[31], g[32],g[33],g[34],g[35],g[36]],\
               [h[28],h[25],h[22], h[19],h[16],h[13],h[10],h[8],h[29],h[30],h[31], h[32],h[33],h[34],h[35],h[36]],\
               [bse[28],bse[25],bse[22], bse[19],bse[16],bse[13],bse[10],bse[8],bse[29],bse[30],bse[31], bse[32],bse[33],bse[34],bse[35],bse[36]],\
               equMDS,  "{}r_MDS".format(r))

    # Returns the output variable of the round
    var_out = [var[29],var[30],var[31],var[32], var[33], var[34], var[35], var[36]]
    con_out = [con[29], con[30], con[31], con[32], con[33], con[34], con[35], con[36]]
    deg_out = [deg[29],deg[30],deg[31], deg[32],deg[33],deg[34],deg[35],deg[36]]
    g_out = [g[29],g[30],g[31], g[32],g[33],g[34],g[35],g[36]]
    h_out = [h[29],h[30],h[31], h[32],h[33],h[34],h[35],h[36]]
    bse_out = [bse[29],bse[30],bse[31], bse[32],bse[33],bse[34],bse[35],bse[36]]

    m.update()
    return var_out, con_out, deg_out, g_out, h_out, bse_out

def parsing_Griffin4_solution(m, R):
    # Initial MDS
    print('------------r = {}-------------'.format(-1))
    print('var_-1r', end=" ")
    for j in range(t):
        print(round(m.getVarByName("var_-1r_{}".format(j)).X), end = " ")
    print('\nvar_0r ', end=" ")
    for j in range(t):
        print(round(m.getVarByName("var_0r_{}".format(j)).X), end=" ")
    print('\ncon_-1r', end=" ")
    for j in range(t):
        print(round(m.getVarByName("con_-1r_{}".format(j)).X), end=" ")
    print('\ncon_0r ', end=" ")
    for j in range(t):
        print(round(m.getVarByName("con_0r_{}".format(j)).X), end=" ")
    print('\ndeg_-1r', end=" ")
    for j in range(t):
        print(round(m.getVarByName("deg_-1r_{}".format(j)).X), end=" ")
    print('\ndeg_0r ', end=" ")
    for j in range(t):
        print(round(m.getVarByName("deg_0r_{}".format(j)).X), end=" ")
    print('\nbse_-1r', end=" ")
    for j in range(t):
        print(round(m.getVarByName("bse_-1r_{}".format(j)).X), end=" ")
    print('\nbse_0r ', end=" ")
    for j in range(t):
        print(round(m.getVarByName("bse_0r_{}".format(j)).X), end=" ")
    print('\ng_-1r', end=" ")
    for j in range(t):
        print(round(m.getVarByName("g_-1r_{}".format(j)).X), end=" ")
    print('\ng_0r ', end=" ")
    for j in range(t):
        print(round(m.getVarByName("g_0r_{}".format(j)).X), end=" ")
    print("\nequMDS_-1r", round(m.getVarByName("equMDS_-1r").X))
    # each round
    for r in range(R):
        print('--------------r = {}-----------------'.format(r))
        print('var_{}r'.format(r), end=" ")
        for j in [28,25,22,19,16,13,9,8]:
            print(round(m.getVarByName("var_{}r_{}".format(r,j)).X), end=" ")
        print('\nvar_{}r'.format(r+1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("var_{}r_{}".format(r+1, j)).X), end=" ")
        print('\ncon_{}r'.format(r), end=" ")
        for j in [28,25,22,19,16,13,9,8]:
            print(round(m.getVarByName("con_{}r_{}".format(r, j)).X), end=" ")
        print('\ncon_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("con_{}r_{}".format(r + 1, j)).X), end=" ")
        print('\ndeg_{}r'.format(r), end=" ")
        for j in [28,25,22,19,16,13,10,8]:
            print(round(m.getVarByName("deg_{}r_{}".format(r, j)).X), end=" ")
        print('\ndeg_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("deg_{}r_{}".format(r + 1, j)).X), end=" ")
        print('\nbse_{}r'.format(r), end=" ")
        for j in [28,25,22,19,16,13,10,8]:
            print(round(m.getVarByName("bse_{}r_{}".format(r, j)).X), end=" ")
        print('\nbse_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("bse_{}r_{}".format(r + 1, j)).X), end=" ")
        print('\ng_{}r'.format(r), end=" ")
        for j in [28,25,22,19,16,13,10,8]:
            print(round(m.getVarByName("g_{}r_{}".format(r, j)).X), end=" ")
        print('\ng_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("g_{}r_{}".format(r + 1, j)).X), end=" ")
        print("\nequMDS_{}r".format(r), round(m.getVarByName("equMDS_{}r".format(r)).X))

        for j in [27,24,21,18,15,12]:
            print("var_{}r_{} = ".format(r, j), end = " ")
            print(round(m.getVarByName("var_{}r_{}".format(r,j)).X))
        for j in [27,24,21,18,15,12]:
            print("con_{}r_{} = ".format(r, j), end = " ")
            print(round(m.getVarByName("con_{}r_{}".format(r,j)).X))
    for var in m.getVars():
        if "equ" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", var.X)
def gen_Griffin4_model(R, t, numVar):
    # R: number of rounds; t: number of branches; numVar: number of introduced variables

    # Creating Models
    m = gp.Model("Griffin4_CICO_MILP_{}r_{}t_{}v".format(R, t, numVar))
    # Create input variables: pending quaternion + input state binary
    var_in={};con_in={};deg_in={};g_in={};h_in={};bse_in={};
    for i in range(t):
        var_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_-1r_{}".format(i))
        con_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_-1r_{}".format(i))
        deg_in[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_-1r_{}".format(i))
        g_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_-1r_{}".format(i))
        h_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_-1r_{}".format(i))
        bse_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_-1r_{}".format(i))
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
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） constraints
    for i in range(t):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # Initial MDS transformation
    var_in, con_in, deg_in, g_in, h_in, bse_in = Griffin4_initM_model(m, t, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # Modeling per round
    for r in range(R):
        var_in, con_in, deg_in, g_in, h_in, bse_in = \
            Griffin4_round_model(R, m, r, t, var_in, con_in, deg_in, g_in, h_in, bse_in)

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
    repeated_con_vars = [] # Take only one constant before and after the S-box and remove the duplicates
    for r in range(R):
        repeated_con_vars.append(m.getVarByName("con_{}r_9".format(r)))
        repeated_con_vars.append(m.getVarByName("con_{}r_8".format(r)))
    sum_cons = sum(all_con_vars) - sum(repeated_con_vars)
    m.addConstr(sum_cons == t) #The number of constants is equal to t, then on average there is 1 solution

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

    # Artificially added constraints
    '''
    m.addConstr(m.getVarByName("var_0r_1")==1)
    m.addConstr(m.getVarByName("var_0r_2") == 1)
    m.addConstr(m.getVarByName("var_0r_5") == 1)

    m.addConstr(m.getVarByName("var_1r_2")==1)
    m.addConstr(m.getVarByName("var_1r_5") == 1)

    m.addConstr(m.getVarByName("var_2r_2") == 1)
    m.addConstr(m.getVarByName("var_2r_5") == 1)

    m.addConstr(m.getVarByName("var_3r_2") == 1)
    m.addConstr(m.getVarByName("var_3r_5") == 1)

    m.addConstr(m.getVarByName("con_4r_2") == 1)
    m.addConstr(m.getVarByName("con_4r_4") == 1)
    m.addConstr(m.getVarByName("con_4r_3") == 1)
    m.addConstr(m.getVarByName("con_4r_5") == 1)

    m.addConstr(m.getVarByName("var_5r_5") == 1)
    '''

    # Model settings
    m.setParam("TimeLimit", 600) # Set the run to stop after 600 seconds and output the currently found optimal solution
    # m.setParam("BestObjStop", 100) #Set the objective function to stop when it falls below 100
    m.write("Griffin4_CICO_MILP_{}rrrrrr_{}t_{}v.lp".format(R, t, numVar))
    m.optimize()
    m.write("Griffin4_MILP_{}r_{}t_{}v.sol".format(R, t, numVar))
    parsing_Griffin4_solution(m, R)
    print("Number of branches：{}，".format(t))
    print('round：{}，'.format(R))
    print('Number of variables：{}，'.format(numVar))
    print('Sum of the number of equations：{}'.format(obj.X))
    print('solution complexity：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))

if __name__ == '__main__':
    gen_Griffin4_model(4, 8, 8)
