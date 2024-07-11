from util4Groebner import *
from math import comb, log2

def Cinimion_first_round_model(m, var_in, con_in ,deg_in, g_in,h_in, bse_in):
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
    for i in range(3):
        var[i] = var_in[i]
        con[i] = con_in[i]
        deg[i] = deg_in[i]
        g[i] = g_in[i]
        h[i] = h_in[i]
        bse[i] = bse_in[i]
    # Position 3 is the base binary + pending quaternion
    var[3], con[3] = gen_branching_vars_vc(m, 0, 3)
    deg[3], g[3], h[3], bse[3] = gen_pending_vars_dghb(m, 0, 3)
    # Positions 4,5,6,7,8 are pendent quaternions.
    for i in [4,5,6,7,8]:
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, 0, i)
    # Position 9 is the vce ternary.
    var[9], con[9], equ[9] = gen_branching_vars_vce(m, 0, 9)
    # Position 10 is the vcde quaternion, named according to the rules of the next round
    var[10], con[10], deg[10], equ[10] = gen_branching_vars_vcde(m, 1, 0)
    # Position 11 is the vcde quaternion, named according to the rules of the next round
    var[11], con[11], deg[11], equ[11] = gen_branching_vars_vcde(m, 1, 1)
    # Position 12 is the PENDING quaternion, named according to the next round of rules
    deg[12], g[12], h[12], bse[12] = gen_pending_vars_dghb(m, 1, 2)
    # Position 13 is the pending quaternion, named according to the rules of the next round. 
    deg[13], g[13], h[13], bse[13] = gen_pending_vars_dghb(m, 1, 14)

    # Add the equ variable to each arithmetic module
    # modules
    equMUL = m.addVar(lb=0, vtype = GRB.INTEGER, name="equMUL_{}r".format(0))
    equ234 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ234_{}r".format(0))
    equ156 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ156_{}r".format(0))
    equ078 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ078_{}r".format(0))
    m.update()

    # Adding constraints to the algorithms

    #0-1-3multiplication
    multiplication_module(m, [deg[0], deg[1]], var[3], con[3], deg[3], equMUL, "0r_MUL")
    # 2-3-4addition
    taddition_module(m, 3, [var[2], var[3], var[10]], \
                     [con[2], con[3], con[10]],\
                     [deg[2], deg[3], deg[4]], [g[2], g[3], g[4]],\
                     [h[2], h[3], h[4]], [bse[2], bse[3], bse[4]],\
                     equ234, "0r_234")
    # 1-5-6addition
    taddition_module(m, 3, [var[1], var[10], var[9]],\
                     [con[1], con[10], con[9]],\
                     [deg[1], deg[5], deg[6]],\
                     [g[1], g[5], g[6]],\
                     [h[1], h[5], h[6]],\
                     [bse[1], bse[5], bse[6]], equ156, "0r_156")
    # 0-7-8 addition
    taddition_module(m, 3, [var[0], var[9], var[11]], \
                     [con[0], con[9], con[11]], \
                     [deg[0], deg[7], deg[8]], \
                     [g[0], g[7], g[8]], \
                     [h[0], h[7], h[8]], \
                     [bse[0], bse[7], bse[8]], equ078, "0r_078")

    # Adding constraints to each branch module
    # 3 down-pending
    down_pending(m, var[3], con[3], g[3], deg[3], "0r_3")
    # 4-5-10-13 pending-pending-pending-up
    pending_pending_pending_up(m, deg[4], deg[5], deg[13], g[4], g[5], g[13], var[10], con[10], equ[10], deg[10], "0r_451013")
    # 6-7-12-9 pending-pending-pending
    pending_pending_pending(m, deg[6], deg[7], deg[12], g[6], g[7], g[12], var[9], con[9], equ[9], "0r_67129")

    # Returns the output variable of the round
    # Branch 0 vcde variables
    out0 = [var[10], con[10], deg[10], equ[10]]
    # Branch 1 vcde variables
    out1 = [var[11], con[11], deg[11], equ[11]]
    # Branch 2 pending quaternion
    out2 = [deg[12], g[12], h[12], bse[12]]
    # 8 , 9 , and 14  positions will also need to be used in the next round and will also need to be returned to the
    out8 = [deg[8], g[8], h[8], bse[8]]
    out9 = [var[9], con[9], equ[9]]
    out13 = [deg[13], g[13], h[13], bse[13]]

    return out0, out1, out2, out8, out9, out13

def Cinimion_middle_round_model(m, r, in0, in1, in2, in8, in9, in13):
    # Initializing the auxiliary variable dictionary
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # Variable Settings
    # Positions 0,1 are vcde quaternions, directly from the current round of input states
    var[0], con[0], deg[0], equ[0] = in0
    var[1], con[1], deg[1], equ[1] = in1
    # Position 2 is the pending quaternion, which comes directly from the current round of input states
    deg[2], g[2], h[2], bse[2] = in2
    # Position 14 is the PENDING quaternion, which comes directly from the input
    deg[14], g[14], h[14], bse[14] = in13
    # Position 3 is the base binary + pending quaternion
    var[3], con[3] = gen_branching_vars_vc(m, r, 3)
    deg[3], g[3], h[3], bse[3] = gen_pending_vars_dghb(m, r, 3)
    # Positions 4,5,6,7,8,10,11 are pending quaternions
    for i in [4,5,6,7,8,15]:
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r, i)
    # Position 9 is the vce ternary.
    var[9], con[9], equ[9] = gen_branching_vars_vce(m, r, 9)
    # Position 10 is the vcde quaternion, named according to the rules of the next round
    var[10], con[10], deg[10], equ[10] = gen_branching_vars_vcde(m, r+1, 0)
    # Position 11 is the vcde quaternion, named according to the rules of the next round
    var[11], con[11], deg[11], equ[11] = gen_branching_vars_vcde(m, r+1, 1)
    # Position 12 is the PENDING quaternion, named according to the next round of rules
    deg[12], g[12], h[12], bse[12] = gen_pending_vars_dghb(m, r+1, 2)
    # Position 13 is the pending quaternion, named according to the rules of the next round. This was added later
    deg[13], g[13], h[13], bse[13] = gen_pending_vars_dghb(m, r+1, 14)

    # Add the equ variable to each arithmetic module
    # Modules of operation
    equMUL = m.addVar(lb=0, vtype = GRB.INTEGER, name="equMUL_{}r".format(r))
    equ234 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ234_{}r".format(r))
    equ1356 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ1356_{}r".format(r))
    equ1478 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ1478_{}r".format(r))
    m.update()

    # Adding constraints to the algorithms


    #0-1-3multiplication
    multiplication_module(m, [deg[0], deg[1]], var[3], con[3], deg[3], equMUL, "{}r_MUL".format(r))
    # 2-3-4addition
    taddition_module(m, 3, [in9[0], var[3], var[10]], \
                     [in9[1], con[3], con[10]],\
                     [deg[2], deg[3], deg[4]], [g[2], g[3], g[4]],\
                     [h[2], h[3], h[4]], [bse[2], bse[3], bse[4]],\
                     equ234, "{}r_234".format(r))
    # 15-5-6addition
    taddition_module(m, 3, [var[1], var[10], var[9]],\
                     [con[1], con[10], con[9]],\
                     [deg[15], deg[5], deg[6]],\
                     [g[15], g[5], g[6]],\
                     [h[15], h[5], h[6]],\
                     [bse[15], bse[5], bse[6]], equ1356, "{}r_1356".format(r))
    # 14-7-8 addition
    taddition_module(m, 3, [var[0], var[9], var[11]], \
                     [con[0], con[9], con[11]], \
                     [deg[14], deg[7], deg[8]], \
                     [g[14], g[7], g[8]], \
                     [h[14], h[7], h[8]], \
                     [bse[14], bse[7], bse[8]], equ1478, "{}r_1478".format(r))
    
    # Adding constraints to each branch module

    # pre8-1-15 pending-pending-up
    pending_pending_up(m, in8[0], deg[15], in8[1], g[15], var[1], con[1], deg[1], equ[1], "{}r_pre8115".format(r))

    # 3 down-pending
    down_pending(m, var[3], con[3], g[3], deg[3], "{}r_3".format(r))
    # 4-5-10-13 pending-pending-pending-up
    pending_pending_pending_up(m, deg[4], deg[5], deg[13], g[4], g[5], g[13], var[10], con[10], equ[10], deg[10], "{}r_451013".format(r))
    # 6-7-12-9 pending-pending-pending
    pending_pending_pending(m, deg[6], deg[7], deg[12], g[6], g[7], g[12], var[9], con[9], equ[9], "{}r_67129".format(r))

    # Returns the output variable of the round
    # Branch 0 vcde variables
    out0 = [var[10], con[10], deg[10], equ[10]]
    # Branch 1 vcde variables
    out1 = [var[11], con[11], deg[11], equ[11]]
    # Branch 2 pending quaternion
    out2 = [deg[12], g[12], h[12], bse[12]]
    # Positions 8, 9, and 4 will also need to be used in the next round and returned as well
    out8 = [deg[8], g[8], h[8], bse[8]]
    out9 = [var[9], con[9], equ[9]]
    out13 = [deg[13], g[13], h[13], bse[13]]

    return out0, out1, out2, out8, out9, out13

def Cinimion_last_round_model(m, r, in0, in1, in2, in8, in9, in13):
    # Initializing the auxiliary variable dictionary
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # Variable Settings
    # Positions 0,1 are vcde quaternions, directly from the current round of input states
    var[0], con[0], deg[0], equ[0] = in0
    var[1], con[1], deg[1], equ[1] = in1
    # Position 2 is the pending quaternion, which comes directly from the current round of input states
    deg[2], g[2], h[2], bse[2] = in2
    # Position 14 is the PENDING quaternion, which comes directly from the input
    deg[14], g[14], h[14], bse[14] = in13
    # Position 3 is the base binary + pending quaternion
    var[3], con[3] = gen_branching_vars_vc(m, r, 3)
    deg[3], g[3], h[3], bse[3] = gen_pending_vars_dghb(m, r, 3)
    # Positions 4,5,6,7,8,10,11 are pending quaternions
    for i in [4, 5, 6, 7, 8, 15]:
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r, i)
    # Addition of base binary at position 8
    var[8], con[8] = gen_branching_vars_vc(m, r, 8)
    # Position 9 is the vce ternary.
    var[9], con[9], equ[9] = gen_branching_vars_vce(m, r, 9)
    # Position 10 is the vce triad
    var[10], con[10], equ[10] = gen_branching_vars_vce(m, r, 10)
    # Positions 11, 12 and 13 are not required.

    # Add the equ variable to each arithmetic module
    # modules
    equMUL = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMUL_{}r".format(r))
    equ234 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ234_{}r".format(r))
    equ1356 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1356_{}r".format(r))
    equ1478 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1478_{}r".format(r))
    m.update()

    # Adding constraints to the algorithms


    # 0-1-3multiplication
    multiplication_module(m, [deg[0], deg[1]], var[3], con[3], deg[3], equMUL, "{}r_MUL".format(r))
    # 2-3-4addition
    taddition_module(m, 3, [in9[0], var[3], var[10]], \
                     [in9[1], con[3], con[10]], \
                     [deg[2], deg[3], deg[4]], [g[2], g[3], g[4]], \
                     [h[2], h[3], h[4]], [bse[2], bse[3], bse[4]], \
                     equ234, "{}r_234".format(r))
    # 15-5-6addition
    taddition_module(m, 3, [var[1], var[10], var[9]], \
                     [con[1], con[10], con[9]], \
                     [deg[15], deg[5], deg[6]], \
                     [g[15], g[5], g[6]], \
                     [h[15], h[5], h[6]], \
                     [bse[15], bse[5], bse[6]], equ1356, "{}r_1556".format(r))
    # 14-7-8 addition
    taddition_module(m, 3, [var[0], var[9], var[8]], \
                     [con[0], con[9], con[8]], \
                     [deg[14], deg[7], deg[8]], \
                     [g[14], g[7], g[8]], \
                     [h[14], h[7], h[8]], \
                     [bse[14], bse[7], bse[8]], equ1478, "{}r_1478".format(r))

    # Adding constraints to each branch module
    # pre8-1-15 pending-pending-up
    pending_pending_up(m, in8[0], deg[15], in8[1], g[15], var[1], con[1], deg[1], equ[1], "{}r_pre8115".format(r))
    # 3 down-pending
    down_pending(m, var[3], con[3], g[3], deg[3], "{}r_3".format(r))
    # 4-5-10 pending-pending
    pending_pending(m, deg[4], deg[5], g[4], g[5], var[10], con[10], equ[10], "{}r_4510".format(r))
    # 6-7-9 pending-pending
    pending_pending(m, deg[6], deg[7], g[6], g[7], var[9], con[9], equ[9], "{}r_679".format(r))

    # Returns the output variable of the round
    var_out = {}
    con_out = {}
    deg_out = {}
    g_out = {}
    # Branch 0 vc variables, and their associated deg, g
    var_out[0] = var[10]
    con_out[0] = con[10]
    deg_out[0] = [deg[4], deg[5]]
    g_out[0] = [g[4], g[5]]
    # The first branch vc variable, and its associated deg, g
    var_out[1] = var[8]
    con_out[1] = con[8]
    deg_out[1] = [deg[8]]
    g_out[1] = [g[8]]
    # branch 2 vc variables, and their associated deg, g
    var_out[2] = var[9]
    con_out[2] = con[9]
    deg_out[2] = [deg[6], deg[7]]
    g_out[2] = [g[6], g[7]]

    return var_out, con_out, deg_out, g_out
def parsing_Cinimion_solution(m, R):

    for var in m.getVars():
        if "var_" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", round(var.X))
    for var in m.getVars():
        if "con_" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", round(var.X))
    for var in m.getVars():
        if "bse_" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", var.X)
                deg_name = var.VarName.replace("bse_", "deg_")
                print(deg_name, " = ", round(m.getVarByName(deg_name).X))
    for var in m.getVars():
        if "equ" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", round(var.X))
def gen_Cinimion_model(R, numVar):
    # Creating Models
    m = gp.Model("Cinimion_CICO_MILP_{}r_{}v".format(R, numVar))
    # Create input variables: pending quaternion + input state binary
    var_in = {};
    con_in = {};
    deg_in = {};
    g_in = {};
    h_in = {};
    bse_in = {};
    for i in range(3):
        var_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_0r_{}".format(i))
        con_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_0r_{}".format(i))
        deg_in[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_0r_{}".format(i))
        g_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_0r_{}".format(i))
        h_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_0r_{}".format(i))
        bse_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_0r_{}".format(i))
    m.update()

    # Input boundary constraints
    # (1/3) Input constraints in CICO
    m.addConstr(con_in[0] == 1)

    # (2/3) Input boundary add var=0,con=0 then g=0 constraints
    varPcon = {}
    for i in range(3):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_MDS_{}".format(-1, i))
    m.update()
    for i in range(3):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3）communal restraint
    for i in range(3):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # First round of constraints
    in0, in1, in2, in8, in9, in13 = Cinimion_first_round_model(m, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # Intermediate wheel constraints
    for r in range(1,R-1):
        in0, in1, in2, in8, in9, in13 = Cinimion_middle_round_model(m, r, in0, in1, in2, in8, in9, in13)

    # Final round of restraints
    var_out, con_out, deg_out, g_out = Cinimion_last_round_model(m, R-1, in0, in1, in2, in8, in9, in13)

    # Constraints on the output boundary
    # （1/3）Output constraints in CICO
    m.addConstr(con_out[0] == 1)
    #m.addConstr(con_out[1] == 1)
    # （2/3）Output bounds add the constraint that var=0,con=0 then g=0
    for i in range(3):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_{}".format(R, i))
    m.update()
    for i in range(3):
        m.addConstr(varPcon[i] == var_out[i] + con_out[i])
        for cg in g_out[i]:
            m.addConstr((varPcon[i] == 0) >> (cg == 0))
    m.update()
    # （3/3） communal restraint
    for i in range(3):
        for cdeg in deg_out[i]:
            common_constraints(m, var_out[i], con_out[i], cdeg)

    # Number of constants
    all_con_vars = [var for var in m.getVars() if "con_" in var.VarName[:4]]
    m.addConstr(sum(all_con_vars) == 3)  # 常数个数等于t，则平均有1个解

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

    #m.addConstr(m.getVarByName("con_0r_0")==1)
    #m.addConstr(m.getVarByName("con_2r_10") == 1)
    #m.addConstr(m.getVarByName("con_2r_8") == 1)

    #m.addConstr(m.getVarByName("var_0r_1")==1)
    #m.addConstr(m.getVarByName("var_0r_2") == 1)



    # Model setup
    m.setParam("TimeLimit", 600)  # Set the run to stop after 600 seconds and output the currently found optimal solution

    obj_bound = get_obj_bound(2*R+1, numVar, R, 2)
    #m.setParam("BestObjStop", obj_bound) #Set the objective function to stop when it falls below 100
    print("It's better to reach obj_bound", obj_bound)
    m.write("Cinimion_CICO_MILP_{}r_{}v.lp".format(R, numVar))
    m.optimize()
    m.write("Cinimion_MILP_{}r_{}v.sol".format(R, numVar))

    parsing_Cinimion_solution(m, R)
    print('round：{}，'.format(R))
    print('Number of variables：{}，'.format(numVar))
    print('Sum of the number of equations：{}'.format(obj.X))
    if numVar == 1:
        print('Solution complexity under univariate：2^{}'.format(univariate_comp(int(obj.X))))
    else:
        print('solution complexity：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))
    
if __name__ == '__main__':
    gen_Cinimion_model(10, 2)
