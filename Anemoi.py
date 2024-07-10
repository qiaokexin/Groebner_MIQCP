from util4Groebner import *
from math import comb, log2

alpha = 3
def Anemoi_round_model(m, r, l, var_in, con_in, deg_in, g_in, h_in, bse_in):
    # Modeling the rth round of the round function, l is the width of the MDS, var_in and the variables that follow are the input variables for this round

    # Initializing the auxiliary variable dictionary
    var={}; con ={}; deg = {}; g = {}; h = {}; bse = {}; equ = {};
    # Variable Settings
    for i in range(l):
        # Add variables to each branch module
        # The 0th point position is the pending quaternion + base binary, obtained directly from the input
        var[26 * i + 0] = var_in[i]
        con[26 * i + 0] = con_in[i]
        deg[26 * i + 0] = deg_in[i]
        g[26 * i + 0] = g_in[i]
        h[26 * i + 0] = h_in[i]
        bse[26 * i + 0] = bse_in[i]
        # The 1st point position is the pending quaternion + the base binary, obtained directly from the input
        var[26 * i + 1] = var_in[l + i]
        con[26 * i + 1] = con_in[l + i]
        deg[26 * i + 1] = deg_in[l + i]
        g[26 * i + 1] = g_in[l + i]
        h[26 * i + 1] = h_in[l + i]
        bse[26 * i + 1] = bse_in[l + i]
        # The second point is the pending quaternion
        deg[26 * i + 2], g[26 * i + 2], h[26 * i + 2], bse[26 * i + 2] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,2))
        # The third point is the pending quaternion
        deg[26 * i + 3], g[26 * i + 3], h[26 * i + 3], bse[26 * i + 3] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,3))
        # Points 4 and 5 are the var,con,equ triples.
        var[26 * i + 4], con[26 * i + 4], equ[26 * i + 4] = gen_branching_vars_vce(m, r, "{}th_{}".format(i,4))
        var[26 * i + 5], con[26 * i + 5], equ[26 * i + 5] = gen_branching_vars_vce(m, r, "{}th_{}".format(i,5))
        # Points 6,7,8,9 are pendent quaternions.
        for j in range(6, 10):
            deg[26 * i + j], g[26 * i + j], h[26 * i + j], bse[26 * i + j] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,j))
        # The 10th point is the var,con,equ triad.
        var[26 * i + 10], con[26 * i + 10], equ[26 * i + 10] = gen_branching_vars_vce(m, r, "{}th_{}".format(i,10))
        # The 11th point is the pending quaternion
        deg[26 * i + 11], g[26 * i + 11], h[26 * i + 11], bse[26 * i + 11] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,11))
        # The 12th point is the var,con,deg,equ quaternion.
        var[26 * i + 12], con[26 * i + 12], deg[26 * i + 12], equ[26 * i + 12] = gen_branching_vars_vcde(m, r, "{}th_{}".format(i,12))
        # The 13th point is the pending quaternion + binary
        deg[26 * i + 13], g[26 * i + 13], h[26 * i + 13], bse[26 * i + 13] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,13))
        var[26 * i + 13], con[26 * i + 13] = gen_branching_vars_vc(m, r, "{}th_{}".format(i,13))
        # Points 14, 15 are pendant quaternions.
        deg[26 * i + 14], g[26 * i + 14], h[26 * i + 14], bse[26 * i + 14] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,14))
        deg[26 * i + 15], g[26 * i + 15], h[26 * i + 15], bse[26 * i + 15] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,15))
        # The 16th point is the pending quaternion + binary
        deg[26 * i + 16], g[26 * i + 16], h[26 * i + 16], bse[26 * i + 16] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,16))
        var[26 * i + 16], con[26 * i + 16] = gen_branching_vars_vc(m, r, "{}th_{}".format(i,16))
        # The 17th point is the var,con,deg,equ quaternion.
        var[26 * i + 17], con[26 * i + 17], deg[26 * i + 17], equ[26 * i + 17] = gen_branching_vars_vcde(m, r, "{}th_{}".format(i,17))
        # The 18th and 19th points are pendant quaternions.
        deg[26 * i + 18], g[26 * i + 18], h[26 * i + 18], bse[26 * i + 18] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,18))
        deg[26 * i + 19], g[26 * i + 19], h[26 * i + 19], bse[26 * i + 19] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,19))
        # The 20th point is the var,con,deg,equ quaternion, but where var, con are named according to the next round of indicators
        deg[26 * i + 20] = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "deg_{}r_{}th_{}".format(r, i, 20))
        equ[26 * i + 20] = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "equ_{}r_{}th_{}".format(r, i, 20))
        var[26 * i + 20], con[26 * i + 20] = gen_branching_vars_vc(m, r + 1, "{}th_{}".format(i, 1))
        # The 21st point is pending quaternion + binary
        deg[26 * i + 21], g[26 * i + 21], h[26 * i + 21], bse[26 * i + 21] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,21))
        var[26 * i + 21], con[26 * i + 21] = gen_branching_vars_vc(m, r, "{}th_{}".format(i,21))
        # The 22nd point is the pending quaternion.
        deg[26 * i + 22], g[26 * i + 22], h[26 * i + 22], bse[26 * i + 22] = gen_pending_vars_dghb(m, r,
                                                                                                   "{}th_{}".format(i,
                                                                                                                    22))
        # The 23rd point is the equ+vc binary, where the vc binary is named according to the next round of indicators
        equ[26 * i + 23] = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ_{}r_{}th_{}".format(r, i, 23))
        var[26 * i + 23], con[26 * i + 23] = gen_branching_vars_vc(m, r + 1, "{}th_{}".format(i, 0))


        # Points 24 and 25 are the inputs for the next round. Pending quaternions, named according to the next round's corners.
        deg[26 * i + 24], g[26 * i + 24], h[26 * i + 24], bse[26 * i + 24] = gen_pending_vars_dghb(m, r+1, "{}th_{}".format(i,0))
        deg[26 * i + 25], g[26 * i + 25], h[26 * i + 25], bse[26 * i + 25] = gen_pending_vars_dghb(m, r + 1,
                                                                                                   "{}th_{}".format(i,1))


        # Add the equ variable to each arithmetic module
        # Modules of operation
        equ_smallMDS= m.addVar(lb=0, vtype=GRB.INTEGER, name="equ_smallMDS_{}r_{}".format(r, i))
        equ111314 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ111314_{}r_{}".format(r, i))
        equ1213 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1213_{}r_{}".format(r, i))
        equ1617 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1617_{}r_{}".format(r, i))
        equ151619 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ151619_{}r_{}".format(r, i))
        equ2021 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ2021_{}r_{}".format(r, i))
        equ182122 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ182122_{}r_{}".format(r, i))

        # Adding constraints to each branch module
        # 2,4,6 pending-pending
        pending_pending(m, deg[26 * i + 2], deg[26 * i + 6], g[26 * i + 2], g[26 * i + 6], var[26 * i + 4],
                        con[26 * i + 4], equ[26 * i + 4], "{}r_{}th_2-4-6".format(r, i))
        # 3,5,7 pending-pending
        pending_pending(m, deg[26 * i + 3], deg[26 * i + 7], g[26 * i + 3], g[26 * i + 7], var[26 * i + 5],
                        con[26 * i + 5], equ[26 * i + 5], "{}r_{}th_3-5-7".format(r, i))
        # 8,10,11 pending-pending
        pending_pending(m, deg[26 * i + 8], deg[26 * i + 11], g[26 * i + 8], g[26 * i + 11], var[26 * i + 10],
                        con[26 * i + 10], equ[26 * i + 10], "{}r_{}th_8-10-11".format(r, i))
        # 9,12,15 pending-pending-up
        pending_pending_up(m, deg[26 * i + 9], deg[26 * i + 15], g[26 * i + 9], g[26 * i + 15], var[26 * i + 12],
                           con[26 * i + 12], deg[26 * i + 12], equ[26 * i + 12], "{}r_{}th_9-12-15".format(r, i))
        # 13 down-pending
        down_pending(m, var[26 * i + 13], con[26 * i + 13], g[26 * i + 13], deg[26 * i + 13],
                     "{}r_{}th_13".format(r, i))
        # 14,17,18 down-pending-pending
        down_pending_pending(m, g[26 * i + 14], deg[26 * i + 14], g[26 * i + 18], deg[26 * i + 18], var[26 * i + 17],
                             con[26 * i + 17], deg[26 * i + 17], equ[26 * i + 17], "{}r_{}th_14-17-18".format(r, i))
        # 16 pengding-up
        pending_up(m, var[26 * i + 16], con[26 * i + 16], g[26 * i + 16], deg[26 * i + 16], "{}r_{}th_16".format(r, i))
        # 19,20,25 pending-pending-up
        pending_pending_up(m, deg[26 * i + 19], deg[26 * i + 25], g[26 * i + 19], g[26 * i + 25],
                           var[26 * i + 20], con[26 * i + 20], deg[26 * i + 20], equ[26 * i + 20], "{}r_{}th_19-20-25".format(r, i))

        # 21 down-pending
        down_pending(m, var[26 * i + 21], con[26 * i + 21], g[26 * i + 21], deg[26 * i + 21],
                     "{}r_{}th_21".format(r, i))
        # 22,23,24 pending-pending
        pending_pending(m, deg[26 * i + 22], deg[26 * i + 24], g[26 * i + 22], g[26 * i + 24], var[26 * i + 23],
                        con[26 * i + 23], equ[26 * i + 23], "{}r_{}th_22-23-24".format(r, i))


        # Adding constraints to each branch module
        # 6,7,8,9 小MDS
        MDS_module(m, 2, [var[26 * i + 4], var[26 * i + 5], var[26 * i + 10], var[26 * i + 12]],\
                   [con[26 * i + 4], con[26 * i + 5], con[26 * i + 10], con[26 * i + 12]],
                   [deg[26 * i + 6], deg[26 * i + 7], deg[26 * i + 8], deg[26 * i + 9]],\
                   [g[26 * i + 6], g[26 * i + 7], g[26 * i + 8], g[26 * i + 9]],
                   [h[26 * i + 6], h[26 * i + 7], h[26 * i + 8], h[26 * i + 9]],
                   [bse[26 * i + 6], bse[26 * i + 7], bse[26 * i + 8], bse[26 * i + 9]], \
                   equ_smallMDS, "_{}r_{}th".format(r, i))
        # 11,13,14 taddition
        taddition_module(m, 3, [var[26 * i + 10], var[26 * i +13], var[26 * i + 17]],
                         [con[26 * i + 10], con[26 * i + 13], con[26 * i + 17]],
                         [deg[26 * i + 11], deg[26 * i +13], deg[26 * i + 14]],
                         [g[26 * i + 11], g[26 * i + 13], g[26 * i + 14]],
                         [h[26 * i + 11], h[26 * i + 13], h[26 * i + 14]],
                         [bse[26 * i + 11], bse[26 * i + 13], bse[26 * i + 14]],
                         equ111314, "{}r_{}th_111314".format(r, i))
        # 12,13 Sbox
        NLUP_module(m, [var[26 * i + 12], var[26 * i + 13]], [con[26 * i + 12], con[26 * i + 13]],
                    [deg[26 * i + 12], deg[26 * i + 13]], equ1213, "{}r_{}th_1213".format(r, i), 2)
        # 15,16,19 taddition
        taddition_module(m, 3, [var[26 * i + 12], var[26 * i + 16], var[26 * i + 20]],
                         [con[26 * i + 12], con[26 * i + 16], con[26 * i + 20]],
                         [deg[26 * i + 15], deg[26 * i + 16], deg[26 * i + 19]],
                         [g[26 * i + 15], g[26 * i + 16], g[26 * i + 19]],
                         [h[26 * i + 15], h[26 * i + 16], h[26 * i + 19]],
                         [bse[26 * i + 15], bse[26 * i + 16], bse[26 * i + 19]],
                         equ151619, "{}r_{}th_151619".format(r, i))
        # 16, 17 Sbox
        NLUP_module(m, [var[26 * i + 16], var[26 * i + 17]], [con[26 * i + 16], con[26 * i + 17]],
                    [deg[26 * i + 16], deg[26 * i + 17]], equ1617, "{}r_{}th_1617".format(r, i), alpha)
        # 18,21,22 tadditon
        taddition_module(m, 3, [var[26 * i + 17], var[26 * i + 21], var[26 * i + 23]],
                         [con[26 * i + 17], con[26 * i + 21], con[26 * i + 23]],
                         [deg[26 * i + 18], deg[26 * i + 21], deg[26 * i + 22]],
                         [g[26 * i + 18], g[26 * i + 21], g[26 * i + 22]],
                         [h[26 * i + 18], h[26 * i + 21], h[26 * i + 22]],
                         [bse[26 * i + 18], bse[26 * i + 21], bse[26 * i + 22]],
                         equ182122, "{}r_{}th_182122".format(r, i))
        # 20,21 Sbox
        NLUP_module(m, [var[26 * i + 20], var[26 * i + 21]], [con[26 * i + 20], con[26 * i + 21]],
                    [deg[26 * i + 20], deg[26 * i + 21]], equ2021, "{}r_{}th_2021".format(r, i), 2)

    # If l is greater than or equal to 2, add the affine layer
    if (l >= 2):
        equMDS1 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS1_{}r_{}".format(r, i))
        equMDS2 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS2_{}r_{}".format(r, i))

        # Prepare the input and output variables for the MDS module
        MDS_var1 = []
        MDS_con1 = []
        MDS_deg1 = []
        MDS_g1 = []
        MDS_h1 = []
        MDS_bse1 = []
        MDS_var2 = []
        MDS_con2 = []
        MDS_deg2 = []
        MDS_g2 = []
        MDS_h2 = []
        MDS_bse2 = []
        for i in range(l):
            # Variables for the first MDS
            MDS_var1.append(var[26 * i])
            MDS_var1.append(var[26 * i + 4])
            MDS_con1.append(con[26 * i])
            MDS_con1.append(con[26 * i + 4])
            MDS_deg1.append(deg[26 * i])
            MDS_deg1.append(deg[26 * i + 2])
            MDS_g1.append(g[26 * i])
            MDS_g1.append(g[26 * i + 2])
            MDS_h1.append(h[26 * i])
            MDS_h1.append(h[26 * i + 2])
            MDS_bse1.append(bse[26 * i])
            MDS_bse1.append(bse[26 * i + 2])
            # Variables for the second MDS
            MDS_var2.append(var[26 * i + 1])
            MDS_var2.append(var[26 * i + 5])
            MDS_con2.append(con[26 * i + 1])
            MDS_con2.append(con[26 * i + 5])
            MDS_deg2.append(deg[26 * i + 1])
            MDS_deg2.append(deg[26 * i + 3])
            MDS_g2.append(g[26 * i + 1])
            MDS_g2.append(g[26 * i + 3])
            MDS_h2.append(h[26 * i + 1])
            MDS_h2.append(h[26 * i + 3])
            MDS_bse2.append(bse[26 * i + 1])
            MDS_bse2.append(bse[26 * i + 3])
            # The first MDS constraint
        MDS_module(m, l, MDS_var1, MDS_con1, MDS_deg1, MDS_g1, MDS_h1, MDS_bse1, equMDS1, "_{}r_MDS1".format(r))
        # The second MDS constraint
        MDS_module(m, l, MDS_var2, MDS_con2, MDS_deg2, MDS_g2, MDS_h2, MDS_bse2, equMDS2, "_{}r_MDS2".format(r))

    # Returns the output variable of the round
    var_out = []
    con_out = []
    deg_out = []
    g_out = []
    h_out = []
    bse_out = []

    for i in range(l):
        var_out.append(var[26 * i + 23])
        con_out.append(con[26 * i + 23])
        deg_out.append(deg[26 * i + 24])
        g_out.append(g[26 * i + 24])
        h_out.append(h[26 * i + 24])
        bse_out.append(bse[26 * i + 24])

    for i in range(l):
        var_out.append(var[26 * i + 20])
        con_out.append(con[26 * i + 20])
        deg_out.append(deg[26 * i + 25])
        g_out.append(g[26 * i + 25])
        h_out.append(h[26 * i + 25])
        bse_out.append(bse[26 * i + 25])

    return var_out, con_out, deg_out, g_out, h_out, bse_out
def last_MDS_model(m,  l, var_in, con_in, deg_in, g_in, h_in, bse_in):
    var_out = []
    con_out = []
    deg_out = []
    g_out = []
    h_out = []
    bse_out = []
    # Adding Output Variables
    for i in range(2*l):
        var_out.append(m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_R_{}".format(i)))
        con_out.append(m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_R_{}".format(i)))
        deg_out.append(m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_R_{}".format(i)))
        g_out.append(m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_R_{}".format(i)))
        h_out.append(m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_R_{}".format(i)))
        bse_out.append(m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_R_{}".format(i)))

    equ_M_0 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS_R_0")
    equ_M_1 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS_R_1")
    m.update()

    # MDS constraints
    MDS_module(m, l, var_in[:l] + var_out[:l], \
               con_in[:l] + con_out[:l], \
               deg_in[:l] + deg_out[:l], g_in[:l] + g_out[:l], h_in[:l] + h_out[:l], \
               bse_in[:l] + bse_out[:l], equ_M_0, "lastMDS_0" )
    MDS_module(m, l, var_in[l:] + var_out[l:], \
               con_in[l:] + con_out[l:], \
               deg_in[l:] + deg_out[l:], g_in[l:] + g_out[l:], h_in[l:] + h_out[l:], \
               bse_in[l:] + bse_out[l:], equ_M_1, "lastMDS_1")

    return var_out, con_out, deg_out, g_out, h_out, bse_out
def gen_Anemoi_model(R, l, numVar):
    # R is the number of rounds, l is half the number of branches, numVar is the number of set variables

    # Creating Models
    m = gp.Model("Anemoi_CICO_MILP_{}r_{}l_{}v".format(R, l, numVar))
    # Create input variables: pending quaternion + input state binary
    var_in = {};
    con_in = {};
    deg_in = {};
    g_in = {};
    h_in = {};
    bse_in = {};
    for i in range(l):
        var_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_0r_{}th_0".format(i))
        con_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_0r_{}th_0".format(i))
        deg_in[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_0r_{}th_0".format(i))
        g_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_0r_{}th_0".format(i))
        h_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_0r_{}th_0".format(i))
        bse_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_0r_{}th_0".format(i))
    for i in range(l):
        var_in[l+i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_0r_{}th_1".format(i))
        con_in[l+i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_0r_{}th_1".format(i))
        deg_in[l+i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_0r_{}th_1".format(i))
        g_in[l+i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_0r_{}th_1".format(i))
        h_in[l+i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_0r_{}th_1".format(i))
        bse_in[l+i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_0r_{}th_1".format(i))
    m.update()

    # Input boundary constraints
    # （1/3） Input constraints in CICO
    m.addConstr(con_in[0] == 1)

    m.update()

    # （2/3）Input boundaries add the constraint that var=0,con=0 then g=0
    varPcon = {}
    for i in range(2*l):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_in_{}".format(i))
    m.update()
    for i in range(2*l):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） communal restraint
    for i in range(2*l):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # Modeling per round
    for r in range(R):
        var_in, con_in, deg_in, g_in, h_in, bse_in = Anemoi_round_model(m, r, l, var_in, con_in, deg_in, g_in, h_in,
                                                                        bse_in)
    # The Last MDS
    var_in, con_in, deg_in, g_in, h_in, bse_in = last_MDS_model(m, l, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # Constraints on the output boundary
    # （1/3）Output constraints in CICO
    m.addConstr(con_in[0] == 1)

    m.update()

    # （2/3）Output bounds add the constraint that var=0,con=0 then g=0
    varPcon = {}
    for i in range(2 * l):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_out_{}".format(i))
    m.update()
    for i in range(2 * l):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） communal restraint
    for i in range(2 * l):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # Artificially added constraints
    '''
    m.addConstr(m.getVarByName("var_0r_0th_1")==1)
    m.addConstr(m.getVarByName("var_0r_1th_1") == 1)
    m.addConstr(m.getVarByName("var_0r_0th_16") == 1)
    m.addConstr(m.getVarByName("var_0r_1th_16") == 1)
    m.addConstr(m.getVarByName("var_1r_0th_16") == 1)
    m.addConstr(m.getVarByName("var_2r_0th_16") == 1)
    m.addConstr(m.getVarByName("var_2r_1th_16") == 1)
    '''
    m.update()

    # Number of constants
    all_con_vars = [var for var in m.getVars() if "con_" in var.VarName[:4]]
    repeated_con_vars = []  # Take only one constant before and after the S-box and remove the duplicates
    for r in range(R):
        for i in range(l):
            repeated_con_vars.append(m.getVarByName("con_{}r_{}th_{}".format(r, i, 13)))
            repeated_con_vars.append(m.getVarByName("con_{}r_{}th_{}".format(r, i, 17)))
            repeated_con_vars.append(m.getVarByName("con_{}r_{}th_{}".format(r, i, 21)))
    sum_cons = sum(all_con_vars) - sum(repeated_con_vars)
    m.addConstr(sum_cons == 2*l)  # The number of constants is equal to 2*l, then on average there is 1 solution

    # Given the number of variables
    all_vars = [var for var in m.getVars() if "var_" in var.VarName]
    m.addConstr(sum(all_vars) == numVar)

    # objective function
    # Sum of the number of equations
    all_equs = [var for var in m.getVars() if "equ" in var.VarName]
    obj = m.addVar(vtype=GRB.INTEGER, name="obj")
    m.addConstr(obj == sum(all_equs))
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()

    # Model setup
    #m.setParam("TimeLimit", 200)# Set the run to stop after 200 seconds and output the currently found optimal solution
    #m.setParam("BestObjStop", 32) #Set the objective function to stop when it falls below 100
    m.write("Anemoi_CICO_MILP_{}r_{}l_{}v.lp".format(R, l, numVar))
    m.optimize()
    m.write("Anemoi_CICO_MILP_{}r_{}l_{}v.sol".format(R, l, numVar))
    parsing_Anemoi_solution(m, R, l)
    print("Number of branches：{}，".format(2*l))
    print('round：{}，'.format(R))
    print('Number of variables：{}，'.format(numVar))
    print('Sum of the number of equations：{}'.format(obj.X))
    print('solution complexity：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))

def print_out_vector(m, l, vname):
    print(vname+"_R", end = " ")
    for j in range(2*l):
        current_var = m.getVarByName(vname+"_R_{}".format(j))
        print(round(current_var.X), end = " ")
    print()
def parsing_Anemoi_solution(m, R, l):
    # input variable
    for r in range(R):
        print('----------r = {}----------'.format(r))
        for i in range(l):
            for j in [0,1,4,5,10,12,13,16,17,21]:
                for vname in ['var', 'con']:
                    current_var = m.getVarByName(vname + '_{}r_{}th_{}'.format(r, i, j))
                    if (round(current_var.X) !=0):
                        print(current_var.VarName, ' = ',round(current_var.X))

            for j in [0,1,2,3,6,7,8,9,11,13,14,15,16,18,19,21,22]:
                for vname in ['deg', 'g', 'h', 'bse']:
                    current_var = m.getVarByName(vname+'_{}r_{}th_{}'.format(r, i, j))
                    print(current_var.VarName, ' = ', round(current_var.X))
            for j in [12,17,20]:
                current_var = m.getVarByName('deg_{}r_{}th_{}'.format(r, i, j))
                print(current_var.VarName, ' = ', round(current_var.X))
    print('----------R = {}----------'.format(R))

    for vname in ['var', 'con','deg', 'g', 'h', 'bse']:
        for j in [0, 1]:
            for i in range(l):
                current_var = m.getVarByName(vname + '_{}r_{}th_{}'.format(R, i, j))
                print(current_var.VarName, ' = ', round(current_var.X))

    # exports
    print('----------output----------'.format(R))
    print_out_vector(m, l, "var")
    print_out_vector(m, l, "con")
    print_out_vector(m, l, "deg")
    print_out_vector(m, l, "g")
    print_out_vector(m, l, "h")
    print_out_vector(m, l, "bse")
    for var in m.getVars():
        if (var.VarName[:3] == "equ"):
            if var.X > 0:
                print(var.VarName, "=", var.X)

def print_vc_type(m, r, l, vname, in1, in2, out1, out2):
    # var
    if r == 0:
        print(vname+'_0r', end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_0r_{}".format(j)).X), end=" ")
        print(end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_0r_{}".format(l + j)).X), end=" ")
    else:
        print(vname+'_{}r'.format(r), end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r - 1, j, in1)).X), end=" ")
        print(end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r - 1, j, in2)).X), end=" ")
    print()
    print(vname+'_{}r'.format(r), end=" ")
    for j in range(l):
        print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, out1)).X), end=" ")
    print(end=" ")
    for j in range(l):
        print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, out2)).X), end=" ")
    print()

def print_dghb_type(m, r, l, vname, in1, in2, out1, out2):
    # input
    if r == 0:
        print(vname+'_0r', end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_0r_{}".format(j)).X), end=" ")
        print(end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_0r_{}".format(l + j)).X), end=" ")
    else:
        print(vname+'_{}r'.format(r), end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, in1)).X), end=" ")
        print(end=" ")
        for j in range(l):
            print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, in2)).X), end=" ")
    print()
    # exports
    print(vname+'_{}r'.format(r), end=" ")
    for j in range(l):
        print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, out1)).X), end=" ")
    print(end=" ")
    for j in range(l):
        print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, out2)).X), end=" ")
    print()
def parsing_Anemoi_solution_v0(m, R, l):
    # each round
    for r in range(R):
        print('--------------r = {}-----------------'.format(r))
        print('-----r = {}-Around the large MDS----------------'.format(r))
        print_vc_type(m, r, l, 'var', 23, 20, 4, 5)
        print_vc_type(m, r, l, 'con', 23, 20, 4,5)
        print_dghb_type(m, r, l, 'deg', 0, 1, 2,3)
        print_dghb_type(m, r, l, 'g', 0, 1, 2, 3)
        print_dghb_type(m, r, l, 'h', 0, 1, 2, 3)
        print_dghb_type(m, r, l, 'bse', 0, 1, 2, 3)
    for var in m.getVars():
        if var.VarName[:3] in ["var", "con", "equ"]:
            if var.X > 0:
                print(var.VarName, "=", var.X)


if __name__ == '__main__':
    gen_Anemoi_model(2, 2, 4)
