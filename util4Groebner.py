import gurobipy as gp
from gurobipy import GRB
from math import comb, log2

def get_obj_bound(exp, numVar, R, alpha):
    tempobj = alpha * R
    tempexp = log2(comb(tempobj + 1, numVar) ** 2)
    while (tempexp > exp):
        tempobj -= 1
        tempexp = log2(comb(tempobj + 1, numVar) ** 2)
    return tempobj

def univariate_comp(d):
    return log2(d * log2(d) * (log2(d) + 64) * log2(log2(d)))

    
def gen_pending_vars_dghb(m, r, i): # generate deg, g, h, bse variables for pending ends
    deg = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "deg_{}r_{}".format(r, i))
    g = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "g_{}r_{}".format(r, i))
    h = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "h_{}r_{}".format(r, i))
    bse =  m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "bse_{}r_{}".format(r, i))
    m.update()
    return deg, g, h, bse

def gen_branching_vars_vc(m, r, i): # var,con
    var = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "var_{}r_{}".format(r, i))
    con = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "con_{}r_{}".format(r, i))
    m.update()
    return var, con
def gen_vars_d(m, r, i): # deg
    deg = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "deg_{}r_{}".format(r, i))
    m.update()
    return deg
def gen_branching_vars_vcd(m, r, i): # var, con, deg
    var = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "var_{}r_{}".format(r, i))
    con = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "con_{}r_{}".format(r, i))
    deg = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "deg_{}r_{}".format(r, i))
    m.update()
    return var, con, deg

def gen_branching_vars_vce(m, r, i): #var, con, equ
    var = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "var_{}r_{}".format(r, i))
    con = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "con_{}r_{}".format(r, i))
    equ = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "equ_{}r_{}".format(r, i))
    m.update()
    return var, con, equ

def gen_branching_vars_vcde(m, r, i): # var, con, deg, equ
    var = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "var_{}r_{}".format(r, i))
    con = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "con_{}r_{}".format(r, i))
    deg = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "deg_{}r_{}".format(r, i))
    equ = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "equ_{}r_{}".format(r, i))
    m.update()
    return var, con, deg, equ


def gen_pending_vars_ghb(m, r, i): # generate g,h,bse variabels for pending ends
    g = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "g_{}r_{}".format(r, i))
    h = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "h_{}r_{}".format(r, i))
    bse =  m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "bse_{}r_{}".format(r, i))
    m.update()
    return g, h, bse
    
def MDS_module(m, t, var, con, deg, g, h, bse, equ, index): # constraints for MDS module
    # index: [round, name] index[0] is the round index; 
    assert(len(var) == 2*t)
    assert (len(con) == 2 * t)

    # add assistant variables
    #varPcon = m.addVars(2 * t, lb = 0, vtype = GRB.INTEGER, name = "varPcon_" + index[1])
    # bseTdeg = bse × deg 
    bseTdeg = m.addVars( 2 * t, lb = 0, vtype = GRB.INTEGER, name = "bseTdeg_" + index)
    # max_bseTdeg is the basis degree
    max_bseTdeg = m.addVar(lb=0, vtype=GRB.INTEGER, name="max_bseTdeg_" + index)
    # hTdeg = h × deg 
    hTdeg = m.addVars(2 * t, lb=0, vtype=GRB.INTEGER, name="hTdeg_" + index)
    # sum_con
    sum_con = m.addVar(lb=0, vtype=GRB.INTEGER, name="sum_con_" + index)
    # bin_max_bseTdeg
    bin_max_bseTdeg = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "bin_max_bseTdeg_" + index)
    m.update()

    for i in range(2 * t):
        m.addQConstr(bseTdeg[i] == bse[i] * deg[i] )
        #m.addConstr(varPcon[i] == var[i] + con[i])

    m.addGenConstrMax(max_bseTdeg, [bseTdeg[i] for i in range(2 * t)])
    m.addConstr(sum_con == sum(con))

    m.addConstr((bin_max_bseTdeg == 0) >> (max_bseTdeg >= 1))
    m.addConstr((bin_max_bseTdeg == 1) >> (max_bseTdeg == 0))

    for i in range(2 * t):
        m.addQConstr(hTdeg[i] == h[i] * deg[i])


    # add constraints

    # If a variable is introduced, the branch determines
    m.addConstrs((var[i] == 1) >> (g[i] == 1) for i in range(2*t))
    # If a constant is introduced, the branch determines
    m.addConstrs((con[i] == 1) >> (g[i] ==1) for i in range(2*t))
    # add sum_con <=t
    m.addConstr(sum_con <= t)
    #The base can only be chosen from states that consume degrees of freedom (non-essential, the value of h already guarantees this)
    #m.addConstrs(bse[i] <= g[i] for i in range(2*t))
    #h is numerically equal to g - bse
    m.addConstrs(h[i] == g[i] - bse[i] for i in range(2*t))
    # There are t branches in the base
    m.addConstr(sum(bse[i] for i in range(2 * t)) == t)
    # The degree of branches not selected into the base should be bigger than or equal to the number of base
    m.addConstrs((h[i] == 1) >> (max_bseTdeg <= deg[i]) for i in range(2*t))
    # The sum of the degrees of the added equations
    m.addLConstr(equ == sum(hTdeg[i] for i in range(2*t)))
    # if a branch is not determined, it should be expressed by basis
    m.addConstrs((g[i] == 0) >> (deg[i] == max_bseTdeg) for i in range(2*t))

    # indicator constraints condition must be a binary variable == 0 or 1
    m.addConstr((bin_max_bseTdeg == 1) >> (sum_con == t))
    m.update()
    return bseTdeg, max_bseTdeg, hTdeg


def taddition_module(m, t, var, con, deg, g, h, bse, equ, index): # Defining t-additon module constraints
    # index: name with number of rounds, position and module name
    assert(len(var) == t)

    # Adding auxiliary variables
    # bseTdeg is the product of bse and deg.
    bseTdeg = m.addVars( t, lb = 0, vtype = GRB.INTEGER, name = "bseTdeg_" + index)
    # max_bseTdeg is the number of affine layer bases
    max_bseTdeg = m.addVar(lb=0, vtype=GRB.INTEGER, name="max_bseTdeg_" + index)
    # hTdeg is the product of h and deg
    hTdeg = m.addVars(t, lb=0, vtype=GRB.INTEGER, name="hTdeg_" + index)
    # Add sum_con
    sum_con = m.addVar(lb=0, vtype=GRB.INTEGER, name="sum_con_" + index)
    # Add bin_max_bseTdeg
    bin_max_bseTdeg = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bin_max_bseTdeg_" + index)
    m.update()

    for i in range(t):
        m.addQConstr(bseTdeg[i] == bse[i] * deg[i])

    m.addGenConstrMax(max_bseTdeg, [bseTdeg[i] for i in range(t)])
    m.addConstr(sum_con == sum(con))

    m.addConstr((bin_max_bseTdeg == 0) >> (max_bseTdeg >= 1))
    m.addConstr((bin_max_bseTdeg == 1) >> (max_bseTdeg == 0))


    for i in range(t):
        m.addQConstr(hTdeg[i] == h[i] * deg[i] )

    # Adding Constraints
    # variable is introduced, the branch determines
    m.addConstrs((var[i] == 1) >> (g[i] == 1) for i in range(t))
    # constant is introduced, the branch determines
    m.addConstrs((con[i] == 1) >> (g[i] ==1) for i in range(t))
    # Add sum_con <= t
    m.addConstr(sum_con <= t - 1)
    # The base can only be selected from states that consume degrees of freedom
    #m.addConstrs(bse[i] <= g[i] for i in range(t))
    # h is numerically equal to g - bse
    m.addConstrs(h[i] == g[i] - bse[i] for i in range(t))
    # There are t-1 branches in the base
    m.addConstr(sum(bse[i] for i in range(t)) == t-1)
    # The degree of branches not selected into the base should be bigger than or equal to the number of base 
    m.addConstrs((h[i] == 1) >> (max_bseTdeg <= deg[i]) for i in range(t))
    # Sum of the degrees of equations added
    m.addConstr(equ == sum(hTdeg[i] for i in range(t)))
    # if a branch is not determined, it should be expressed by basis
    m.addConstrs((g[i] == 0) >> (deg[i] == max_bseTdeg) for i in range(t))

    # indicator constraints condition must be a binary variable == 0 or 1
    m.addConstr((bin_max_bseTdeg == 1) >> (sum_con == t-1))
    m.update()
    return bseTdeg, max_bseTdeg, hTdeg

def NLUP_module(m, var, con, deg, equ, index, alpha): # Univariate nonlinear modules
    # var: input and output variables, var[0] is upstream, var[1] is downstream. con, deg are also
    # index: name with number of rounds, position and module name
    assert(len(var) == 2)
    # Adding Constraints
    # If the constant is in the downstream position, its upstream state remains constant and vice versa
    m.addConstr((con[0] == 1) >> (con[1] == 1))
    m.addConstr((con[1] == 1) >> (con[0] == 1))
    # Introducing variables downstream introduces the equation
    m.addConstr((var[1]==1) >> (equ == alpha * deg[0]))
    # No variables are introduced downstream, then the expression from upstream
    m.addConstr((var[1]==0) >> (deg[1] == alpha * deg[0]))
    m.addConstr((var[1]==0) >> (equ == 0))
    m.update()
    
def NLMP_module(m, deg_in, var_out, con_out, deg_out, equ, index, alpha):# Nonlinear multivariate polynomial constraints
    # m: model
    # var_in, ... g_in: input variables
    # var_out, ... g_out: output variables
    # index: name with number of rounds, position and module name
    # alpha: number of polynomials
    
    # Adding auxiliary variables
    max_deg_in = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "max_deg_in" + index)
    varPcon_out = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_out_" + index)
    m.update()

    m.addGenConstrMax(max_deg_in, deg_in)

    m.addConstr(varPcon_out == var_out + con_out)

    # Adding Constraints
    
    # Introducing variables downstream adds the equation
    m.addConstr((var_out == 1) >> (equ == alpha * max_deg_in))

    # downstream is a constant, the equation is introduced
    m.addConstr((con_out == 1) >> (equ == alpha * max_deg_in))
    # The downstream does not introduce variables and is non-constant, it is expressed by the upstream    
    m.addConstr((varPcon_out == 0) >> (deg_out == alpha * max_deg_in))
    m.addConstr((varPcon_out == 0) >> (equ == 0))
    m.update()
    return max_deg_in, varPcon_out
    

def multiplication_module(m, deg_in, var_out, con_out, deg_out, equ, index):# Multiplication Module Constraints
    # m: model
    # var_in, ... g_in: input variables
    # var_out, ... g_out: output variables
    # index: name with number of rounds, position and module name
    # alpha: number of polynomials
    assert(len(deg_in) == 2)
    
    # Adding auxiliary variables
    varPcon_out = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "varPcon_out_" + index)
    m.update()
    m.addConstr(varPcon_out == var_out + con_out)

    # Adding Constraints
    
    # Introducing variables downstream adds the equation
    m.addConstr((var_out == 1) >> (equ == deg_in[0] + deg_in[1]))

    # downstream is a constant, the equation is introduced
    m.addConstr((con_out == 1) >> (equ == deg_in[0] + deg_in[1]))
    # The downstream does not introduce variables and is non-constant, it is expressed by the upstream    
    m.addConstr((varPcon_out == 0) >> (deg_out == deg_in[0] + deg_in[1]))
    m.addConstr((varPcon_out == 0) >> (equ == 0))
    m.update()
    return varPcon_out

def common_constraints(m, var, con, deg):
    m.addConstr((var == 1) >> (deg == 1))
    m.addConstr((con == 1) >> (deg == 0))
    m.update()

def down_up(m, var, con, deg):
    common_constraints(m, var, con, deg)

def down_pending(m, var, con, g, deg, index):
    # Adding auxiliary variables
    varPcon = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "varPcon_" + index)
    m.update()
    m.addConstr(varPcon == var + con)

    # Adding Constraints
    m.addConstr((varPcon == 0) >> (g==1))
    common_constraints(m, var, con, deg)
    m.update()
    return varPcon

def down_pending_pending(m, g1, deg1, g2, deg2, var, con, deg, equ, index):
    # Adding auxiliary variables
    #g1Pg2 = m.addVar(lb = 0, ub = 2, vtype = GRB.INTEGER, name = "g1Pg2_" + index)
    BIN_g1Pg2eq2 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_g1Pg2eq2_" + index)
    BIN_g1Pg2eq1 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_g1Pg2eq1_" + index)
    BIN_g1Pg2eq0 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_g1Pg2eq0_" + index)
    
    min_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="min_deg_" + index)
    max_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="max_deg_" + index)
    


    m.update()

    
    # (g1,g2, BIN_g1Pg2eq0)
    m.addLConstr(-1 * g1 - g2 - 2 * BIN_g1Pg2eq0 >= -2)
    m.addLConstr(2 * g1 + 2 * g2 + 2 * BIN_g1Pg2eq0 >= 2)

    # (g1,g2, BIN_g1Pg2eq1)
    m.addLConstr(2 * g1 + 2 * g2 - 2 * BIN_g1Pg2eq1 >= 0)
    m.addLConstr(-1 * g1 - g2 - BIN_g1Pg2eq1 >= -2)
    m.addLConstr(g1 - g2 + BIN_g1Pg2eq1 >= 0)
    m.addLConstr(-1 * g1 + g2 + BIN_g1Pg2eq1 >= 0)

    # (g1, g2, BIN_g1Pg2eq2)
    m.addLConstr(g1 + g2 - 2 * BIN_g1Pg2eq2 >= 0)
    m.addLConstr(-1 * g1 - g2 + BIN_g1Pg2eq2 >= -1)


    m.addGenConstrMin(min_deg, [deg1, deg2, deg])

    m.addGenConstrMax(max_deg, [deg1, deg2, deg])

    
    
    # Adding Constraints
    #g1 + g2 = 2 
    m.addConstr((BIN_g1Pg2eq2 == 1) >> (deg1 == deg) )
    m.addConstr((BIN_g1Pg2eq2 == 1) >> (deg2 == deg) )
    m.addConstr((BIN_g1Pg2eq2 == 1) >> (equ == 0) )
    # g1 + g2 = 0
    m.addConstr((BIN_g1Pg2eq0 == 1) >> (equ == deg1 + deg2 + deg - min_deg))
    # g1 + g2 = 1
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (deg1 + deg2 + deg - max_deg == 2 * min_deg))
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (equ == max_deg))

    # communal restraint
    common_constraints(m, var, con, deg)
    common_constraints(m, var, con, deg1)
    common_constraints(m, var, con, deg2)

    m.update()
    
def pending_up(m, var, con, g, deg, index):
    # Adding auxiliary variables
    varPcon = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "varPcon_" + index)
    m.update()
    m.addConstr(varPcon == var + con)

    # Adding Constraints
    m.addConstr((varPcon == 0) >> (g == 0))

    # communal restraint
    common_constraints(m, var, con, deg)
    m.update()
    return varPcon
    

def down_pending_ups(m, deg_d, deg_p, g, deg_us, var, con, equ, index):
    # Adding auxiliary variables
   
    min_deg_dp = m.addVar(lb=0, vtype=GRB.INTEGER, name="min_degdp_" + index)

    m.update()


 

    m.addGenConstrMin(min_deg_dp, [deg_d, deg_p])

    # Add constraints, organized by g
    m.addConstr((g ==1) >> (deg_p == deg_d))
    m.addConstr((g ==0) >> (equ == deg_d +deg_p - min_deg_dp))
   
    for deg_u in deg_us:
        m.addConstr((g ==1) >> (deg_u == deg_d))
        m.addConstr((g == 0) >> (deg_u == min_deg_dp))
        common_constraints(m, var, con, deg_u)
    common_constraints(m, var, con, deg_d)
    m.update()
    return min_deg_dp
    
def pending_pending(m, deg1, deg2, g1, g2, var, con, equ, index):
    # Adding auxiliary variables
    varPcon = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + index)
    BIN_g1Pg2eq1 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_g1Pg2eq1_" + index)
    BIN_g1Pg2eq0 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_g1Pg2eq0_" + index)

    max_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="max_deg_" + index)
    m.update()

    m.addConstr(varPcon == var + con)
   

    # (g1,g2, BIN_g1Pg2eq0)
    m.addLConstr(-1 * g1 - g2 - 2 * BIN_g1Pg2eq0 >= -2)
    m.addLConstr(2 * g1 + 2 * g2 + 2 * BIN_g1Pg2eq0 >= 2)

    # (g1,g2, BIN_g1Pg2eq1)
    m.addLConstr(2 * g1 + 2 * g2 - 2 * BIN_g1Pg2eq1 >= 0)
    m.addLConstr(-1 * g1 - g2 - BIN_g1Pg2eq1 >= -2)
    m.addLConstr(g1 - g2 + BIN_g1Pg2eq1 >= 0)
    m.addLConstr(-1 * g1 + g2 + BIN_g1Pg2eq1 >= 0)

    m.addGenConstrMax(max_deg, [deg1, deg2])
    
    # Adding Constraints
    m.addConstr((varPcon == 0) >> (g1 + g2 <= 1))
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (deg1 == deg2) )
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (equ == 0))

    m.addConstr((BIN_g1Pg2eq0 == 1) >> (equ == max_deg))

    # communal restraint
    common_constraints(m, var, con, deg1)
    common_constraints(m, var, con, deg2)
    m.update()
  




def pending_pending_up(m, deg1, deg2, g1, g2, var, con, deg, equ, index):
    # Adding auxiliary variables
    min_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="min_deg_" + index)
    m.update()



    m.addGenConstrMin(min_deg, [deg1, deg2])

    # Adding Constraints
    pending_pending(m, deg1, deg2, g1, g2, var, con, equ, index)

    m.addConstr(deg == min_deg)

    # communal restraint
    common_constraints(m, var, con, deg)

    m.update()
   




def pending_pending_pending(m, deg1, deg2, deg3, g1, g2, g3, var, con, equ, index):
    # Adding auxiliary variables
  
    varPcon = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + index)
    BIN_sumg0 = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "BIN_sumg0_" + index)
    BIN_sumg1 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_sumg1_" + index)
    BIN_sumg2 = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="BIN_sumg2_" + index)
    max_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="max_deg_" + index)
    min_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="min_deg_" + index)
   
    m.update()

    


    # (g1,g2,g3, BIN_sumg0)
    m.addLConstr(-g1 - g3 - 2 * BIN_sumg0 >= -2)
    m.addLConstr(-g1 - g2 - 2 * BIN_sumg0 >= -2)
    m.addLConstr(g1 + g2 + g3 + BIN_sumg0 >= 1)


    # (g1,g2,g3, BIN_sumg1)
    m.addLConstr(-1 * g1 + g2 + g3 + 2 * BIN_sumg1 >= 0)
    m.addLConstr(g1 - g2 + g3 + 2* BIN_sumg2 >= 0)
    m.addLConstr(2 * g1 + 2 * g2+ 2 * g3 - 2 * BIN_sumg1 >= 0 )
    m.addLConstr(g1 + g2 - 2 * g3 + BIN_sumg1 >= -1)
    m.addLConstr(- 2 * g1 - 2* g2 - 2 * g3 - 3 * BIN_sumg1 >= -6)

    # (g1,g2,g3, BIN_sumg2)
    m.addLConstr(2 * g1 + g2 + 2 * g3 - 3 * BIN_sumg2 >= 0)
    m.addLConstr(- 1 * g1 - g2 + 2* g3 + BIN_sumg2 >= -1)
    m.addLConstr(g1 - g2 - g3 + 2 * BIN_sumg2 >= -1)
    m.addLConstr(-1 * g1 + g2 - 2 * g3 + BIN_sumg2 >= -2)
    m.addLConstr(-1 * g1 - 2 * g2 - g3 - BIN_sumg2 >= -4)

    m.addConstr(varPcon == var + con)

    m.addGenConstrMax(max_deg, [deg1, deg2, deg3])

    m.addGenConstrMin(min_deg, [deg1, deg2, deg3])

    

    # Adding Constraints
    m.addConstr((varPcon == 0) >> (g1 + g2 + g3 <= 2))
    # g1 + g2 + g3 = 2.
    m.addConstr((BIN_sumg2 == 1) >> (deg1 == deg2))
    m.addConstr((BIN_sumg2 == 1) >> (deg1 == deg3))
    m.addConstr((BIN_sumg2 == 1) >> (equ == 0))

    # g1 + g2 + g3 = 1.
    m.addConstr((BIN_sumg1 == 1) >> (deg1 + deg2 + deg3 -max_deg - 2* min_deg == 0))
    m.addConstr((BIN_sumg1 == 1) >> (equ == max_deg))

    # g1 + g2 + g3 = 0.
    m.addConstr((BIN_sumg0 == 1) >> (equ == deg1 + deg2 + deg3 - min_deg))

    # communal restraint
    common_constraints(m, var, con, deg1)
    common_constraints(m, var, con, deg2)
    common_constraints(m, var, con, deg3)

    m.update()
    return min_deg

def pending_pending_pending_up(m, deg1, deg2, deg3, g1, g2, g3, var, con, equ, deg, index):

    min_deg = pending_pending_pending(m, deg1, deg2, deg3, g1, g2, g3, var, con, equ, index)

    #The upstream count is the smallest of the three pending
    m.addLConstr(deg == min_deg)

    # communal restraint
    common_constraints(m, var, con, deg)

    m.update()

