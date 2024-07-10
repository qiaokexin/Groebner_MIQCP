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

    
def gen_pending_vars_dghb(m, r, i): # generate variable quadruple for pending ends
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

    # 引入了变量，则该分支确定
    m.addConstrs((var[i] == 1) >> (g[i] == 1) for i in range(2*t))
    # 引入了常数，则该分支确定
    m.addConstrs((con[i] == 1) >> (g[i] ==1) for i in range(2*t))
    # 添加sum_con <=t
    m.addConstr(sum_con <= t)
    # 基只能从消耗自由度的状态中选取(没必要，h的取值已经保证了这一点)
    #m.addConstrs(bse[i] <= g[i] for i in range(2*t))
    # h 在数值上等于 g - bse
    m.addConstrs(h[i] == g[i] - bse[i] for i in range(2*t))
    # 基中有t个分支
    m.addConstr(sum(bse[i] for i in range(2 * t)) == t)
    # 未选入基的分支的次数应该大于等于基的次数
    m.addConstrs((h[i] == 1) >> (max_bseTdeg <= deg[i]) for i in range(2*t))
    # 添加的方程次数之和
    m.addLConstr(equ == sum(hTdeg[i] for i in range(2*t)))
    # if a branch is not determined, it should be expressed by basis
    m.addConstrs((g[i] == 0) >> (deg[i] == max_bseTdeg) for i in range(2*t))

    # indicator constraints 条件必须是一个binary变量==0或1
    m.addConstr((bin_max_bseTdeg == 1) >> (sum_con == t))
    m.update()
    return bseTdeg, max_bseTdeg, hTdeg


def taddition_module(m, t, var, con, deg, g, h, bse, equ, index): # 定义t-additon模块约束
    # index: 名字，包含轮数、位置和模块名称
    assert(len(var) == t)

    # 添加辅助变量
    # bseTdeg 为 bse 与 deg 的乘积
    bseTdeg = m.addVars( t, lb = 0, vtype = GRB.INTEGER, name = "bseTdeg_" + index)
    # max_bseTdeg 为仿射层基的次数
    max_bseTdeg = m.addVar(lb=0, vtype=GRB.INTEGER, name="max_bseTdeg_" + index)
    # hTdeg 为 h 与 deg 的乘积
    hTdeg = m.addVars(t, lb=0, vtype=GRB.INTEGER, name="hTdeg_" + index)
    # 添加sum_con
    sum_con = m.addVar(lb=0, vtype=GRB.INTEGER, name="sum_con_" + index)
    # 添加bin_max_bseTdeg
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

    # 添加约束
    # 引入了变量，则该分支确定
    m.addConstrs((var[i] == 1) >> (g[i] == 1) for i in range(t))
    # 引入了常数，则该分支确定
    m.addConstrs((con[i] == 1) >> (g[i] ==1) for i in range(t))
    # 添加sum_con <=t
    m.addConstr(sum_con <= t - 1)
    # 基只能从消耗自由度的状态中选取
    #m.addConstrs(bse[i] <= g[i] for i in range(t))
    # h 在数值上等于 g - bse
    m.addConstrs(h[i] == g[i] - bse[i] for i in range(t))
    # 基中有t-1个分支
    m.addConstr(sum(bse[i] for i in range(t)) == t-1)
    # 未选入基的分支的次数应该大于等于基的次数 
    m.addConstrs((h[i] == 1) >> (max_bseTdeg <= deg[i]) for i in range(t))
    # 添加的方程次数之和
    m.addConstr(equ == sum(hTdeg[i] for i in range(t)))
    # if a branch is not determined, it should be expressed by basis
    m.addConstrs((g[i] == 0) >> (deg[i] == max_bseTdeg) for i in range(t))

    # indicator constraints 条件必须是一个binary变量==0或1
    m.addConstr((bin_max_bseTdeg == 1) >> (sum_con == t-1))
    m.update()
    return bseTdeg, max_bseTdeg, hTdeg

def NLUP_module(m, var, con, deg, equ, index, alpha): # 单变量非线性模块
    # var: 输入输出变量, var[0]是上游，var[1]是下游。 con, deg也是
    # index: 名字，包含轮数、位置和模块名称
    assert(len(var) == 2)
    # 添加约束
    # 若常量在下游位置，则其上游状态仍为常量，反之亦然
    m.addConstr((con[0] == 1) >> (con[1] == 1))
    m.addConstr((con[1] == 1) >> (con[0] == 1))
    # 下游引入变量，则引入方程
    m.addConstr((var[1]==1) >> (equ == alpha * deg[0]))
    # 下游不引入变量，则由上游表达
    m.addConstr((var[1]==0) >> (deg[1] == alpha * deg[0]))
    m.addConstr((var[1]==0) >> (equ == 0))
    m.update()
    
def NLMP_module(m, deg_in, var_out, con_out, deg_out, equ, index, alpha):# 非线性多变元多项式约束
    # m: model
    # var_in, ... g_in: input variables
    # var_out, ... g_out: output variables
    # index：名字，包含轮数、位置和模块名称
    # alpha: 多项式次数
    
    # 添加辅助变量
    max_deg_in = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "max_deg_in" + index)
    varPcon_out = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_out_" + index)
    m.update()

    m.addGenConstrMax(max_deg_in, deg_in)

    m.addConstr(varPcon_out == var_out + con_out)

    # 添加约束
    
    # 下游引入变量，则添加方程
    m.addConstr((var_out == 1) >> (equ == alpha * max_deg_in))

    # 下游为常数，则引入方程
    m.addConstr((con_out == 1) >> (equ == alpha * max_deg_in))
    # 下游不引入变量且非常数，则由上游表达    
    m.addConstr((varPcon_out == 0) >> (deg_out == alpha * max_deg_in))
    m.addConstr((varPcon_out == 0) >> (equ == 0))
    m.update()
    return max_deg_in, varPcon_out
    

def multiplication_module(m, deg_in, var_out, con_out, deg_out, equ, index):# 乘法模块约束
    # m: model
    # var_in, ... g_in: input variables
    # var_out, ... g_out: output variables
    # index：名字，包含轮数、位置和模块名称
    # alpha: 多项式次数
    assert(len(deg_in) == 2)
    
    # 添加辅助变量
    varPcon_out = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "varPcon_out_" + index)
    m.update()
    m.addConstr(varPcon_out == var_out + con_out)

    # 添加约束
    
    # 下游引入变量，则添加方程
    m.addConstr((var_out == 1) >> (equ == deg_in[0] + deg_in[1]))

    # 下游为常数，则引入方程
    m.addConstr((con_out == 1) >> (equ == deg_in[0] + deg_in[1]))
    # 下游不引入变量且非常数，则由上游表达    
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
    # 添加辅助变量
    varPcon = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "varPcon_" + index)
    m.update()
    m.addConstr(varPcon == var + con)

    # 添加约束
    m.addConstr((varPcon == 0) >> (g==1))
    common_constraints(m, var, con, deg)
    m.update()
    return varPcon

def down_pending_pending(m, g1, deg1, g2, deg2, var, con, deg, equ, index):
    # 添加辅助变量
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

    
    
    # 添加约束
    #g1 + g2 = 2 
    m.addConstr((BIN_g1Pg2eq2 == 1) >> (deg1 == deg) )
    m.addConstr((BIN_g1Pg2eq2 == 1) >> (deg2 == deg) )
    m.addConstr((BIN_g1Pg2eq2 == 1) >> (equ == 0) )
    # g1 + g2 = 0
    m.addConstr((BIN_g1Pg2eq0 == 1) >> (equ == deg1 + deg2 + deg - min_deg))
    # g1 + g2 = 1
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (deg1 + deg2 + deg - max_deg == 2 * min_deg))
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (equ == max_deg))

    # 共有约束
    common_constraints(m, var, con, deg)
    common_constraints(m, var, con, deg1)
    common_constraints(m, var, con, deg2)

    m.update()
    
def pending_up(m, var, con, g, deg, index):
    # 添加辅助变量
    varPcon = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "varPcon_" + index)
    m.update()
    m.addConstr(varPcon == var + con)

    # 添加约束
    m.addConstr((varPcon == 0) >> (g == 0))

    # 共有约束
    common_constraints(m, var, con, deg)
    m.update()
    return varPcon
    

def down_pending_ups(m, deg_d, deg_p, g, deg_us, var, con, equ, index):
    # 添加辅助变量
   
    min_deg_dp = m.addVar(lb=0, vtype=GRB.INTEGER, name="min_degdp_" + index)

    m.update()


 

    m.addGenConstrMin(min_deg_dp, [deg_d, deg_p])

    # 添加约束，按 g 划分
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
    # 添加辅助变量
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
    
    # 添加约束
    m.addConstr((varPcon == 0) >> (g1 + g2 <= 1))
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (deg1 == deg2) )
    m.addConstr((BIN_g1Pg2eq1 == 1) >> (equ == 0))

    m.addConstr((BIN_g1Pg2eq0 == 1) >> (equ == max_deg))

    # 共有约束
    common_constraints(m, var, con, deg1)
    common_constraints(m, var, con, deg2)
    m.update()
  




def pending_pending_up(m, deg1, deg2, g1, g2, var, con, deg, equ, index):
    # 添加辅助变量
    min_deg = m.addVar(lb=0, vtype=GRB.INTEGER, name="min_deg_" + index)
    m.update()



    m.addGenConstrMin(min_deg, [deg1, deg2])

    # 添加约束
    pending_pending(m, deg1, deg2, g1, g2, var, con, equ, index)

    m.addConstr(deg == min_deg)

    # 共有约束
    common_constraints(m, var, con, deg)

    m.update()
   




def pending_pending_pending(m, deg1, deg2, deg3, g1, g2, g3, var, con, equ, index):
    # 添加辅助变量
  
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

    

    # 添加约束
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

    # 共有约束
    common_constraints(m, var, con, deg1)
    common_constraints(m, var, con, deg2)
    common_constraints(m, var, con, deg3)

    m.update()
    return min_deg

def pending_pending_pending_up(m, deg1, deg2, deg3, g1, g2, g3, var, con, equ, deg, index):

    min_deg = pending_pending_pending(m, deg1, deg2, deg3, g1, g2, g3, var, con, equ, index)

    #upstream的次数是三个pending中最小的那个
    m.addLConstr(deg == min_deg)

    # 共有约束
    common_constraints(m, var, con, deg)

    m.update()

