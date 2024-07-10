from util4Groebner import *
from math import comb, log2

t = 4 #分支数
d = 5 #S盒次数

def Griffin4_initM_model(m, t, var_in, con_in ,deg_in, g_in,h_in, bse_in):
    var_out = {}
    con_out = {}
    deg_out = {}
    g_out = {}
    h_out = {}
    bse_out = {}
    # 添加输出变量
    for i in range(t):
        var_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_0r_{}".format(i))
        con_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_0r_{}".format(i))
        deg_out[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_0r_{}".format(i))
        g_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_0r_{}".format(i))
        h_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_0r_{}".format(i))
        bse_out[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_0r_{}".format(i))

    equ_M = m.addVar(lb=0, vtype = GRB.INTEGER, name="equMDS_-1r")
    m.update()

    # MDS 约束
    MDS_module(m, t, [var_in[0], var_in[1], var_in[2], var_in[3], var_out[0], var_out[1], var_out[2], var_out[3]], \
               [con_in[0], con_in[1], con_in[2], con_in[3], con_out[0], con_out[1], con_out[2], con_out[3]], \
               [deg_in[0], deg_in[1], deg_in[2], deg_in[3], deg_out[0], deg_out[1], deg_out[2], deg_out[3]], \
               [g_in[0], g_in[1], g_in[2], g_in[3], g_out[0], g_out[1], g_out[2], g_out[3]], \
               [h_in[0], h_in[1], h_in[2], h_in[3], h_out[0], h_out[1], h_out[2], h_out[3]],\
               [bse_in[0], bse_in[1], bse_in[2], bse_in[3], bse_out[0], bse_out[1], bse_out[2], bse_out[3]],\
               equ_M,  "{}r_MDS".format(-1))


    return var_out, con_out, deg_out, g_out, h_out, bse_out

def Griffin4_round_model(R, m, r, t, var_in, con_in ,deg_in, g_in,h_in, bse_in):
    # 对第r轮轮函数建模，t是分支数

    # 初始化辅助变量字典
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # 变量设置
    # 第0,1,2,3点位是pending四元组 + 基础二元组，直接从本轮输入状态中来
    for i in range(t):
        var[i] = var_in[i]
        con[i] = con_in[i]
        deg[i] = deg_in[i]
        g[i] = g_in[i]
        h[i] = h_in[i]
        bse[i] = bse_in[i]
    # 第4点位是vcde元组
    var[4], con[4], deg[4], equ[4]= gen_branching_vars_vcde(m, r, 4)
    # 第5点位是基础二元组 + pending四元组
    var[5], con[5] = gen_branching_vars_vc(m, r, 5)
    deg[5], g[5], h[5], bse[5] = gen_pending_vars_dghb(m, r, 5)
    # 第6点位是deg单个变量
    deg[6] = gen_vars_d(m, r, 6)
    # 第7点位是三元组
    var[7], con[7], deg[7] = gen_branching_vars_vcd(m, r, 7)
    # 第8点位是pending四元组 + 基础二元组
    var[8], con[8] = gen_branching_vars_vc(m, r, 8)
    deg[8], g[8], h[8], bse[8] = gen_pending_vars_dghb(m, r, 8)
    # 第9点位是deg单个变量
    deg[9] = gen_vars_d(m, r, 9)
    # 第10点位是三元组
    var[10], con[10], deg[10] = gen_branching_vars_vcd(m, r, 10)
    # 第11点位是pending四元组 + 基础二元组
    var[11], con[11] = gen_branching_vars_vc(m, r, 11)
    deg[11], g[11], h[11],bse[11] = gen_pending_vars_dghb(m, r, 11)
    # 第12点位是pending四元组
    deg[12], g[12], h[12], bse[12] = gen_pending_vars_dghb(m, r, 12)
    # 添加输出变量：下一轮第13， 14， 15， 16 点位是pending四元组 + 基础二元组
    for i in [13,14,15,16]:
        var[i], con[i] = gen_branching_vars_vc(m, r+1, i-13)
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r+1, i-13)


    # 为每个运算模块添加equ变量
    # 各运算模块
    equ078 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ078_{}r".format(r))
    equG2 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equG2_{}r".format(r))
    equG1 = m.addVar(lb = 0, vtype=GRB.INTEGER, name="equG1_{}r".format(r))
    equ11011 = m.addVar(lb = 0, vtype=GRB.INTEGER, name="equ11011_{}r".format(r))
    equ24 = m.addVar(lb = 0, vtype=GRB.INTEGER, name="equ24_{}r".format(r))
    equ53 = m.addVar(lb = 0, vtype=GRB.INTEGER, name="equ35_{}r".format(r))
    equMDS = m.addVar(lb=0, vtype=GRB.INTEGER, name = "equMDS_{}r".format(r))
    m.update()

    # 对各分支模块添加约束
    # 0 pending-up; 1 pending-up(s),同pending-up; 2 pending-up;
    # 5 pending-up(s),同pending-up
    for i in [0, 1, 2, 5]:
        pending_up(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r,i))
    # 3 down-pending; 8 down-pending; 11 down-pending
    for i in [3, 8, 11]:
        down_pending(m, var[i], con[i], g[i], deg[i], "{}r_{}".format(r, i))
    # 7, 10 down-up
    for i in [7, 10]:
        down_up(m, var[i], con[i], deg[i])
    # 4-6-9-12 down-pending-ups
    down_pending_ups(m, deg[4],deg[12], g[12], [deg[6], deg[9]], var[4], con[4],equ[4], "{}r_4".format(r))

    # 对各运算模块添加约束
    #0,7 8 multiplication
    multiplication_module(m, [deg[0], deg[7]], var[8], con[8], deg[8], equ078, "{}r_078M".format(r))
    # 1,10,11 multiplication
    multiplication_module(m, [deg[1], deg[10]], var[11], con[11], deg[11], equ11011, "{}r_11011M".format(r))
    # 2,4 单变量非线性模块
    NLUP_module(m, [var[2], var[4]],[con[2], con[4]], [deg[2], deg[4]], equ24, "{}r_S1".format(r), d)
    # 5,3 单变量非线性模块
    NLUP_module(m, [var[5], var[3]], [con[5], con[3]], [deg[5], deg[3]], equ53, "{}r_S2".format(r), d)
    # G2
    NLMP_module(m, [deg[5], deg[6], deg[1]], var[7], con[7], deg[7], equG2, "{}r_G2".format(r), 2)
    # G1
    NLMP_module(m, [deg[5], deg[9]], var[10], con[10], deg[10], equG1, "{}r_G1".format(r), 2)
    # MDS

    MDS_module(m, 4, [var[8],var[11],var[4],var[5], var[13], var[14], var[15], var[16]], \
               [con[8], con[11], con[4], con[5], con[13], con[14], con[15], con[16]],\
               [deg[8],deg[11],deg[12], deg[5],deg[13],deg[14],deg[15],deg[16]],\
               [g[8],g[11],g[12], g[5],g[13],g[14],g[15],g[16]],\
               [h[8],h[11],h[12], h[5],h[13],h[14],h[15],h[16]],\
               [bse[8],bse[11],bse[12], bse[5],bse[13],bse[14],bse[15],bse[16]],\
               equMDS,  "{}r_MDS".format(r))

    # 将本轮的输出变量返回
    var_out = [var[13], var[14], var[15], var[16]]
    con_out = [con[13], con[14], con[15], con[16]]
    deg_out = [deg[13],deg[14],deg[15],deg[16]]
    g_out = [g[13],g[14],g[15],g[16]]
    h_out = [h[13],h[14],h[15],h[16]]
    bse_out = [bse[13],bse[14],bse[15],bse[16]]

    m.update()
    return var_out, con_out, deg_out, g_out, h_out, bse_out

def parsing_Griffin4_solution(m, R):
    # 初始MDS
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
    # 每轮
    for r in range(R):
        print('--------------r = {}-----------------'.format(r))
        print('var_{}r'.format(r), end=" ")
        for j in [8,11,4,5]:
            print(round(m.getVarByName("var_{}r_{}".format(r,j)).X), end=" ")
        print('\nvar_{}r'.format(r+1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("var_{}r_{}".format(r+1, j)).X), end=" ")
        print('\ncon_{}r'.format(r), end=" ")
        for j in [8, 11, 4, 5]:
            print(round(m.getVarByName("con_{}r_{}".format(r, j)).X), end=" ")
        print('\ncon_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("con_{}r_{}".format(r + 1, j)).X), end=" ")
        print('\ndeg_{}r'.format(r), end=" ")
        for j in [8, 11, 12, 5]:
            print(round(m.getVarByName("deg_{}r_{}".format(r, j)).X), end=" ")
        print('\ndeg_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("deg_{}r_{}".format(r + 1, j)).X), end=" ")
        print('\nbse_{}r'.format(r), end=" ")
        for j in [8, 11, 12, 5]:
            print(round(m.getVarByName("bse_{}r_{}".format(r, j)).X), end=" ")
        print('\nbse_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("bse_{}r_{}".format(r + 1, j)).X), end=" ")
        print('\ng_{}r'.format(r), end=" ")
        for j in [8, 11, 12, 5]:
            print(round(m.getVarByName("g_{}r_{}".format(r, j)).X), end=" ")
        print('\ng_{}r'.format(r + 1), end=" ")
        for j in range(t):
            print(round(m.getVarByName("g_{}r_{}".format(r + 1, j)).X), end=" ")
        print("\nequMDS_{}r".format(r), round(m.getVarByName("equMDS_{}r".format(r)).X))

        for j in [7, 10]:
            print("var_{}r_{} = ".format(r, j), end = " ")
            print(round(m.getVarByName("var_{}r_{}".format(r,j)).X))
        for j in [7, 10]:
            print("con_{}r_{} = ".format(r, j), end = " ")
            print(round(m.getVarByName("con_{}r_{}".format(r,j)).X))
    for var in m.getVars():
        if "equ" in var.VarName[:4]:
            if var.X > 0:
                print(var.VarName, "=", var.X)
        if "max_bseTdeg" in var.VarName:
            print(var.VarName, "=", var.X)
        if "sum_con" in var.VarName:
            print(var.VarName, "=", var.X)
def gen_Griffin4_model(R, t, numVar):
    # R: 轮数; t: 分支数； numVar: 引入变量个数

    # 创建模型
    m = gp.Model("Griffin4_CICO_MILP_{}r_{}t_{}v".format(R, t, numVar))
    # 创建输入变量：pending四元组+输入状态二元组
    var_in={};con_in={};deg_in={};g_in={};h_in={};bse_in={};
    for i in range(t):
        var_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_-1r_{}".format(i))
        con_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="con_-1r_{}".format(i))
        deg_in[i] = m.addVar(lb=0, vtype=GRB.INTEGER, name="deg_-1r_{}".format(i))
        g_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="g_-1r_{}".format(i))
        h_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="h_-1r_{}".format(i))
        bse_in[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="bse_-1r_{}".format(i))
    m.update()

    # 输入边界的约束
    # （1/3） CICO中的输入约束
    m.addConstr(con_in[0] == 1)


    # （2/3）输入边界添加 var=0,con=0 则 g=0的约束
    varPcon = {}
    for i in range(t):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_MDS_{}".format(-1, i))
    m.update()
    for i in range(t):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） 共有约束
    for i in range(t):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # 初始MDS变换
    var_in, con_in, deg_in, g_in, h_in, bse_in = Griffin4_initM_model(m, t, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # 每轮建立模型
    for r in range(R):
        var_in, con_in, deg_in, g_in, h_in, bse_in = \
            Griffin4_round_model(R, m, r, t, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # 输出边界的约束
    # （1/3）CICO中的输出约束
    m.addConstr(con_in[t-1] == 1)
    # （2/3）输出边界添加 var=0,con=0 则 g=0的约束
    for i in range(t):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_{}".format(R, i))
    m.update()
    for i in range(t):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） 共有约束
    for i in range(t):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])


    # 常数个数
    all_con_vars = [var for var in m.getVars() if "con_" in var.VarName[:4]]
    repeated_con_vars = [] # S盒前后的常数只取一个，把重复的去掉
    for r in range(R):
        repeated_con_vars.append(m.getVarByName("con_{}r_5".format(r)))
        repeated_con_vars.append(m.getVarByName("con_{}r_4".format(r)))
    sum_cons = sum(all_con_vars) - sum(repeated_con_vars)
    m.addConstr(sum_cons == t) #常数个数等于t，则平均有1个解

    # 给定变量个数
    all_vars = [var for var in m.getVars() if "var_" in var.VarName]
    m.addConstr(sum(all_vars) == numVar)

    # 目标函数
    # 方程次数之和
    all_equs = [var for var in m.getVars() if "equ" in var.VarName]
    obj = m.addVar(vtype=GRB.INTEGER, name="obj")
    m.update()
    m.addConstr(obj == sum(all_equs))
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()

    # 人为添加约束
    '''
    m.addConstr(m.getVarByName("con_0r_3")==1)
    m.addConstr(m.getVarByName("con_0r_2") == 1)
    m.addConstr(m.getVarByName("con_0r_5") == 1)

    m.addConstr(m.getVarByName("var_-1r_3")==1)
    m.addConstr(m.getVarByName("var_1r_5") == 1)
    '''



    # 模型设置
    m.setParam("TimeLimit", 600) # 设定跑600秒后就停止，输出当前找到的最优解
    # m.setParam("BestObjStop", 100) #设定目标函数低于100时就停止
    m.write("Griffin4_CICO_MILP_{}rrrrrr_{}t_{}v.lp".format(R, t, numVar))
    m.optimize()
    m.write("Griffin4_MILP_{}r_{}t_{}v.sol".format(R, t, numVar))
    parsing_Griffin4_solution(m, R)
    print("分支数：{}，".format(t))
    print('轮数：{}，'.format(R))
    print('变量个数：{}，'.format(numVar))
    print('方程次数之和：{}'.format(obj.X))
    print('求解复杂度：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))

if __name__ == '__main__':
    gen_Griffin4_model(6, 4, 12)
