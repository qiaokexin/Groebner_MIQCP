from util4Groebner import *
from math import comb, log2

t = 8 #分支数
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
    MDS_module(m, t, [var_in[0], var_in[1], var_in[2], var_in[3],var_in[4], var_in[5], var_in[6], var_in[7], var_out[0], var_out[1], var_out[2], var_out[3],var_out[4], var_out[5], var_out[6], var_out[7],], \
               [con_in[0], con_in[1], con_in[2], con_in[3],con_in[4], con_in[5], con_in[6], con_in[7], con_out[0], con_out[1], con_out[2], con_out[3],con_out[4], con_out[5], con_out[6], con_out[7]], \
               [deg_in[0], deg_in[1], deg_in[2], deg_in[3],deg_in[4], deg_in[5], deg_in[6], deg_in[7], deg_out[0], deg_out[1], deg_out[2], deg_out[3],deg_out[4], deg_out[5], deg_out[6], deg_out[7]], \
               [g_in[0], g_in[1], g_in[2], g_in[3],g_in[4], g_in[5], g_in[6], g_in[7], g_out[0], g_out[1], g_out[2], g_out[3],g_out[4], g_out[5], g_out[6], g_out[7]], \
               [h_in[0], h_in[1], h_in[2], h_in[3],h_in[4], h_in[5], h_in[6], h_in[7], h_out[0], h_out[1], h_out[2], h_out[3],h_out[4], h_out[5], h_out[6], h_out[7]],\
               [bse_in[0], bse_in[1], bse_in[2], bse_in[3],bse_in[4], bse_in[5], bse_in[6], bse_in[7], bse_out[0], bse_out[1], bse_out[2], bse_out[3],bse_out[4], bse_out[5], bse_out[6], bse_out[7]],\
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
    # 第0,1,2,3,4,5,6,7点位是pending四元组 + 基础二元组，直接从本轮输入状态中来
    for i in range(t):
        var[i] = var_in[i]
        con[i] = con_in[i]
        deg[i] = deg_in[i]
        g[i] = g_in[i]
        h[i] = h_in[i]
        bse[i] = bse_in[i]
    # 第8点位是vcde元组
    var[8], con[8] = gen_branching_vars_vc(m, r, 8)
    deg[8], g[8], h[8], bse[8] = gen_pending_vars_dghb(m, r, 8)
    # 第9点位是vcde元组
    var[9], con[9], deg[9], equ[9]= gen_branching_vars_vcde(m, r, 9)
    # 第10点位是 pending四元组
    deg[10], g[10], h[10], bse[10] = gen_pending_vars_dghb(m, r, 10)
    # 第11点位是deg单个变量
    deg[11] = gen_vars_d(m, r, 11)
    # 第12点位是三元组
    var[12], con[12], deg[12] = gen_branching_vars_vcd(m, r, 12)
    # 第13点位是vcde元组
    var[13], con[13] = gen_branching_vars_vc(m, r, 13)
    deg[13], g[13], h[13], bse[13] = gen_pending_vars_dghb(m, r, 13)
    # 第14点位是deg单个变量
    deg[14] = gen_vars_d(m, r, 14)
    # 第15点位是三元组
    var[15], con[15], deg[15] = gen_branching_vars_vcd(m, r, 15)
    # 第16点位是vcde元组
    var[16], con[16] = gen_branching_vars_vc(m, r, 16)
    deg[16], g[16], h[16], bse[16] = gen_pending_vars_dghb(m, r, 16)
     # 第17点位是deg单个变量
    deg[17] = gen_vars_d(m, r, 17)
    # 第18点位是三元组
    var[18], con[18], deg[18] = gen_branching_vars_vcd(m, r, 18)
    # 第19点位是vcde元组
    var[19], con[19] = gen_branching_vars_vc(m, r, 19)
    deg[19], g[19], h[19], bse[19] = gen_pending_vars_dghb(m, r, 19)
    #第20点位是deg单个变量
    deg[20] = gen_vars_d(m, r, 20)
    # 第21点位是三元组
    var[21], con[21], deg[21] = gen_branching_vars_vcd(m, r, 21)
    # 第22点位是vcde元组
    var[22], con[22] = gen_branching_vars_vc(m, r, 22)
    deg[22], g[22], h[22], bse[22] = gen_pending_vars_dghb(m, r, 22)
    # 第23点位是deg单个变量
    deg[23] = gen_vars_d(m, r, 23)
    # 第24点位是三元组
    var[24], con[24], deg[24] = gen_branching_vars_vcd(m, r, 24)
    # 第25点位是vcde元组
    var[25], con[25] = gen_branching_vars_vc(m, r, 25)
    deg[25], g[25], h[25], bse[25] = gen_pending_vars_dghb(m, r, 25)
    # 第26点位是deg单个变量
    deg[26] = gen_vars_d(m, r, 26)
    # 第27点位是三元组
    var[27], con[27], deg[27] = gen_branching_vars_vcd(m, r, 27)
    # 第28点位是vcde元组
    var[28], con[28] = gen_branching_vars_vc(m, r, 28)
    deg[28], g[28], h[28], bse[28] = gen_pending_vars_dghb(m, r, 28)
    # 添加输出变量：下一轮第29,30,31,32,33,34,35,36 点位是pending四元组 + 基础二元组
    for i in [29,30,31,32,33,34,35,36]:
        var[i], con[i] = gen_branching_vars_vc(m, r+1, i-29)
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r+1, i-29)


    # 为每个运算模块添加equ变量
    # 各运算模块
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

    # 对各分支模块添加约束
    # 0 pending-up; 1 pending-up(s),同pending-up; 2 pending-up;
    # 5 pending-up(s),同pending-up
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

    # 对各运算模块添加约束
    #0,7 8 multiplication
    multiplication_module(m, [deg[0], deg[27]], var[28], con[28], deg[28], equ02728, "{}r_02728M".format(r))
    # 1,10,11 multiplication
    multiplication_module(m, [deg[1], deg[24]], var[25], con[25], deg[25], equ12425, "{}r_12425M".format(r))
    multiplication_module(m, [deg[2], deg[21]], var[22], con[22], deg[22], equ22122, "{}r_22122M".format(r))
    multiplication_module(m, [deg[3], deg[18]], var[19], con[19], deg[19], equ31819, "{}r_31819M".format(r))
    multiplication_module(m, [deg[4], deg[15]], var[16], con[16], deg[16], equ41516, "{}r_41516M".format(r))
    multiplication_module(m, [deg[5], deg[12]], var[13], con[13], deg[13], equ51213, "{}r_51213M".format(r))
    # 2,4 单变量非线性模块
    NLUP_module(m, [var[6], var[9]],[con[6], con[9]], [deg[6], deg[9]], equ69, "{}r_S1".format(r), d)
    # 5,3 单变量非线性模块
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

    # 将本轮的输出变量返回
    var_out = [var[29],var[30],var[31],var[32], var[33], var[34], var[35], var[36]]
    con_out = [con[29], con[30], con[31], con[32], con[33], con[34], con[35], con[36]]
    deg_out = [deg[29],deg[30],deg[31], deg[32],deg[33],deg[34],deg[35],deg[36]]
    g_out = [g[29],g[30],g[31], g[32],g[33],g[34],g[35],g[36]]
    h_out = [h[29],h[30],h[31], h[32],h[33],h[34],h[35],h[36]]
    bse_out = [bse[29],bse[30],bse[31], bse[32],bse[33],bse[34],bse[35],bse[36]]

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

    # 输出边间的约束
    # （1/3）CICO中的输出约束
    m.addConstr(con_in[t - 1] == 1)
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
        repeated_con_vars.append(m.getVarByName("con_{}r_9".format(r)))
        repeated_con_vars.append(m.getVarByName("con_{}r_8".format(r)))
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
    gen_Griffin4_model(4, 8, 8)