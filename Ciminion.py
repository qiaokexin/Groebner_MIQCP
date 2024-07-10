from util4Groebner import *
from math import comb, log2

def Cinimion_first_round_model(m, var_in, con_in ,deg_in, g_in,h_in, bse_in):
    # 初始化辅助变量字典
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # 变量设置
    # 第0,1,2点位是pending四元组 + 基础二元组，直接从本轮输入状态中来
    for i in range(3):
        var[i] = var_in[i]
        con[i] = con_in[i]
        deg[i] = deg_in[i]
        g[i] = g_in[i]
        h[i] = h_in[i]
        bse[i] = bse_in[i]
    # 第3点位是基础二元组 + pending四元组
    var[3], con[3] = gen_branching_vars_vc(m, 0, 3)
    deg[3], g[3], h[3], bse[3] = gen_pending_vars_dghb(m, 0, 3)
    # 第4,5,6,7,8点位是pending四元组
    for i in [4,5,6,7,8]:
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, 0, i)
    # 第9点位是vce三元组
    var[9], con[9], equ[9] = gen_branching_vars_vce(m, 0, 9)
    # 第10点位是vcde四元组，按照下一轮规则命名
    var[10], con[10], deg[10], equ[10] = gen_branching_vars_vcde(m, 1, 0)
    # 第11点位是vcde四元组，按照下一轮规则命名
    var[11], con[11], deg[11], equ[11] = gen_branching_vars_vcde(m, 1, 1)
    # 第12点位是pending四元组，按照下一轮规则命名
    deg[12], g[12], h[12], bse[12] = gen_pending_vars_dghb(m, 1, 2)
    # 第13 点位是pending四元组，按照下一轮规则命名. 这个是后来添上的
    deg[13], g[13], h[13], bse[13] = gen_pending_vars_dghb(m, 1, 14)

    # 为每个运算模块添加equ变量
    # 各运算模块
    equMUL = m.addVar(lb=0, vtype = GRB.INTEGER, name="equMUL_{}r".format(0))
    equ234 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ234_{}r".format(0))
    equ156 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ156_{}r".format(0))
    equ078 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ078_{}r".format(0))
    m.update()

    # 对各运算模块添加约束

    #0-1-3乘法
    multiplication_module(m, [deg[0], deg[1]], var[3], con[3], deg[3], equMUL, "0r_MUL")
    # 2-3-4加法
    taddition_module(m, 3, [var[2], var[3], var[10]], \
                     [con[2], con[3], con[10]],\
                     [deg[2], deg[3], deg[4]], [g[2], g[3], g[4]],\
                     [h[2], h[3], h[4]], [bse[2], bse[3], bse[4]],\
                     equ234, "0r_234")
    # 1-5-6加法
    taddition_module(m, 3, [var[1], var[10], var[9]],\
                     [con[1], con[10], con[9]],\
                     [deg[1], deg[5], deg[6]],\
                     [g[1], g[5], g[6]],\
                     [h[1], h[5], h[6]],\
                     [bse[1], bse[5], bse[6]], equ156, "0r_156")
    # 0-7-8 加法
    taddition_module(m, 3, [var[0], var[9], var[11]], \
                     [con[0], con[9], con[11]], \
                     [deg[0], deg[7], deg[8]], \
                     [g[0], g[7], g[8]], \
                     [h[0], h[7], h[8]], \
                     [bse[0], bse[7], bse[8]], equ078, "0r_078")

    # 对各分支模块添加约束
    # 3点位 down-pending
    down_pending(m, var[3], con[3], g[3], deg[3], "0r_3")
    # 4-5-10-13 pending-pending-pending-up
    pending_pending_pending_up(m, deg[4], deg[5], deg[13], g[4], g[5], g[13], var[10], con[10], equ[10], deg[10], "0r_451013")
    # 6-7-12-9 pending-pending-pending
    pending_pending_pending(m, deg[6], deg[7], deg[12], g[6], g[7], g[12], var[9], con[9], equ[9], "0r_67129")

    # 将本轮的输出变量返回
    # 第0分支vcde变量
    out0 = [var[10], con[10], deg[10], equ[10]]
    # 第1分支vcde变量
    out1 = [var[11], con[11], deg[11], equ[11]]
    # 第2分支pending四元组
    out2 = [deg[12], g[12], h[12], bse[12]]
    # 8点位，9点位，14点位在下一轮也需要使用，也需要返回
    out8 = [deg[8], g[8], h[8], bse[8]]
    out9 = [var[9], con[9], equ[9]]
    out13 = [deg[13], g[13], h[13], bse[13]]

    return out0, out1, out2, out8, out9, out13

def Cinimion_middle_round_model(m, r, in0, in1, in2, in8, in9, in13):
    # 初始化辅助变量字典
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # 变量设置
    # 第0,1点位是vcde四元组，直接从本轮输入状态中来
    var[0], con[0], deg[0], equ[0] = in0
    var[1], con[1], deg[1], equ[1] = in1
    # 第2点位是pending四元组，直接从本轮输入状态中来
    deg[2], g[2], h[2], bse[2] = in2
    # 第14点位是pending四元组，直接从输入中来
    deg[14], g[14], h[14], bse[14] = in13
    # 第3点位是基础二元组 + pending四元组
    var[3], con[3] = gen_branching_vars_vc(m, r, 3)
    deg[3], g[3], h[3], bse[3] = gen_pending_vars_dghb(m, r, 3)
    # 第4,5,6,7,8,10,11点位是pending四元组
    for i in [4,5,6,7,8,15]:
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r, i)
    # 第9点位是vce三元组
    var[9], con[9], equ[9] = gen_branching_vars_vce(m, r, 9)
    # 第10点位是vcde四元组，按照下一轮规则命名
    var[10], con[10], deg[10], equ[10] = gen_branching_vars_vcde(m, r+1, 0)
    # 第11点位是vcde四元组，按照下一轮规则命名
    var[11], con[11], deg[11], equ[11] = gen_branching_vars_vcde(m, r+1, 1)
    # 第12点位是pending四元组，按照下一轮规则命名
    deg[12], g[12], h[12], bse[12] = gen_pending_vars_dghb(m, r+1, 2)
    # 第13 点位是pending四元组，按照下一轮规则命名. 这个是后来添上的
    deg[13], g[13], h[13], bse[13] = gen_pending_vars_dghb(m, r+1, 14)

    # 为每个运算模块添加equ变量
    # 各运算模块
    equMUL = m.addVar(lb=0, vtype = GRB.INTEGER, name="equMUL_{}r".format(r))
    equ234 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ234_{}r".format(r))
    equ1356 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ1356_{}r".format(r))
    equ1478 = m.addVar(lb=0, vtype = GRB.INTEGER, name="equ1478_{}r".format(r))
    m.update()

    # 对各运算模块添加约束


    #0-1-3乘法
    multiplication_module(m, [deg[0], deg[1]], var[3], con[3], deg[3], equMUL, "{}r_MUL".format(r))
    # 2-3-4加法
    taddition_module(m, 3, [in9[0], var[3], var[10]], \
                     [in9[1], con[3], con[10]],\
                     [deg[2], deg[3], deg[4]], [g[2], g[3], g[4]],\
                     [h[2], h[3], h[4]], [bse[2], bse[3], bse[4]],\
                     equ234, "{}r_234".format(r))
    # 15-5-6加法
    taddition_module(m, 3, [var[1], var[10], var[9]],\
                     [con[1], con[10], con[9]],\
                     [deg[15], deg[5], deg[6]],\
                     [g[15], g[5], g[6]],\
                     [h[15], h[5], h[6]],\
                     [bse[15], bse[5], bse[6]], equ1356, "{}r_1356".format(r))
    # 14-7-8 加法
    taddition_module(m, 3, [var[0], var[9], var[11]], \
                     [con[0], con[9], con[11]], \
                     [deg[14], deg[7], deg[8]], \
                     [g[14], g[7], g[8]], \
                     [h[14], h[7], h[8]], \
                     [bse[14], bse[7], bse[8]], equ1478, "{}r_1478".format(r))
    
    # 对各分支模块添加约束

    # pre8-1-15 pending-pending-up
    pending_pending_up(m, in8[0], deg[15], in8[1], g[15], var[1], con[1], deg[1], equ[1], "{}r_pre8115".format(r))

    # 3点位 down-pending
    down_pending(m, var[3], con[3], g[3], deg[3], "{}r_3".format(r))
    # 4-5-10-13 pending-pending-pending-up
    pending_pending_pending_up(m, deg[4], deg[5], deg[13], g[4], g[5], g[13], var[10], con[10], equ[10], deg[10], "{}r_451013".format(r))
    # 6-7-12-9 pending-pending-pending
    pending_pending_pending(m, deg[6], deg[7], deg[12], g[6], g[7], g[12], var[9], con[9], equ[9], "{}r_67129".format(r))

    # 将本轮的输出变量返回
    # 第0分支vcde变量
    out0 = [var[10], con[10], deg[10], equ[10]]
    # 第1分支vcde变量
    out1 = [var[11], con[11], deg[11], equ[11]]
    # 第2分支pending四元组
    out2 = [deg[12], g[12], h[12], bse[12]]
    # 8点位，9点位，4点位在下一轮也需要使用，也需要返回
    out8 = [deg[8], g[8], h[8], bse[8]]
    out9 = [var[9], con[9], equ[9]]
    out13 = [deg[13], g[13], h[13], bse[13]]

    return out0, out1, out2, out8, out9, out13

def Cinimion_last_round_model(m, r, in0, in1, in2, in8, in9, in13):
    # 初始化辅助变量字典
    var = {}
    con = {}
    deg = {}
    g = {}
    h = {}
    bse = {}
    equ = {}

    # 变量设置
    # 第0,1点位是vcde四元组，直接从本轮输入状态中来
    var[0], con[0], deg[0], equ[0] = in0
    var[1], con[1], deg[1], equ[1] = in1
    # 第2点位是pending四元组，直接从本轮输入状态中来
    deg[2], g[2], h[2], bse[2] = in2
    # 第14点位是pending四元组，直接从输入中来
    deg[14], g[14], h[14], bse[14] = in13
    # 第3点位是基础二元组 + pending四元组
    var[3], con[3] = gen_branching_vars_vc(m, r, 3)
    deg[3], g[3], h[3], bse[3] = gen_pending_vars_dghb(m, r, 3)
    # 第4,5,6,7,8,10,11点位是pending四元组
    for i in [4, 5, 6, 7, 8, 15]:
        deg[i], g[i], h[i], bse[i] = gen_pending_vars_dghb(m, r, i)
    # 第8点位增加基础二元组
    var[8], con[8] = gen_branching_vars_vc(m, r, 8)
    # 第9点位是vce三元组
    var[9], con[9], equ[9] = gen_branching_vars_vce(m, r, 9)
    # 第10点位是vce三元组
    var[10], con[10], equ[10] = gen_branching_vars_vce(m, r, 10)
    # 第11,12,13点位不需要

    # 为每个运算模块添加equ变量
    # 各运算模块
    equMUL = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMUL_{}r".format(r))
    equ234 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ234_{}r".format(r))
    equ1356 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1356_{}r".format(r))
    equ1478 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1478_{}r".format(r))
    m.update()

    # 对各运算模块添加约束


    # 0-1-3乘法
    multiplication_module(m, [deg[0], deg[1]], var[3], con[3], deg[3], equMUL, "{}r_MUL".format(r))
    # 2-3-4加法
    taddition_module(m, 3, [in9[0], var[3], var[10]], \
                     [in9[1], con[3], con[10]], \
                     [deg[2], deg[3], deg[4]], [g[2], g[3], g[4]], \
                     [h[2], h[3], h[4]], [bse[2], bse[3], bse[4]], \
                     equ234, "{}r_234".format(r))
    # 15-5-6加法
    taddition_module(m, 3, [var[1], var[10], var[9]], \
                     [con[1], con[10], con[9]], \
                     [deg[15], deg[5], deg[6]], \
                     [g[15], g[5], g[6]], \
                     [h[15], h[5], h[6]], \
                     [bse[15], bse[5], bse[6]], equ1356, "{}r_1556".format(r))
    # 14-7-8 加法
    taddition_module(m, 3, [var[0], var[9], var[8]], \
                     [con[0], con[9], con[8]], \
                     [deg[14], deg[7], deg[8]], \
                     [g[14], g[7], g[8]], \
                     [h[14], h[7], h[8]], \
                     [bse[14], bse[7], bse[8]], equ1478, "{}r_1478".format(r))

    # 对各分支模块添加约束
    # pre8-1-15 pending-pending-up
    pending_pending_up(m, in8[0], deg[15], in8[1], g[15], var[1], con[1], deg[1], equ[1], "{}r_pre8115".format(r))
    # 3点位 down-pending
    down_pending(m, var[3], con[3], g[3], deg[3], "{}r_3".format(r))
    # 4-5-10 pending-pending
    pending_pending(m, deg[4], deg[5], g[4], g[5], var[10], con[10], equ[10], "{}r_4510".format(r))
    # 6-7-9 pending-pending
    pending_pending(m, deg[6], deg[7], g[6], g[7], var[9], con[9], equ[9], "{}r_679".format(r))

    # 将本轮的输出变量返回
    var_out = {}
    con_out = {}
    deg_out = {}
    g_out = {}
    # 第0分支vc变量,及与之关联的deg，g
    var_out[0] = var[10]
    con_out[0] = con[10]
    deg_out[0] = [deg[4], deg[5]]
    g_out[0] = [g[4], g[5]]
    # 第1分支vc变量, 及与之关联的deg，g
    var_out[1] = var[8]
    con_out[1] = con[8]
    deg_out[1] = [deg[8]]
    g_out[1] = [g[8]]
    # 第2分支vc变量, 及与之关联的deg，g
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
    # 创建模型
    m = gp.Model("Cinimion_CICO_MILP_{}r_{}v".format(R, numVar))
    # 创建输入变量：pending四元组+输入状态二元组
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

    # 输入边界的约束
    # （1/3） CICO中的输入约束
    m.addConstr(con_in[0] == 1)

    # （2/3）输入边界添加 var=0,con=0 则 g=0的约束
    varPcon = {}
    for i in range(3):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_MDS_{}".format(-1, i))
    m.update()
    for i in range(3):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） 共有约束
    for i in range(3):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # 首轮约束
    in0, in1, in2, in8, in9, in13 = Cinimion_first_round_model(m, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # 中间轮约束
    for r in range(1,R-1):
        in0, in1, in2, in8, in9, in13 = Cinimion_middle_round_model(m, r, in0, in1, in2, in8, in9, in13)

    # 最后一轮约束
    var_out, con_out, deg_out, g_out = Cinimion_last_round_model(m, R-1, in0, in1, in2, in8, in9, in13)

    # 输出边界的约束
    # （1/3）CICO中的输出约束
    m.addConstr(con_out[0] == 1)
    #m.addConstr(con_out[1] == 1)
    # （2/3）输出边界添加 var=0,con=0 则 g=0的约束
    for i in range(3):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_" + "{}r_{}".format(R, i))
    m.update()
    for i in range(3):
        m.addConstr(varPcon[i] == var_out[i] + con_out[i])
        for cg in g_out[i]:
            m.addConstr((varPcon[i] == 0) >> (cg == 0))
    m.update()
    # （3/3） 共有约束
    for i in range(3):
        for cdeg in deg_out[i]:
            common_constraints(m, var_out[i], con_out[i], cdeg)

    # 常数个数
    all_con_vars = [var for var in m.getVars() if "con_" in var.VarName[:4]]
    m.addConstr(sum(all_con_vars) == 3)  # 常数个数等于t，则平均有1个解

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

    #m.addConstr(m.getVarByName("con_0r_0")==1)
    #m.addConstr(m.getVarByName("con_2r_10") == 1)
    #m.addConstr(m.getVarByName("con_2r_8") == 1)

    #m.addConstr(m.getVarByName("var_0r_1")==1)
    #m.addConstr(m.getVarByName("var_0r_2") == 1)



    # 模型设置
    m.setParam("TimeLimit", 600)  # 设定跑600秒后就停止，输出当前找到的最优解

    obj_bound = get_obj_bound(2*R+1, numVar, R, 2)
    #m.setParam("BestObjStop", obj_bound) #设定目标函数低于100时就停止
    print("It's better to reach obj_bound", obj_bound)
    m.write("Cinimion_CICO_MILP_{}r_{}v.lp".format(R, numVar))
    m.optimize()
    m.write("Cinimion_MILP_{}r_{}v.sol".format(R, numVar))

    parsing_Cinimion_solution(m, R)
    print('轮数：{}，'.format(R))
    print('变量个数：{}，'.format(numVar))
    print('方程次数之和：{}'.format(obj.X))
    print('求解复杂度：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))
    print('用于单变量下：2^{}'.format(univariate_comp(int(obj.X))))

if __name__ == '__main__':
    gen_Cinimion_model(10, 2)
