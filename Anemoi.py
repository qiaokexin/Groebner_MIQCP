from util4Groebner import *
from math import comb, log2

alpha = 3
def Anemoi_round_model(m, r, l, var_in, con_in, deg_in, g_in, h_in, bse_in):
    # 对第 r 轮轮函数建模，l是MDS的宽度, var_in及后面的变量是本轮的输入变量

    # 初始化辅助变量字典
    var={}; con ={}; deg = {}; g = {}; h = {}; bse = {}; equ = {};
    # 变量设置
    for i in range(l):
        # 对各个分支模块添加变量
        # 第 0个点位是pending四元组+基础二元组，直接从输入中获得
        var[26 * i + 0] = var_in[i]
        con[26 * i + 0] = con_in[i]
        deg[26 * i + 0] = deg_in[i]
        g[26 * i + 0] = g_in[i]
        h[26 * i + 0] = h_in[i]
        bse[26 * i + 0] = bse_in[i]
        # 第 1 个点位是pending四元组+基础二元组，直接从输入中获得
        var[26 * i + 1] = var_in[l + i]
        con[26 * i + 1] = con_in[l + i]
        deg[26 * i + 1] = deg_in[l + i]
        g[26 * i + 1] = g_in[l + i]
        h[26 * i + 1] = h_in[l + i]
        bse[26 * i + 1] = bse_in[l + i]
        # 第 2 个点位是pending四元组
        deg[26 * i + 2], g[26 * i + 2], h[26 * i + 2], bse[26 * i + 2] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,2))
        # 第 3 个点位是pending四元组
        deg[26 * i + 3], g[26 * i + 3], h[26 * i + 3], bse[26 * i + 3] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,3))
        # 第 4，5 个点位是var,con,equ三元组
        var[26 * i + 4], con[26 * i + 4], equ[26 * i + 4] = gen_branching_vars_vce(m, r, "{}th_{}".format(i,4))
        var[26 * i + 5], con[26 * i + 5], equ[26 * i + 5] = gen_branching_vars_vce(m, r, "{}th_{}".format(i,5))
        # 第6,7,8,9点位是pending四元组
        for j in range(6, 10):
            deg[26 * i + j], g[26 * i + j], h[26 * i + j], bse[26 * i + j] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,j))
        # 第 10 个点位是var,con,equ三元组
        var[26 * i + 10], con[26 * i + 10], equ[26 * i + 10] = gen_branching_vars_vce(m, r, "{}th_{}".format(i,10))
        # 第 11 个点位是pending四元组
        deg[26 * i + 11], g[26 * i + 11], h[26 * i + 11], bse[26 * i + 11] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,11))
        # 第 12 个点位是var,con,deg,equ四元组
        var[26 * i + 12], con[26 * i + 12], deg[26 * i + 12], equ[26 * i + 12] = gen_branching_vars_vcde(m, r, "{}th_{}".format(i,12))
        # 第 13 个点位是pending四元组+二元组
        deg[26 * i + 13], g[26 * i + 13], h[26 * i + 13], bse[26 * i + 13] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,13))
        var[26 * i + 13], con[26 * i + 13] = gen_branching_vars_vc(m, r, "{}th_{}".format(i,13))
        # 第 14, 15个点位是pending四元组
        deg[26 * i + 14], g[26 * i + 14], h[26 * i + 14], bse[26 * i + 14] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,14))
        deg[26 * i + 15], g[26 * i + 15], h[26 * i + 15], bse[26 * i + 15] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,15))
        # 第 16 个点位是pending四元组+二元组
        deg[26 * i + 16], g[26 * i + 16], h[26 * i + 16], bse[26 * i + 16] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,16))
        var[26 * i + 16], con[26 * i + 16] = gen_branching_vars_vc(m, r, "{}th_{}".format(i,16))
        # 第 17 个点位是var,con,deg,equ四元组
        var[26 * i + 17], con[26 * i + 17], deg[26 * i + 17], equ[26 * i + 17] = gen_branching_vars_vcde(m, r, "{}th_{}".format(i,17))
        # 第 18, 19个点位是pending四元组
        deg[26 * i + 18], g[26 * i + 18], h[26 * i + 18], bse[26 * i + 18] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,18))
        deg[26 * i + 19], g[26 * i + 19], h[26 * i + 19], bse[26 * i + 19] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,19))
        # 第 20 个点位是var,con,deg,equ四元组，但其中var, con按照下一轮的指标命名
        deg[26 * i + 20] = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "deg_{}r_{}th_{}".format(r, i, 20))
        equ[26 * i + 20] = m.addVar(lb = 0, vtype = GRB.INTEGER, name = "equ_{}r_{}th_{}".format(r, i, 20))
        var[26 * i + 20], con[26 * i + 20] = gen_branching_vars_vc(m, r + 1, "{}th_{}".format(i, 1))
        # 第 21 个点位是pending四元组+二元组
        deg[26 * i + 21], g[26 * i + 21], h[26 * i + 21], bse[26 * i + 21] = gen_pending_vars_dghb(m, r, "{}th_{}".format(i,21))
        var[26 * i + 21], con[26 * i + 21] = gen_branching_vars_vc(m, r, "{}th_{}".format(i,21))
        # 第 22 个点位是pending四元组
        deg[26 * i + 22], g[26 * i + 22], h[26 * i + 22], bse[26 * i + 22] = gen_pending_vars_dghb(m, r,
                                                                                                   "{}th_{}".format(i,
                                                                                                                    22))
        # 第 23 个点位是equ+vc二元组，其中vc二元组按照下一轮的指标命名
        equ[26 * i + 23] = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ_{}r_{}th_{}".format(r, i, 23))
        var[26 * i + 23], con[26 * i + 23] = gen_branching_vars_vc(m, r + 1, "{}th_{}".format(i, 0))


        # 第 24,25点位是下一轮的输入,pending四元组,按照下一轮的角标取名
        deg[26 * i + 24], g[26 * i + 24], h[26 * i + 24], bse[26 * i + 24] = gen_pending_vars_dghb(m, r+1, "{}th_{}".format(i,0))
        deg[26 * i + 25], g[26 * i + 25], h[26 * i + 25], bse[26 * i + 25] = gen_pending_vars_dghb(m, r + 1,
                                                                                                   "{}th_{}".format(i,1))


        # 为每个运算模块添加equ变量
        # 各运算模块
        equ_smallMDS= m.addVar(lb=0, vtype=GRB.INTEGER, name="equ_smallMDS_{}r_{}".format(r, i))
        equ111314 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ111314_{}r_{}".format(r, i))
        equ1213 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1213_{}r_{}".format(r, i))
        equ1617 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ1617_{}r_{}".format(r, i))
        equ151619 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ151619_{}r_{}".format(r, i))
        equ2021 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ2021_{}r_{}".format(r, i))
        equ182122 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equ182122_{}r_{}".format(r, i))

        # 对各分支模块添加约束
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


        # 对各运算模块添加约束
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

    # 若l大于等于2，添加仿射层
    if (l >= 2):
        equMDS1 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS1_{}r_{}".format(r, i))
        equMDS2 = m.addVar(lb=0, vtype=GRB.INTEGER, name="equMDS2_{}r_{}".format(r, i))

        # 准备MDS模块的输入和输出变量
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
            # 第一个MDS的变量
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
            # 第二个MDS的变量
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
            # 第一个MDS的约束
        MDS_module(m, l, MDS_var1, MDS_con1, MDS_deg1, MDS_g1, MDS_h1, MDS_bse1, equMDS1, "_{}r_MDS1".format(r))
        # 第二个MDS的约束
        MDS_module(m, l, MDS_var2, MDS_con2, MDS_deg2, MDS_g2, MDS_h2, MDS_bse2, equMDS2, "_{}r_MDS2".format(r))

    # 将本轮的输出变量返回
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
    # 添加输出变量
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

    # MDS 约束
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
    # R是轮数，l是一半的分支数, numVar是设定的变量的个数

    # 创建模型
    m = gp.Model("Anemoi_CICO_MILP_{}r_{}l_{}v".format(R, l, numVar))
    # 创建输入变量：pending四元组+输入状态二元组
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

    # 输入边界的约束
    # （1/3） CICO中的输入约束
    m.addConstr(con_in[0] == 1)

    m.update()

    # （2/3）输入边界添加 var=0,con=0 则 g=0的约束
    varPcon = {}
    for i in range(2*l):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_in_{}".format(i))
    m.update()
    for i in range(2*l):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） 共有约束
    for i in range(2*l):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # 每轮建立模型
    for r in range(R):
        var_in, con_in, deg_in, g_in, h_in, bse_in = Anemoi_round_model(m, r, l, var_in, con_in, deg_in, g_in, h_in,
                                                                        bse_in)
    # 最后的MDS
    var_in, con_in, deg_in, g_in, h_in, bse_in = last_MDS_model(m, l, var_in, con_in, deg_in, g_in, h_in, bse_in)

    # 输出边界的约束
    # （1/3）CICO中的输出约束
    m.addConstr(con_in[0] == 1)

    m.update()

    # （2/3）输出边界添加 var=0,con=0 则 g=0的约束
    varPcon = {}
    for i in range(2 * l):
        varPcon[i] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="varPcon_out_{}".format(i))
    m.update()
    for i in range(2 * l):
        m.addConstr(varPcon[i] == var_in[i] + con_in[i])
        m.addConstr((varPcon[i] == 0) >> (g_in[i] == 0))
    m.update()
    # （3/3） 共有约束
    for i in range(2 * l):
        common_constraints(m, var_in[i], con_in[i], deg_in[i])

    # 人为添加约束
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

    # 常数个数
    all_con_vars = [var for var in m.getVars() if "con_" in var.VarName[:4]]
    repeated_con_vars = []  # S盒前后的常数只取一个，把重复的去掉
    for r in range(R):
        for i in range(l):
            repeated_con_vars.append(m.getVarByName("con_{}r_{}th_{}".format(r, i, 13)))
            repeated_con_vars.append(m.getVarByName("con_{}r_{}th_{}".format(r, i, 17)))
            repeated_con_vars.append(m.getVarByName("con_{}r_{}th_{}".format(r, i, 21)))
    sum_cons = sum(all_con_vars) - sum(repeated_con_vars)
    m.addConstr(sum_cons == 2*l)  # 常数个数等于2*l，则平均有1个解

    # 给定变量个数
    all_vars = [var for var in m.getVars() if "var_" in var.VarName]
    m.addConstr(sum(all_vars) == numVar)

    # 目标函数
    # 方程次数之和
    all_equs = [var for var in m.getVars() if "equ" in var.VarName]
    obj = m.addVar(vtype=GRB.INTEGER, name="obj")
    m.addConstr(obj == sum(all_equs))
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()

    # 模型设置
    #m.setParam("TimeLimit", 200)# 设定跑200秒后就停止，输出当前找到的最优解
    #m.setParam("BestObjStop", 32) #设定目标函数低于100时就停止
    m.write("Anemoi_CICO_MILP_{}r_{}l_{}v.lp".format(R, l, numVar))
    m.optimize()
    m.write("Anemoi_CICO_MILP_{}r_{}l_{}v.sol".format(R, l, numVar))
    parsing_Anemoi_solution(m, R, l)
    print("分支数：{}，".format(2*l))
    print('轮数：{}，'.format(R))
    print('变量个数：{}，'.format(numVar))
    print('方程次数之和：{}'.format(obj.X))
    print('求解复杂度：2 ^ {}'.format(log2(comb(int(obj.X) + 1, numVar) ** 2)))

def print_out_vector(m, l, vname):
    print(vname+"_R", end = " ")
    for j in range(2*l):
        current_var = m.getVarByName(vname+"_R_{}".format(j))
        print(round(current_var.X), end = " ")
    print()
def parsing_Anemoi_solution(m, R, l):
    # 输入变量
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

    # 输出
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
    # 输入
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
    # 输出
    print(vname+'_{}r'.format(r), end=" ")
    for j in range(l):
        print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, out1)).X), end=" ")
    print(end=" ")
    for j in range(l):
        print(round(m.getVarByName(vname+"_{}r_{}th_{}".format(r, j, out2)).X), end=" ")
    print()
def parsing_Anemoi_solution_v0(m, R, l):
    # 每轮
    for r in range(R):
        print('--------------r = {}-----------------'.format(r))
        print('-----r = {}-大MDS周围----------------'.format(r))
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
