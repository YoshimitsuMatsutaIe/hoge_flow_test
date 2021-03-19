"""baxterロボットのシミュレーションテスト
・左腕のみ
"""

import numpy as np
from math import cos, sin, tan, pi
#import itertools
import csv
import matplotlib.pyplot as plt
import matplotlib.animation as anm
from mpl_toolkits.mplot3d import Axes3D 
#import matplotlib.patches as patches
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname=r'C:\WINDOWS\Fonts\Arial.ttf', size=14)
import time

# from baxter_utils_2 import *
# from baxter_utils import *
from baxter_utils_3 import *
from rmp import *


### パラメータ ###
BaxterKinema = BaxterKinematics3(L = 278e-3,
                                 h = 64e-3,
                                 H = 1104e-3,
                                 L0 = 270.35e-3,
                                 L1 = 69e-3,
                                 L2 = 364.35e-3,
                                 L3 = 69e-3,
                                 L4 = 374.29e-3,
                                 L5 = 10e-3,
                                 L6 = 368.3e-3)

# joint limits
q1_min, q1_max = -141, 51
q2_min, q2_max = -123, 60
q3_min, q3_max = -173, 173
q4_min, q4_max = -3, 150
q5_min, q5_max = -175, 175
q6_min, q6_max = -90, 120
q7_min, q7_max = -175, 175
q_min = np.array([[q1_min, q2_min, q3_min, q4_min, q5_min, q6_min, q7_min]]).T * pi / 180
q_max = np.array([[q1_max, q2_max, q3_max, q4_max, q5_max, q6_max, q7_max]]).T * pi / 180


### 初期値 ###
# baxterの初期姿勢（規定）
q1_init = 0
q2_init = -31
q3_init = 0
q4_init = 43
q5_init = 0
q6_init = 72
q7_init = 0
q = np.array([[q1_init, q2_init, q3_init, q4_init, q5_init, q6_init, q7_init]]).T * pi / 180 # ジョイント角度ベクトル
dq = np.array([[0, 0, 0, 0, 0, 0, 0]]).T  # ジョイント角速度ベクトル
origins = BaxterKinema.origins(q)  # 局所座標系の原点を計算
dorigins = []
for i in range(0, 11):
    dorigins.append(np.zeros((3, 1)))
J_all = BaxterKinema.jacobi_all(q)
dJ_all = BaxterKinema.djacobi_all(q, dq)

# print("origins = ", origins)
# print("J_all = ", J_all)
# print("dJ_all = ", dJ_all)

### 初期値格納 ###
q_init_list = (q.ravel()).tolist()
dq_init_list = (dq.ravel()).tolist()
origins_init_list = []
for i in range(0, 11):
    origins_init_list.extend((origins[i].ravel()).tolist())


q_temp, dq_temp, origins_temp = [0], [0], [0]

q_temp.extend(list(q_init_list))
dq_temp.extend(list(dq_init_list))
origins_temp.extend(list(origins_init_list))

q_his = [["time", "q_1", "q_2", "q_3", "q_4", "q_5", "q_6", "q_7"]]
dq_his = [["time", "dq_1", "dq_2", "dq_3", "dq_4", "dq_5", "dq_6", "dq_7"]]
origins_his = [["time", \
    "origin_Wo_x", "origin_Wo_y", "origin_Wo_z", \
        "origin_BL_x", "origin_BL_y", "origin_BL_z", \
            "origin_0_x", "origin_0_y", "origin_0_z", \
                "origin_1_x", "origin_1_y", "origin_1_z", \
                    "origin_2_x", "origin_2_y", "origin_2_z", \
                        "origin_3_x", "origin_3_y", "origin_3_z", \
                            "origin_4_x", "origin_4_y", "origin_4_z", \
                                "origin_5_x", "origin_5_y", "origin_5_z", \
                                    "origin_6_x", "origin_6_y", "origin_6_z", \
                                        "origin_7_x", "origin_7_y", "origin_7_z", \
                                            "origin_GL_x", "origin_GL_y", "origin_GL_z"]]

q_his.append(q_temp)
dq_his.append(dq_temp)
origins_his.append(origins_temp)


time_span = 10  # シミュレーション時間[sec]
time_interval = 0.1  # 刻み時間[sec]

# 目標位置
goal_posi = np.array([[0.2, -0.6, 1]])

# 障害物位置
#obs_posi = np.array([[0.4, -0.6, 1]])
obs_posi = np.array([[0.6, -0.6, 1],
                     [0.6, -0.6, 1.05],
                     [0.6, -0.6, 0.95],
                     [0.6, -0.6, 1.10]])
oend = obs_posi.shape[0]

time_sim_start = time.time()

# RMP豆乳
RMP = RMP1(attract_max_speed = 1, 
           attract_gain = 10, 
           attract_a_damp_r = 0.3,
           attract_sigma_W = 1, 
           attract_sigma_H = 1, 
           attract_A_damp_r = 0.3, 
           obs_scale_rep = 0.2,
           obs_scale_damp = 1,
           obs_ratio = 0.5,
           obs_rep_gain = 0.5,
           obs_r = 15,
           jl_gamma_p = 0.05,
           jl_gamma_d = 0.1,
           jl_lambda = 0.7)

RMPGDS = RMP2(attract_max_speed = 0.1, 
              attract_gain = 5,
              attract_alpha_f = 0.3,
              attract_sigma_alpha = 1,
              attract_sigma_gamma = 1,
              attract_w_u = 10,
              attract_w_l = 1,
              attract_alpha = 0.1,
              attract_epsilon = 1e-5)

result = []

### シミュレーション本体 ###
for t in np.arange(time_interval, time_span + time_interval, time_interval):
    
    if np.linalg.norm(goal_posi.T - origins[10]) < 5e-3:  #手先誤差5mm以下になったら成功
        print("目標到達！")
        result.append("succes!")
        break
    
    if (q < q_min).any() or (q_max < q).any():  # ジョイント制限判定
        print("ジョイント制限を突破")
        result.append("exceed angle limit")
        break
    
    else:
        pull_f_all = []
        pull_M_all = []
        
        for i in range(10, 11, 1):
            ## RMP計算
            for j in range(0, oend, 1):
                # 障害物回避あり
                a = RMP.a_obs(origins[i], dorigins[i], obs_posi[j:j+1, :].T)
                M = RMP.metric_obs(origins[i], dorigins[i], obs_posi[j:j+1].T, a)
                f = M @ a
                
                # # 障害物回避なし
                # f = np.zeros((3, 1))
                # M = np.zeros((3, 3))
                
                ## pull演算
                J = J_all[i]
                dJ = dJ_all[i]
                pull_f = J.T @ (f - M @ dJ @ dq)
                pull_M = J.T @ M @ J
                pull_f_all.append(pull_f)
                pull_M_all.append(pull_M)
            
            if i == 10:
                # # RMP1
                # a_GL = RMP.a_attract(origins[10], dorigins[10], goal_posi.T)
                # M_GL = RMP.metric_attract(origins[10], dorigins[10], goal_posi.T, a_GL)
                # f_GL = M_GL @ a_GL
                
                # RMP2
                M_GL = RMPGDS.M_attract(origins[10], dorigins[10], goal_posi.T, np.zeros((3, 1)))
                f_GL = RMPGDS.f_attract(origins[10], dorigins[10], goal_posi.T, np.zeros((3, 1)), M_GL)
                
                pull_f = J.T @ (f_GL - M_GL @ dJ @ dq)
                pull_M = J.T @ M_GL @ J
                pull_f_all.append(pull_f)
                pull_M_all.append(pull_M)
        
        pull_f_all = np.sum(pull_f_all, axis = 0)
        pull_M_all = np.sum(pull_M_all, axis = 0)
        
        # ジョイント制限処理RMPを配置空間で追加
        a_jl = RMP.a_joint_limit(q, dq, q_min, q_max)
        M_jl = RMP.metric_joint_limit(q)
        f_jl = M_jl @ a_jl
        
        ## resolve演算
        a = np.linalg.pinv(pull_M_all) @ pull_f_all  # 制御指令
        #print("a = ", a)
        
        # オイラー積分する
        
        dq = dq + a * time_interval
        q = q + dq * time_interval
        
        # push演算
        origins = BaxterKinema.origins(q)
        J_all = BaxterKinema.jacobi_all(q)
        dJ_all = BaxterKinema.djacobi_all(q, dq)
        dorigins = []
        for i in range(0, 11):
            dorigins.append(J_all[i] @ dq)
        
        
        # データ格納
        q_list = (q.ravel()).tolist()
        dq_list = (dq.ravel()).tolist()
        origins_list = []
        for i in range(0, 11):
            origins_list.extend((origins[i].ravel()).tolist())
        
        q_temp, dq_temp, origins_temp = [t], [t], [t]
        q_temp.extend(list(q_list))
        dq_temp.extend(list(dq_list))
        origins_temp.extend(list(origins_list))
        
        q_his.append(q_temp)
        dq_his.append(dq_temp)
        origins_his.append(origins_temp)
        
        # ターミナルに計算値を表示（やると遅い）
        # print("t = ", t)
        # print("ee = ", origins[10].T)
        # print("error = ", np.linalg.norm(goal_posi.T - origins[10], ord=2))

tend = t
if len(result) == 0:
    print("時間切れ")
    result.append("timeout")

print("シミュレーション終了")
print("シミュレーション時間", time.time() - time_sim_start)

### データ保存 ###
# f = open('out3.csv', 'w', newline = "")
# writer = csv.writer(f)
# writer.writerows(origins_his)
# f.close()


# ### 静止グラフで確認 ###

# origins_his_T = [list(x) for x in zip(*origins_his)]  # 位置履歴を転置させる

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlabel("X[m]")
# ax.set_ylabel("Y[m]")
# ax.set_zlabel("Z[m]")
# ax.plot(origins_his_T[31][1:], origins_his_T[32][1:], origins_his_T[33][1:], ".-", label = "ee trajectory")
# ax.scatter(goal_posi[0, 0], goal_posi[1, 0], goal_posi[2, 0], marker = "*", label = "goal point")
# ax.scatter(origins[10][0, :], origins[10][1, :], origins[10][2, :], label = "ee last position")

# ax.legend()

# plt.show()


### 動画で確認 ###
time_ani_start = time.time()

origins_his_T = [list(x) for x in zip(*origins_his)]

fig_ani = plt.figure()
ax = fig_ani.gca(projection = '3d')
ax.grid(True)
ax.set_xlabel('X[m]')
ax.set_ylabel('Y[m]')
ax.set_zlabel('Z[m]')


## 三軸のスケールを揃える
# 使用するデータを指定
list_x = []  # x軸配列
list_y = []  # y軸配列
list_z = []  # z軸配列
for i in range(0, 11, 1):
    list_x.extend(origins_his_T[1 + 3 * i][1:])
    list_y.extend(origins_his_T[2 + 3 * i][1:])
    list_z.extend(origins_his_T[3 + 3 * i][1:])
# 軸をセット
max_range = np.array([max(list_x) - min(list_x),
                      max(list_y) - min(list_y),
                      max(list_z) - min(list_z)]).max() * 0.5
mid_x = (max(list_x) + min(list_x)) * 0.5
mid_y = (max(list_y) + min(list_y)) * 0.5
mid_z = (max(list_z) + min(list_z)) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)


# 目標点
ax.scatter(goal_posi[0, 0], goal_posi[0, 1], goal_posi[0, 2],
           s = 100, label = 'goal point', marker = '*', color = '#ff7f00', alpha = 1, linewidths = 1.5, edgecolors = 'red')
ax.scatter(obs_posi[:, 0], obs_posi[:, 1], obs_posi[:, 2],
           s = 100, label = 'obstacle point', marker = '+', color = 'k', alpha = 1)


# アーム全体
body_x, body_y, body_z = [], [], []
for j in range(0, 11, 1):
    # body用の配列作成
    body_x.append(origins_his[1][1 + 3 * j])
    body_y.append(origins_his[1][2 + 3 * j])
    body_z.append(origins_his[1][3 + 3 * j])

bodys = []
bodys.append(ax.plot(body_x, body_y, body_z, "o-", color = "blue")[0])

# グリッパー（エンドエフェクター）軌跡
gl = []
gl.append(ax.plot(origins_his_T[31][1:4], 
                  origins_his_T[32][1:4], 
                  origins_his_T[33][1:4], "-", label = "gl", color = "#ff7f00")[0])

# 局所座標系の原点の番号
t_all = []
name = ("o_Wo", "o_BL", "o_0", "o_1", "o_2", "o_3", "o_4", "o_5", "o_6", "o_7", "o_GL")
for i in range(0, 11, 1):
    t_all.append([ax.text(origins_his_T[3 * i + 1][1],
                          origins_his_T[3 * i + 2][1],
                          origins_his_T[3 * i + 3][1], name[i])])

#ax.legend()

# 時刻表示
timeani = [ax.text(0.8, 0.2, 0.01, "time = 0.0 [s]", size = 10)]
time_template = 'time = %s [s]'

# 結果表示
ax.text(0.8, 0.3, 0.01, result[0], color = "r", size = 14)

ax.set_box_aspect((1,1,1))

def update(i):
    """アニメーションの関数"""
    i = i + 1
    
    body_x, body_y, body_z = [], [], []
    for j in range(0, 11, 1):
        # body用の配列作成
        body_x.append(origins_his[i][1 + 3 * j])
        body_y.append(origins_his[i][2 + 3 * j])
        body_z.append(origins_his[i][3 + 3 * j])
        
        # 原点番号
        t_all[j].pop().remove()
        t_, = [ax.text(origins_his_T[3 * j + 1][i],
                       origins_his_T[3 * j + 2][i],
                       origins_his_T[3 * j + 3][i], name[j])]
        t_all[j].append(t_)
    
    item1 = bodys.pop(0)
    ax.lines.remove(item1)
    bodys.append(ax.plot(body_x, body_y, body_z, "o-", color = "blue")[0])
    
    item2 = gl.pop(0)
    ax.lines.remove(item2)
    gl.append(ax.plot(origins_his_T[31][1:i], 
                      origins_his_T[32][1:i], 
                      origins_his_T[33][1:i], "-", color = "#ff7f00")[0])
    
    # 時刻表示
    timeani.pop().remove()
    timeani_, = [ax.text(0.8, 0.12, 0.01, time_template % (i * time_interval), size = 10)]
    timeani.append(timeani_)
    return None

ani = anm.FuncAnimation(fig = fig_ani, 
                        func = update, 
                        frames = int(tend / time_interval),
                        interval = time_interval * 0.001)

print("アニメ化完了")
print("アニメ化時間", time.time() - time_ani_start)

# ani.save(filename = "hogehooge.gif", 
#             fps = 1 / time_interval, 
#             writer='pillow')



plt.show()








