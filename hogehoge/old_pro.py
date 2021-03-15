"""使ってないプログラム"""



## グラフで確認 ###
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("X[m]")
ax.set_ylabel("Y[m]")
ax.set_zlabel("Z[m]")

#ax.plot(origin_positions[0:2, 0], origin_positions[0:2, 1], origin_positions[0:2, 2], "-", label = "base")
# ax.plot(origin_positions[3:5, 0], origin_positions[3:5, 1], origin_positions[3:5, 2], "-", label = "shoulder")
# ax.plot(origin_positions[4:7, 0], origin_positions[4:7, 1], origin_positions[4:7, 2], "-", label = "arm1")
# ax.plot(origin_positions[6:9, 0], origin_positions[6:9, 1], origin_positions[6:9, 2], "-", label = "arm2")
# ax.plot(origin_positions[8:10, 0], origin_positions[8:10, 1], origin_positions[8:10, 2], "-", label = "arm3")
# ax.plot(origin_positions[9:11, 0], origin_positions[9:11, 1], origin_positions[9:11, 2], "-", label = "endeffecter")

ax.plot(origin_positions[:, 0], origin_positions[:, 1], origin_positions[:, 2], "-", label = "arm")
ax.scatter(origin_positions[0, 0], origin_positions[0, 1], origin_positions[0, 2], label = "World origin", s = 100)
ax.scatter(origin_positions[1, 0], origin_positions[1, 1], origin_positions[1, 2], label = "Base", s = 100)
ax.scatter(origin_positions[4, 0], origin_positions[4, 1], origin_positions[4, 2], label = "S1(joint2)", s = 100)
ax.scatter(origin_positions[6, 0], origin_positions[6, 1], origin_positions[6, 2], label = "E1(joint4)", s = 100)
ax.scatter(origin_positions[8, 0], origin_positions[8, 1], origin_positions[8, 2], label = "W1(joint6)", s = 100)
ax.scatter(origin_positions[10, 0], origin_positions[10, 1], origin_positions[10, 2], marker = "*", label = "Gripper", s = 100)

ax.scatter(origin_positions[0, 0], origin_positions[0, 1], origin_positions[0, 2], label = "${o_{WL}}$", s = 10)
ax.scatter(origin_positions[1, 0], origin_positions[1, 1], origin_positions[1, 2], label = "${o_{BL}}$", s = 10)
ax.scatter(origin_positions[2, 0], origin_positions[2, 1], origin_positions[2, 2], label = "${o_{0}}$", s = 10)
ax.scatter(origin_positions[3, 0], origin_positions[3, 1], origin_positions[3, 2], label = "${o_{1}}$", s = 10)
ax.scatter(origin_positions[4, 0], origin_positions[4, 1], origin_positions[4, 2], label = "${o_{2}}$", s = 10)
ax.scatter(origin_positions[5, 0], origin_positions[5, 1], origin_positions[5, 2], label = "${o_{3}}$", s = 10)
ax.scatter(origin_positions[6, 0], origin_positions[6, 1], origin_positions[6, 2], label = "${o_{4}}$", s = 10)
ax.scatter(origin_positions[7, 0], origin_positions[7, 1], origin_positions[7, 2], label = "${o_{5}}$", s= 10)
ax.scatter(origin_positions[8, 0], origin_positions[8, 1], origin_positions[8, 2], label = "${o_{6}}$", s = 10)
ax.scatter(origin_positions[9, 0], origin_positions[9, 1], origin_positions[9, 2], label = "${o_{7}}$", s = 10)
ax.scatter(origin_positions[10, 0], origin_positions[10, 1], origin_positions[10, 2], marker = "*", label = "${o_{GL}}$", s = 10)

ax.text(origin_positions[0, 0], origin_positions[0, 1], origin_positions[0, 2], "${o_{WL}}$", fontsize = 10)
ax.text(origin_positions[1, 0], origin_positions[1, 1], origin_positions[1, 2], "${o_{BL}}$", fontsize = 10)
ax.text(origin_positions[2, 0], origin_positions[2, 1], origin_positions[2, 2], "${o_{0}}$", fontsize = 10)
ax.text(origin_positions[3, 0], origin_positions[3, 1], origin_positions[3, 2], "${o_{1}}$", fontsize = 10)
ax.text(origin_positions[4, 0], origin_positions[4, 1], origin_positions[4, 2], "${o_{2}}$", fontsize = 10)
ax.text(origin_positions[5, 0], origin_positions[5, 1], origin_positions[5, 2], "${o_{3}}$", fontsize = 10)
ax.text(origin_positions[6, 0], origin_positions[6, 1], origin_positions[6, 2], "${o_{4}}$", fontsize = 10)
ax.text(origin_positions[7, 0], origin_positions[7, 1], origin_positions[7, 2], "${o_{5}}$", fontsize = 10)
ax.text(origin_positions[8, 0], origin_positions[8, 1], origin_positions[8, 2], "${o_{6}}$", fontsize = 10)
ax.text(origin_positions[9, 0], origin_positions[9, 1], origin_positions[9, 2], "${o_{7}}$", fontsize = 10)
ax.text(origin_positions[10, 0], origin_positions[10, 1], origin_positions[10, 2], "${o_{GL}}$", fontsize = 10)

ax.legend()

# 三軸のスケールを揃える
max_range = np.array([origin_positions[:, 0:1].max()-origin_positions[:, 0:1].min(), \
    origin_positions[:, 1:2].max()-origin_positions[:, 1:2].min(), \
        origin_positions[:, 2:3].max()-origin_positions[:, 2:3].min()]).max() * 0.5
mid_x = (origin_positions[:, 0:1].max()+origin_positions[:, 0:1].min()) * 0.5
mid_y = (origin_positions[:, 1:2].max()+origin_positions[:, 1:2].min()) * 0.5
mid_z = (origin_positions[:, 2:3].max()+origin_positions[:, 2:3].min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()

