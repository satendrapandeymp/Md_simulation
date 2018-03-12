import numpy as np
import random

pos = []
vel = []
sigma = 3.4
box_len = 10.229 * sigma
num_of_particles = 864
cutoff = 2.25 * sigma

# Creating position
#x_pos = np.random.uniform(0, box_len, 10).sort()
#y_pos = np.random.uniform(0, box_len, 10).sort()
#z_pos = np.random.uniform(0, box_len, 10).sort()

x_pos = np.linspace(0.1, box_len, 10)
y_pos = np.linspace(0.1, box_len, 10)
z_pos = np.linspace(0.1, box_len, 10)

for x in x_pos:
    for y in y_pos:
        for z in z_pos:
            dic = {"x":x, "y":y, "z":z}
            pos.append(dic)

# Veifying for distance of sigma
# This can also be use to calculate radial dist function
# Removing extra
while(len(pos) > 864):
    num = random.randint(0,len(pos))
    pos.remove(pos[num])

# Creating velocity
x_vel = np.random.normal(10, 4, 864)
y_vel = np.random.normal(10, 4, 864)
z_vel = np.random.normal(20, 4, 864)

for i in xrange(0,864):
    dic = {"vx":x_vel[i], "vy":y_vel[i], "vz":z_vel[i]}
    vel.append(dic)

print len(pos), len(vel)
