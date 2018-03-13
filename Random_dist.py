import numpy as np
import random, time

pos = []
vel = []
sigma = 3.4
box_len = 34.778
num_of_particles = 864
cutoff = 7.65
kb = 1.381  #10 -23
mass = 6.69  #10 -23

# Temp = 94.4 K


# Creating position
#x_pos = np.random.uniform(0, box_len, 10).sort()
#y_pos = np.random.uniform(0, box_len, 10).sort()
#z_pos = np.random.uniform(0, box_len, 10).sort()

x_pos = np.linspace(0.1, 31.220, 10)
y_pos = np.linspace(0.1, 31.220, 10)
z_pos = np.linspace(0.1, 31.220, 10)

for x in x_pos:
    for y in y_pos:
        for z in z_pos:
            dic = {"x":x, "y":y, "z":z}
            pos.append(dic)

# Veifying for distance of sigma
# This can also be use to calculate radial dist function
# Removing extra
while(len(pos) > 864):
    num = random.randint(0,len(pos)-1)
    pos.remove(pos[num])

# Creating velocity
# I used numpy in python to get Normal distribution with V2 summation
x_vel = np.random.normal(138, 30, 864)
y_vel = np.random.normal(138, 30, 864)
z_vel = np.random.normal(138, 30, 864)

for i in xrange(0,864):
    dic = {"vx":x_vel[i], "vy":y_vel[i], "vz":z_vel[i]}
    vel.append(dic)

# Summetion of velocity square needs to be 50509.67246636772 so by the variable I am just checking value of Vel ** 2, For analysys I am using vel**2 sum as 50524 which is pretty close

temperature = 0
for num in range(864):
    temp = vel[num]
    velocity = temp["vx"]**2 + temp["vy"]**2 + temp["vz"]**2
    temperature += velocity

print len(pos), len(vel), temperature

# scaling vel
scale = (50509664/temperature) ** .5

for num in range(864):
    temp = vel[num]
    num = random.randint(0,1)
    if num == 0:
        num = -1
    temp["vx"] = num * scale * temp["vx"]
    num = random.randint(0,1)
    if num == 0:
        num = -1
    temp["vy"] = num * scale * temp["vy"]
    num = random.randint(0,1)
    if num == 0:
        num = -1
    temp["vz"] = num * scale * temp["vz"]

temperature = 0
for num in range(864):
    temp = vel[num]
    velocity = temp["vx"]**2 + temp["vy"]**2 + temp["vz"]**2
    temperature += velocity

print len(pos), len(vel), temperature


time.sleep(5)

with open("data.txt", 'w') as data:
    for num in range(864):
        temp_vel = vel[num]
        temp_pos = pos[num]
        data.write(str(temp_pos["x"]) + " " + str(temp_pos["y"]) + " " + str(temp_pos["z"]) + " " + str(temp_vel["vx"]) + " " + str(temp_vel["vy"]) + " " + str(temp_vel["vz"]) + "\n")
