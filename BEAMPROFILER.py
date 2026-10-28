# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 11:13:06 2026

@author: elois
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
#from LinearScanGathering import ExecuteandGather, ReadGathering

#BEAM AXIS FINDER

#takes inputs
beam_radius = c = 10
focal_length = b = 40
lense_radius = a = 15

limits = np.array([[-100, 100],[-75, 75],[-120, 120]])

focus_centre = np.array([10,10, 60]) #these to be found experimentally
far_centre = np.array([15,25,120])

beam_axis = far_centre - focus_centre
axis_length = np.linalg.norm(beam_axis)
axis_unit = beam_axis/ axis_length

'''
I want to define an function that will find the s values of the planes that first touch the limits in each direction
'''

def s_value_intersect(limit,centre,axis,a,b,c):
    A = axis**2 - (1 - axis**2) * a**2 / b**2
    B = -2 * (limit-centre) * axis
    C = (limit-centre)**2 - (1 - axis**2) * c**2
    return np.roots([A, B, C])

def s_limit(limits, focus_centre, axis_unit,a,b,c):

    results = []

    names = [
        ["x_min", "x_max"],
        ["y_min", "y_max"],
        ["z_min", "z_max"]]

    for i in range(3):
        for j in range(2):
            roots = s_value_intersect(
                limits[i][j],
                focus_centre[i],
                axis_unit[i],
                a, b, c)

            results.append({
                "plane": names[i][j],
                "limit": limits[i][j],
                "s_values": roots})
            
    negative_candidates = []
    positive_candidates = []
    
    for result in results:
        for s in result["s_values"]:
    
            if s < 0:
                negative_candidates.append(
                    (s, result["plane"]))
    
            elif s > 0:
                positive_candidates.append(
                    (s, result["plane"]))

    negative_limit = max(negative_candidates)
    positive_limit = min(positive_candidates)
    
    return negative_limit, positive_limit

#print(s_limit(limits, focus_centre, axis_unit, a, b, c))


min_s = np.trunc(s_limit(limits, focus_centre, axis_unit, a, b, c)[0][0])
max_s = np.trunc(s_limit(limits, focus_centre, axis_unit, a, b, c)[1][0])
'''Is rounding to the integer closest to zero enough of a margin for error?'''
nslices = 20

s_values = np.linspace(min_s,max_s, nslices)

'''convert from s back to x y and z coordinates'''

focus_separation = 1
revolutions = (beam_radius / focus_separation)

def perpendicular_basis(axis):
    axis = axis / np.linalg.norm(axis)

    # Pick a vector that isn't parallel to the axis
    if abs(axis[0]) < 0.9:
        temp = np.array([1.0, 0.0, 0.0])
    else:
        temp = np.array([0.0, 1.0, 0.0])

    u = np.cross(axis, temp)
    u /= np.linalg.norm(u)

    v = np.cross(axis, u)
    v /= np.linalg.norm(v)

    return u, v


#hyperboloid surface: (x-x_0)**2+(y-y_0)**2+(z-z_0)**2 -(1+a**2/b**2)*(axis_unit[0]*(x-x_0)+axis_unit[1](y-y_0)+axis_unit[2](z-z_0))**2= c**2


beam_centre = []
for i in range(int(nslices)):
    beam_centre.append(focus_centre + s_values[i]* axis_unit) #This is an array
  
e1 = perpendicular_basis(axis_unit)[0]
e2 = perpendicular_basis(axis_unit)[1]


'''
#############################################################
origin = np.array([0, 0, 0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the three unit vectors
ax.quiver(
    *origin,
    *axis_unit,
    color='r',
    arrow_length_ratio=0.15,
    label='Beam axis'
)

ax.quiver(
    *origin,
    *e1,
    color='b',
    arrow_length_ratio=0.15,
    label='XZ perpendicular'
)

ax.quiver(
    *origin,
    *e2,
    color='g',
    arrow_length_ratio=0.15,
    label='Perpendicular 2'
)

# Set equal scaling
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.legend()

plt.show()

def axi_radius(L):
    return (L+240)*0.09805 + 10
'''


def hyp_radius(s):
    return c * np.sqrt(1 + (s**2 * a**2) / (b**2 * c**2))

trajectory_dict = {}

for i,s_val in enumerate(s_values):
    
    radius_end = hyp_radius(s_val)
    separation = radius_end/ revolutions
    npoints = max(2, int(revolutions*100))
    theta_values = np.linspace(0,revolutions * 2 * np.pi, npoints)
    radius_values = np.linspace(0, radius_end, npoints)
    trajectory = []
    
    for theta, r in zip(theta_values, radius_values):
        point = (beam_centre[i] + r * np.cos(theta) * e2 + r * np.sin(theta) * e1)
        trajectory.append(point)
        
    def intergrand(x): #problem!!!
        return np.sqrt(np.square(separation*x) + np.square(separation))
    
    trajectory = np.array(trajectory)
    displacements = trajectory - beam_centre[i]
    
    trajectory_dict[i] = {
    's': s_val,
    'centre': beam_centre[i],
    'radius': radius_end,
    'revolutions': revolutions,
    'npoints': npoints,
    'trajectory': trajectory, 
    'text': "\n".join(f"{x:.8f}, {y:.8f}, {z:.8f}" for x, y, z in displacements),
    'length': quad(intergrand, 0, revolutions * 2 * np.pi) } #returning an incorrect value!

#print(trajectory_dict[0]['radius'], trajectory_dict[0]['revolutions'], trajectory_dict[0]['length'][0])

trajectory_0 = trajectory_dict[0]
print(trajectory_0)
#relative = trajectory - beam_centre[i]
#print(np.max(np.abs(relative @ axis_unit)))

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

for i in trajectory_dict:
    trajectory = trajectory_dict[i]['trajectory']

    # Plot as Z, Y, X
    ax.plot(
        trajectory[:, 2],   # screen horizontal = Z
        trajectory[:, 1],   # screen vertical   = Y
        trajectory[:, 0]    # depth             = X
    )

# Beam axis
ax.plot(
    [focus_centre[2], far_centre[2]],
    [focus_centre[1], far_centre[1]],
    [focus_centre[0], far_centre[0]],
    linestyle='--',
    label='Beam axis'
)

# Focus
ax.scatter(
    focus_centre[2],
    focus_centre[1],
    focus_centre[0],
    label='Focus'
)

# Far centre
ax.scatter(
    far_centre[2],
    far_centre[1],
    far_centre[0],
    label='Far centre'
)

ax.set_xlabel('Z')
ax.set_ylabel('Y')
ax.set_zlabel('X')

centre = np.array([0, 0, 0])
half_range = 150

ax.set_xlim(-half_range, half_range)
ax.set_ylim(-half_range, half_range)
ax.set_zlim(-half_range, half_range)

ax.set_box_aspect([1, 1, 1])

ax.legend()
plt.show()

'''
#LIMITS CONSIDERATION
x_limits_adj = x_limits - focus_centre[0]
y_limits_adj = y_limits - focus_centre[1]
z_limits_adj = z_limits - focus_centre[2]

print(x_limits_adj, y_limits_adj, z_limits_adj)

x_result = min(abs(x_limits_adj[0]), abs(x_limits_adj[1]))
y_result = min(abs(y_limits_adj[0]), abs(y_limits_adj[1]))
z_result = min(abs(z_limits_adj[0]), abs(z_limits_adj[1]))

#this gives furthest points that can be travelled from the focus point

max_radius = min(x_result, y_result)

# this gives the maximum radius of a spiral that can be traced in the x and y plane centred at the z co-ordinate of the focus

z = np.linspace(z_limits_adj[0], z_limits_adj[1], 20)
r = hyp_radius(z)

plt.scatter(z,r)
plt.show

print(max(r))
    
if not (x_limits_adj[0] <= max(r) <= x_limits_adj[1] and
        y_limits_adj[0] <= max(r) <= y_limits_adj[1]):
    raise ValueError(
        f"Maximum radius value {max(r)} lies outside the limits from focus: "
        f"x={x_limits}, y={y_limits}"
    )
 
print(len(z))

#moves to centre of focus

#executes spiral at centre and then spiral at furthest point without pinhole (4 or 5mm spacing)

#plots gathering data and locates two 'centres' coordinates
'''

