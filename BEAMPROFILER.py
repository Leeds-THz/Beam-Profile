# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 11:13:06 2026

@author: elois
"""
import XPS_C8_drivers
import numpy as np
from newportxps import NewportXPS
import matplotlib.pyplot as plt
#import requests
import os
#import csv
import ftplib
from datetime import datetime
from XPSHelper import GetGatheringFile, InitXPS
from scipy.interpolate import griddata
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


xps = InitXPS('192.168.0.254')  

name = 'spiral_scan'
filename = f"{name}.txt"



def GenerateSplineTrajectory(filename, separation, radius, upload=True, plot=True):
    
    revolutions = radius/separation
    npoints = int(revolutions*100)
    
    theta = np.linspace(0, revolutions*2*np.pi, npoints)
    radius_values = np.linspace(0, radius, npoints)
    
    positions = {
        'X': radius_values*np.cos(theta),
        'Y': radius_values*np.sin(theta)}
    
    text_lines = []

    for i in range(npoints-1):
        x = positions['X'][i]
        y = positions['Y'][i]
        text_lines.append(f"{x:.8f}, {y:.8f}, {0.00}")

    text = "\n".join(text_lines)
    
    first_point = np.array([
        positions['X'][0],
        positions['Y'][0],
        0.00
    ])
    
    if upload:
        xps.upload_trajectory(filename, text)
        print(f"Uploaded : {filename}")
        
    if plot:
        plt.close('all')
        plt.scatter(positions['X'], positions['Y'])
        plt.xlabel("Horizontal position (mm)")
        plt.ylabel("Vertical position (mm)")
        plt.title("Spiral Trajectory")
        plt.axis('equal')
        plt.show()
    
    return first_point


def PrepareTrajectory(start_position):
    err, ret = xps._xps.GroupMoveAbsolute(
        xps._sid,
        'XYZ1',
        start_position)

    if err != 0:
        raise RuntimeError(f'GroupMoveAbsolute error: {ret}')
        
        
def ExecuteandGather(filename, timestep, velocity, start_position):
    PrepareTrajectory(start_position)

    err, ret = xps._xps.GatheringReset(xps._sid)
    if err != 0:
        raise RuntimeError(f'GatheringReset error: {ret}')
        
    err, ret = xps._xps.GatheringConfigurationSet(xps._sid,
        ['GatheringTime','XYZ1.X.CurrentPosition',
        'XYZ1.Y.CurrentPosition', 'XYZ1.Z.CurrentPosition','GPIO4.ADC1'])
    if err != 0:
        raise RuntimeError(f'GatheringConfigurationSet error: {ret}')
        
    err, ret = xps._xps.EventExtendedConfigurationTriggerSet(xps._sid,
        ['XYZ1.Spline.TrajectoryStart'],['0'],['0'],['0'],['0'])
    if err != 0:
        raise RuntimeError(f'EventExtendedConfigurationTriggerSet error: {ret}')

    err, ret = xps._xps.EventExtendedConfigurationActionSet(xps._sid,
        ['GatheringRun'], ['8000'], [timestep], ['0'], ['0'])
    if err != 0:
        raise RuntimeError(f'EventExtendedConfigurationActionSet error: {ret}')

    err, ret = xps._xps.EventExtendedStart(xps._sid)
    if err != 0:
        raise RuntimeError(f'EventExtendedStart error: {ret}')
        
    err, ret = xps._xps.XYZSplineExecution(xps._sid,
        'XYZ1',filename,velocity,100)
    if err != 0:
        raise RuntimeError(f'XYZSplineExecution error: {ret}')

    print(xps._xps.GatheringCurrentNumberGet(xps._sid))
    
    err, ret = xps._xps.GatheringStopAndSave(xps._sid)
    if err != 0:
        raise RuntimeError(f'GatheringStopAndSave error: {ret}')

    xps._xps.EventExtendedRemove(xps._sid, '0')
    
    
def ReadGathering(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    header_tokens = []
    data_start_idx = None
    for i, line in enumerate(lines[1:], start=1):
        tokens = line.split()
        if not tokens:
            continue
        try:
            float(tokens[0])
            data_start_idx = i
            break
        except ValueError:
            header_tokens.extend(tokens)

    columns = header_tokens
    data = {col: [] for col in columns}

    for line in lines[data_start_idx:]:
        tokens = line.split()
        if not tokens:
            continue
        for col, val in zip(columns, tokens):
            data[col].append(float(val))

    return data
    
    
def FindBeamAxis(limits, approx_focus, waist_radius, far_radius):
    """
    Performs one spiral scan at the approximate focus and one at the far
    z limit. Finds the point of maximum intensity from each scan and
    returns the two points defining the beam axis.
    """

    approx_focus = np.array(approx_focus, dtype=float)

    far_focus = np.array([
        approx_focus[0],
        approx_focus[1],
        limits['z'][0] + 1])

    safe_limits = {
        'x': np.array([limits['x'][0] + 1, limits['x'][1] - 1]),
        'y': np.array([limits['y'][0] + 1, limits['y'][1] - 1]),
        'z': np.array([limits['z'][0] + 1, limits['z'][1] - 1])}

    if approx_focus[0] < safe_limits['x'][0] or approx_focus[0] > safe_limits['x'][1]:
        raise ValueError('Approximate focus X position is outside the safe stage limits')

    if approx_focus[1] < safe_limits['y'][0] or approx_focus[1] > safe_limits['y'][1]:
        raise ValueError('Approximate focus Y position is outside the safe stage limits')

    if approx_focus[2] < safe_limits['z'][0] or approx_focus[2] > safe_limits['z'][1]:
        raise ValueError('Approximate focus Z position is outside the safe stage limits')

    if far_focus[0] < safe_limits['x'][0] or far_focus[0] > safe_limits['x'][1]:
        raise ValueError('Far focus X position is outside the safe stage limits')

    if far_focus[1] < safe_limits['y'][0] or far_focus[1] > safe_limits['y'][1]:
        raise ValueError('Far focus Y position is outside the safe stage limits')

    if far_focus[2] < safe_limits['z'][0] or far_focus[2] > safe_limits['z'][1]:
        raise ValueError('Far focus Z position is outside the safe stage limits')

    scan_results = {}

    scans = [
        ('focus', approx_focus, waist_radius, 'focus_spiral.txt'),
        ('far', far_focus, far_radius, 'far_spiral.txt')]

    for name, centre, radius, filename in scans:

        print(f'Running {name} spiral')
        print(f'Centre: {centre}')
        print(f'Radius: {radius}')

        first_point = GenerateSplineTrajectory(
            filename,
            1,
            radius,
            upload=True,
            plot=False)

        ExecuteandGather(
            filename,
            '2000',
            8,
            centre)

        GetGatheringFile(xps)

        gathering_file = 'Gathering.dat'
        saved_file = f'{name}_gathering.dat'

        if os.path.exists(saved_file):
            os.remove(saved_file)

        os.rename(gathering_file, saved_file)

        result = ReadGathering(saved_file)

        X = np.array(result['XYZ1.X.CurrentPosition'])
        Y = np.array(result['XYZ1.Y.CurrentPosition'])
        Z = np.array(result['XYZ1.Z.CurrentPosition'])
        values = np.array(result['GPIO4.ADC1'])

        max_index = np.argmax(values)

        maximum_position = np.array([
            X[max_index],
            Y[max_index],
            Z[max_index]
        ])

        maximum_value = values[max_index]

        scan_results[name] = {
            'X': X,
            'Y': Y,
            'Z': Z,
            'values': values,
            'maximum_position': maximum_position,
            'maximum_value': maximum_value,
            'file': saved_file
        }

        print(f'Maximum intensity: {maximum_value}')
        print(f'Position: {maximum_position}')

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for column, name in enumerate(['focus', 'far']):

        result = scan_results[name]

        X = result['X']
        Y = result['Y']
        values = result['values']
        maximum_position = result['maximum_position']

        axes[0, column].scatter(
            X,
            Y,
            c=values,
            s=5
        )

        axes[0, column].scatter(
            maximum_position[0],
            maximum_position[1],
            marker='x',
            s=100
        )

        axes[0, column].set_xlabel("Horizontal position (mm)")
        axes[0, column].set_ylabel("Vertical position (mm)")
        axes[0, column].set_title(f'{name.capitalize()} scan - data points')
        axes[0, column].axis('equal')

        x = np.linspace(X.min(), X.max(), 500)
        y = np.linspace(Y.min(), Y.max(), 500)

        xg, yg = np.meshgrid(x, y)

        zg = griddata(
            (X, Y),
            values,
            (xg, yg),
            method='cubic'
        )

        axes[1, column].pcolormesh(
            xg,
            yg,
            zg,
            shading='auto'
        )

        axes[1, column].scatter(
            maximum_position[0],
            maximum_position[1],
            marker='x',
            s=100
        )

        axes[1, column].set_xlabel("Horizontal position (mm)")
        axes[1, column].set_ylabel("Vertical position (mm)")
        axes[1, column].set_title(f'{name.capitalize()} scan - interpolated data')
        axes[1, column].axis('equal')

    plt.tight_layout()
    plt.show()

    beam_axis = {
        'focus': scan_results['focus']['maximum_position'],
        'far': scan_results['far']['maximum_position']
    }

    print()
    print(f'Focus position: {beam_axis["focus"]}')
    print(f'Far position: {beam_axis["far"]}')

    return beam_axis


GenerateSplineTrajectory('focus_spiral.txt', 4,40)
ExecuteandGather('focus_spiral.txt', 2000, 8, [-15,-5,120])

GetGatheringFile(xps)

gathering_file = 'Gathering.dat'
saved_file = f'{name}_gathering.dat'

if os.path.exists(saved_file):
    os.remove(saved_file)

os.rename(gathering_file, saved_file)

result = ReadGathering(saved_file)

X = np.array(result['XYZ1.X.CurrentPosition'])
Y = np.array(result['XYZ1.Y.CurrentPosition'])
Z = np.array(result['XYZ1.Z.CurrentPosition'])
values = np.array(result['GPIO4.ADC1'])

max_index = np.argmax(values)

maximum_position = np.array([
    X[max_index],
    Y[max_index],
    Z[max_index]
])

maximum_value = values[max_index]

scan_results = {
    'X': X,
    'Y': Y,
    'Z': Z,
    'values': values,
    'maximum_position': maximum_position,
    'maximum_value': maximum_value,
    'file': saved_file
}

print(f'Maximum intensity: {maximum_value}')
print(f'Position: {maximum_position}')

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

