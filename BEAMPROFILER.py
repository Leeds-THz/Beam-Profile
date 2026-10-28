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

xps = InitXPS('192.168.0.254') 

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

focus_separation = 1
revolutions = (beam_radius / focus_separation)
nslices = 20

velocity = 12
timestep = '2000'



def s_value_intersect(limit,centre,axis,a,b,c):
    A = axis**2 - (1 - axis**2) * a**2 / b**2
    B = -2 * (limit-centre) * axis
    C = (limit-centre)**2 - (1 - axis**2) * c**2
    return np.roots([A, B, C])

def s_values(limits, focus_centre, axis_unit, nslices, a,b,c):
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

    min_s = np.trunc(negative_limit[0])
    max_s = np.trunc(positive_limit[0])

    s_values = np.linspace(min_s,max_s, nslices)
    
    return s_values

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

def hyp_radius(s):
    return c * np.sqrt(1 + (s**2 * a**2) / (b**2 * c**2))


def trajectory_generator(limits, focus_centre, axis_unit, nslices, a,b,c, plot = True, upload = True):

    def intergrand(x): #problem!!! supposed to give the length of the spiral
        return np.sqrt(np.square(separation*x) + np.square(separation))
    
    values = s_values(limits, focus_centre, axis_unit, nslices, a,b,c)
    
    beam_centre = []
    for i in range(int(nslices)):
        beam_centre.append(focus_centre + values[i]* axis_unit) #This is an array
      
    e1 = perpendicular_basis(axis_unit)[0]
    e2 = perpendicular_basis(axis_unit)[1]

    trajectory_dict = {}
    
    for i,s_val in enumerate(values):
        
        radius_end = hyp_radius(s_val)
        separation = radius_end/ revolutions
        npoints = max(2, int(revolutions*100))
        theta_values = np.linspace(0,revolutions * 2 * np.pi, npoints)
        radius_values = np.linspace(0, radius_end, npoints)
        trajectory = []
        
        for theta, r in zip(theta_values, radius_values):
            point = (beam_centre[i] + r * np.cos(theta) * e2 + r * np.sin(theta) * e1)
            trajectory.append(point)
        
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

    if plot:
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
        
        #centre = np.array([0, 0, 0])
        half_range = 150 ###
        
        ax.set_xlim(-half_range, half_range)
        ax.set_ylim(-half_range, half_range)
        ax.set_zlim(-half_range, half_range)
        
        ax.set_box_aspect([1, 1, 1])
        
        ax.legend()
        plt.show()

    if upload: 
        for i in trajectory_dict:
            xps.upload_trajectory("trajectory{0}".format(i), trajectory_dict[i]['text'])
            
    return trajectory_dict
            

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

def total_trajectory_execute(trajectory_dict):

    for i in trajectory_dict:
        ExecuteandGather("trajectory{0}".format(i), timestep, velocity, trajectory_dict[i]['centre'])
        
        saved_file = "gathering{0}.dat".format(i)
        GetGatheringFile(xps, saved_file)
        
        # gathering_file = 'Gathering.dat'
        # saved_file = "gathering{0}.dat".format(i)
    
        #if os.path.exists(saved_file):
         #   os.remove(saved_file)
    
        # os.rename(gathering_file, saved_file)
    
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
        


if __name__ == "__main__":
    trajectory = trajectory_generator(limits, focus_centre, axis_unit, nslices, a,b,c, plot = False, upload = True)        
    total_trajectory_execute(trajectory)
    

