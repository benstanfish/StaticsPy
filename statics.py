"""Library of functions for structural beam statics applications."""
__author__ = "Ben Fisher"
__version__ = "0.0.1"
__license__ = "GPL"
__credits__ = ["Ben Fisher"]
__status__ = "Development"
__maintainer__ = "Ben Fisher"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uuid

from math import sqrt, copysign
from tkinter.filedialog import askopenfilename, asksaveasfilename

default_station_count = 101
stations = np.linspace(0.0,1.0,default_station_count)
stations

class Beam:
    def __init__(self, span_length):
        self.id = str(uuid.uuid4())
        self.name = "default_beam_name"
        self.boundaries = np.array([0,0])     # -1 is cantilevered, 0 is pinned, 1 is fixed
        self.span_length = span_length
        self.default_stations = span_length * np.copy(stations)
        self.added_stations = np.array([])
        self.all_stations = self.default_stations
        self.load_types = []
        self.load_params = []
        self.messages = []
        self.V = []  # See discussion below - this needs to be list to hold arrays of vectors
        self.M = []
        self.Total_V = np.array([])
        self.Total_V = np.array([])

    def Add_Stations(self, new_stations):
        for x in new_stations:
            if (x >= 0) & (x <= self.span_length) & (x not in self.added_stations):
                self.added_stations = np.append(self.added_stations, x)
        self.Update_Stations()
        
    def Remove_Stations(self, remove_stations):
        for x in remove_stations:
            if (x >= 0) & (x <= self.span_length) & (x in self.added_stations):
                self.added_stations = np.setdiff1d(self.added_stations, np.array(remove_stations))
        self.Update_Stations()
    
    def Update_Stations(self):
        self.all_stations = np.union1d(self.default_stations, self.added_stations)

    def Append_Load_Type(self, args):
        self.load_types.append(args)

    def Append_Load_Params(self, args):
        self.load_params.append(args)

    def Append_Shears(self, args):
        self.V.append(args)
        
    def Append_Moments(self, args):
        self.M.append(args)

    def Build_Loads(self):
        self.V.clear()
        self.M.clear()
        for j in range(0,len(self.load_types)):
            m = globals()[self.load_types[j]]
            func = getattr(m,'Get_Load_Effects')
            func(self, self.load_params[j])
        
    def Combine_Loads(self):
        self.Total_V = sum(self.V)
        self.Total_M = sum(self.M)
        
    def Show_Shear(self):
        x = self.all_stations
        y = self.Total_V
        zeros = np.zeros(self.all_stations.size)
        fig, ax = plt.subplots()
        ax.set_xlabel("Distance along Span")
        ax.set_ylabel("Shear")
        plt.plot(x,np.zeros(x.size),color="lightgray")
        for i in range(0,len(self.load_types)):
            plt.fill_between(x,zeros,self.V[i],color="skyblue",alpha=0.25)
        plt.plot(x,y,color="dodgerblue",linewidth=2)
        plt.savefig("{}-shear.png".format(self.name),dpi=300,pad_inches=0.1)

    def Show_Moment(self):
        x = self.all_stations
        y = self.Total_M
        zeros = np.zeros(self.all_stations.size)
        fig, ax = plt.subplots()
        ax.set_xlabel("Distance along Span")
        ax.set_ylabel("Moment")
        plt.plot(x,np.zeros(x.size),color="lightgray")
        for i in range(0,len(self.load_types)):
            plt.fill_between(x,zeros,self.M[i],color="skyblue",alpha=0.25)
        plt.plot(x,y,color="dodgerblue",linewidth=2)
        plt.savefig("{}-moment.png".format(self.name),dpi=300,pad_inches=0.1)



class Simple_Point:
    """
        Static class that provides functions relating to a single concentrated load at a point x = a, measured from the left beam support.
    """
    def Add_Load(beam: Beam, magnitude, location):
        if (beam.boundaries[0] == 0) & (beam.boundaries[1] == 0):   # Prevents registering this load on a non-pin pin beam.      
            beam.Add_Stations([location])
            beam.Append_Load_Type('Simple_Point')
            beam.Append_Load_Params([magnitude, location])

    def Get_Load_Effects(beam: Beam, args):
        a = args[0]
        P = args[1]
        locs = np.copy(beam.all_stations)
        shears = np.zeros(locs.size)
        moments = np.zeros(locs.size)
        length = locs[np.argmax(locs)]
        b = length - a
        shears.fill(P/length)
        for i in range(0,np.argmax(locs) + 1):
            if locs[i] <= a:
                shears[i] *= b
            else:
                shears[i] *= -a
        for i in range(0,np.argmax(locs) + 1):
            mmax =  P * b * a / length
            if locs[i] <= a:
                moments[i] = P * b * locs[i] / length
            else:
                moments[i] = P * a * (length - locs[i]) / length
        beam.Append_Shears(shears)
        beam.Append_Moments(moments)
        
        
        
class Simple_UDL:
    """
        Static class that provides functions relating to a uniform distribued load from x = a to x = a + b distance from the left support.
    """
    def Add_Load(beam: Beam, magnitude, a = 0, b):
        if (beam.boundaries[0] == 0) & (beam.boundaries[1] == 0):   # Prevents registering this load on a non-pin pin beam.      
            beam.Add_Stations([a,b])
            beam.Append_Load_Type('Simple_UDL')
            beam.Append_Load_Params([magnitude, a, b])

    def Get_Load_Effects(beam: Beam, args):
        w = args[0]
        a = args[1]
        b = args[2]
        locs = np.copy(beam.all_stations)
        shears = np.zeros(locs.size)
        moments = np.zeros(locs.size)
        length = locs[np.argmax(locs)]
        c = length - a - b
        R_left = w * b / 2 / length * (2 * c + b)
        R_right = w * b / 2 / length * (2 * a + b)
        for i in range(0,np.argmax(locs) + 1):
            if locs[i] <= a:
                shears[i] = R_left
            elif locs[i] >= a + b:
                shears[i] = -1 * R_right
            else:
                shears[i] = R_left - w * (w - a)
        for i in range(0,np.argmax(locs) + 1):
            if locs[i] <= a:
                moments[i] = R_left * locs[i]
            elif locs[i] >= a + b:
                moments[i] = R_right * (length - locs[i])
            else:
                moments[i] = R_left * locs[i] - w / 2 * (locs[i]-a) ** 2
        beam.Append_Shears(shears)
        beam.Append_Moments(moments)