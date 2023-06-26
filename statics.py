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
from tkinter.filedialog import askopenfilename

default_station_count = 101
stations = np.linspace(0.0,1.0,default_station_count)
stations

class Beam:
    def __init__(self, span_length):
        self.id = str(uuid.uuid4())
        self.name = "default_name"
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
        for i in range(0,len(self.load_types)):
            m = globals()[self.load_types[i]]
            func = getattr(m,'Create_Load_Vectors')
            func(self, self.load_params[i])
        
    def Combine_Loads(self):
        self.Total_V = sum(self.V)
        self.Total_M = sum(self.M)


class Simple_Point:
    """
        Static class that provides functions relating to a single concentrated load at a point x = a, measured from the left beam support.
    """
    def Register_Load(beam: Beam, location, magnitude):
        if (beam.boundaries[0] == 0) & (beam.boundaries[1] == 0):   # Prevents registering this load on a non-pin pin beam.      
            beam.Add_Stations([location])
            beam.Append_Load_Type(Simple_Point)
            beam.Append_Load_Params([location, magnitude])

    def Create_Load_Vectors(beam: Beam, args):
        a = args[0]
        P = args[1]
        locs = np.copy(beam.all_stations)
        shears = np.zeros(locs.size)
        moments = np.zeros(locs.size)
        length = locs[np.argmax(locs)]
        b = length - a
        length
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