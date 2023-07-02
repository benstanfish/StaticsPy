"""Library of functions for structural beam statics applications."""
__author__ = "Ben Fisher"
__version__ = "0.0.1"
__license__ = "GPL"
__credits__ = ["Ben Fisher"]
__status__ = "Development"
__maintainer__ = "Ben Fisher"

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uuid
import os
from enum import IntEnum

from math import sqrt, copysign
from tkinter.filedialog import askopenfilename, asksaveasfile

default_station_count = 101
stations = np.linspace(0.0,1.0,default_station_count)
docs_folder = os.path.join(os.path.expanduser('~'),'Documents')
if not os.path.exists(os.path.join(docs_folder,'Statics')):
    os.makedirs(os.path.join(docs_folder,'Statics'))
    save_folder = os.path.join(docs_folder,'Statics')
else:
    save_folder = docs_folder
# create_folder = os.mkdir(os.path.join(os.path.expanduser('~'),'Documents','Statics'))
# save_folder = os.path.expanduser('~\Documents\Statics')

class supportType(IntEnum):
    """Enumeration that represents the boundary condition.
    """
    free = -1
    cant = -1
    cantilever = -1
    pin = 0
    pinned = 0
    roller = 0
    fix = 1
    fixed = 1

class Beam:
    """Creates a single-span beam object.
    """
    def __init__(self, span_length, left_support: IntEnum=supportType.pin, right_support: IntEnum=supportType.pin):
        self.id = str(uuid.uuid4())
        self.name = "default_beam_name"
        self.boundaries = np.array([left_support,right_support])     # -1 is cantilevered, 0 is pinned, 1 is fixed
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
        """Adds stations to the added_stations array.

        Args:
            new_stations (float): location of the station in length units.
        """
        for x in new_stations:
            if (x >= 0) & (x <= self.span_length) & (x not in self.added_stations):
                self.added_stations = np.append(self.added_stations, x)
        self.Update_Stations()
        
    def Remove_Stations(self, remove_stations):
        """Removes the station from the added_stations array.

        Args:
            remove_stations (_type_): _description_
        """
        for x in remove_stations:
            if (x >= 0) & (x <= self.span_length) & (x in self.added_stations):
                self.added_stations = np.setdiff1d(self.added_stations, np.array(remove_stations))
        self.Update_Stations()
    
    def Update_Stations(self):
        """Combines the default_stations and added_stations arrays into the all_stations array.
        """
        self.all_stations = np.union1d(self.default_stations, self.added_stations)

    def Append_Load_Type(self, func_name):
        """Appends func_name to the end of the load_types list.

        Args:
            func_name (str): Name of the load type static class.
        """
        self.load_types.append(func_name)

    def Append_Load_Params(self, load_params):
        """Appends load_params list to the end of the load_params list.

        Args:
            load_params ([]): List of parameters feed from the load class, e.g. load magnitude, a and b distances, etc.
        """
        self.load_params.append(load_params)

    def Append_Shears(self, shears):
        """Append shears vector to the V[] list

        Args:
            shears (float[]): Numpy array of shear values at each point in all_stations.
        """
        self.V.append(shears)
        
    def Append_Moments(self, moments):
        """Append moments vector to the M[] list

        Args:
            moments (float[]): Numpy array of moments values at each point in all_stations.
        """
        self.M.append(moments)

    def Build_Loads(self):
        """Method that sends the load_params out the the func specified in load_types; the load effects are calcuated at all_stations and fed back to the beam object's V[] and M[] lists. In order to combine all the V and M vectors, use Combine_Loads()
        """
        self.V.clear()
        self.M.clear()
        for j in range(0,len(self.load_types)):
            m = globals()[self.load_types[j]]
            func = getattr(m,'Get_Load_Effects')
            func(self, self.load_params[j])
        
    def Combine_Loads(self):
        """Index-wise combination of individual shear and moment vectors into a Total_V and Total_M vector.
        """
        self.Total_V = sum(self.V)
        self.Total_M = sum(self.M)
        
    def Get_Loads(self):
        """Function that runs the Build_Load() and Combine_Load() methods in order. This is the preferred method for UIs.
        """
        self.Build_Loads()
        self.Combine_Loads()
        
    def Show_Loads(self):
        """Returns a matplotlib plot of the loading in the load_params list.
        """
        x = self.all_stations
        y = self.Total_V
        zeros = np.zeros(self.all_stations.size)
        fig, ax = plt.subplots()
        ax.set_xlabel("Distance along Span")
        ax.set_ylabel("Loading")
        plt.plot(x,np.zeros(x.size),color="lightgray")
        loadArrows = list()
        for i in range(0,len(self.load_types)):
            mag = self.load_params[i][0]
            start = self.load_params[i][1]
            if self.load_types[i].split("_")[1] != "Point":
                end = self.load_params[i][2]
                dist_arrows = Load_Arrows.Draw(mag,start,end)
                for x in dist_arrows:
                    loadArrows.append(x)
            else:
                loadArrows.append(Load_Arrows.Draw(mag,start))
        for i in range(0,len(loadArrows)):
            ax.add_patch(loadArrows[i])
        imgPath = "{}\\{}-loading.png".format(save_folder, self.name)
        plt.savefig(imgPath,dpi=300,pad_inches=0.1)
        plt.show()
        
        
    def Show_Shear(self, show_each = True):
        """Returns a matplotlib plot of the shear diagram.

        Args:
            show_each (bool, optional): Shows the individual shears in the V[] list if true. Defaults to True.
        """
        x = self.all_stations
        y = self.Total_V
        zeros = np.zeros(self.all_stations.size)
        fig, ax = plt.subplots()
        ax.set_xlabel("Distance along Span")
        ax.set_ylabel("Shear")
        plt.plot(x,np.zeros(x.size),color="lightgray")
        if show_each == True:  
            for i in range(0,len(self.load_types)):
                if self.load_types[i] == "Simple_Point":
                    myColor = "pink"
                else:
                    myColor = "skyblue"
                plt.fill_between(x,zeros,self.V[i],color=myColor,alpha=0.25)
        plt.plot(x,y,color="dodgerblue",linewidth=2)
        imgPath = "{}\\{}-shear.png".format(save_folder, self.name)
        plt.savefig(imgPath,dpi=300,pad_inches=0.1)
        plt.show()

    def Show_Moment(self, show_each = True):
        """Returns a matplotlib plot of the moment diagram.

        Args:
            show_each (bool, optional): Shows the individual moments in the M[] list if true. Defaults to True.
        """
        x = self.all_stations
        y = self.Total_M
        zeros = np.zeros(self.all_stations.size)
        fig, ax = plt.subplots()
        ax.set_xlabel("Distance along Span")
        ax.set_ylabel("Moment")
        plt.plot(x,np.zeros(x.size),color="lightgray")
        if show_each == True:
            for i in range(0,len(self.load_types)):
                if self.load_types[i] == "Simple_Point":
                    myColor = "pink"
                else:
                    myColor = "skyblue"
                plt.fill_between(x,zeros,self.M[i],color=myColor,alpha=0.25)
        plt.plot(x,y,color="dodgerblue",linewidth=2)
        imgPath = "{}\\{}-moment.png".format(save_folder, self.name)
        plt.savefig(imgPath,dpi=300,pad_inches=0.1)
        plt.show()

    def Show_All(self, show_each = True):
        self.Show_Loads()
        self.Show_Shear(show_each)
        self.Show_Moment(show_each)

class Simple_Point:
    """Load class for concentrated load on a simply supported beam.
    """
    def Add_Load(beam: Beam, P, a):
        """Adds class name to beam load_type and loading parameters to beam load_params. 
        These are used later to populate the beam's V[] and M[] lists.

        Args:
            beam (Beam): Beam object to apply load information to
            P (float): Magnitude of the load. Note, positive P act downwards.
            a (float): Distance measured from the left support.
        """
        if (beam.boundaries[0] == 0) & (beam.boundaries[1] == 0):   # Prevents registering this load on a non-pin pin beam.      
            beam.Add_Stations([a])
            beam.Append_Load_Type('Simple_Point')
            beam.Append_Load_Params([P, a])

    def Get_Load_Effects(beam: Beam, load_params):
        """Calculates the shears and moments from the load_params stored in the beam object, then sends the values back to the beam object's V[] and M[] lists.

        Args:
            beam (Beam): Beam object to operate on
            load_params ([]): list of load_params stored in the beam load_params[] list. For this class, it will be [P,a] per Add_Load.
        """
        P = load_params[0]
        a = load_params[1]
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
    """Load class for uniform distributed load (UDL) on a simply supported beam.
    """
    def Add_Load(beam: Beam, w, a, b):
        """Adds class name to beam load_type and loading parameters to beam load_params. 
        These are used later to populate the beam's V[] and M[] lists. 
        Note that if c = length - b - a. If a = c = 0 or b = length, UDL over entire beam. Otherwise partially loaded.

        Args:
            beam (Beam): Beam object to apply load information to
            w (float): Intensity of the load. Note, positive w act downwards.
            a (float): Distance of unloaded beam segment, measured from the left support.
            b (float): Distance of the loaded segment of the beam.
        """
        if (beam.boundaries[0] == 0) & (beam.boundaries[1] == 0):   # Prevents registering this load on a non-pin pin beam.      
            beam.Add_Stations([a,b])
            beam.Append_Load_Type('Simple_UDL')
            beam.Append_Load_Params([w, a, b])

    def Get_Load_Effects(beam: Beam, load_params):
        """Calculates the shears and moments from the load_params stored in the beam object, then sends the values back to the beam object's V[] and M[] lists.

        Args:
            beam (Beam): Beam object to operate on
            load_params ([]): list of load_params stored in the beam load_params[] list. For this class, it will be [w,a,b] per Add_Load.
        """
        w = load_params[0]
        a = load_params[1]
        b = load_params[2]
        locs = np.copy(beam.all_stations)
        shears = np.zeros(locs.size)
        moments = np.zeros(locs.size)
        length = locs[np.argmax(locs)]
        R1 = w * b / (2 * length) * (2 * (length - a - b) + b)
        R2 = w * b / (2 * length)  * (2 * a + b)
        for i in range(0,np.argmax(locs) + 1):
            if locs[i] <= a:
                shears[i] = R1
                moments[i] = R1 * locs[i]
            elif (locs[i]>a) & (locs[i]<a+b):
                shears[i] = R1 - w * (locs[i]-a)
                moments[i] = R1 * locs[i] - w/2 * (locs[i]-a)**2
            else:
                shears[i] = R1 - b * w
                moments[i] = R2 * (length - locs[i])
        beam.Append_Shears(shears)
        beam.Append_Moments(moments)


class Load_Arrows:
    """Static class that creates load arrow(s) for a MatPlotLib plot.
    """
    def Draw(mag, start, end = 0, qty = 4):
        if end == 0:
            return mpatches.FancyArrowPatch((start, mag), (start, 0),mutation_scale=20,arrowstyle="->",shrinkA=0,shrinkB=0)
        else:
            arrows = list()
            xs = np.linspace(start,end,qty)
            ys = np.zeros(qty)
            ys.fill(mag)
            for i in range(0, qty):
                arrows.append(mpatches.FancyArrowPatch((xs[i], ys[i]), (xs[i], 0),mutation_scale=20,arrowstyle="->",shrinkA=0,shrinkB=0))
            xyA = (start,mag)
            xyB = (end,mag)
            connector = mpatches.ConnectionPatch(xyA,xyB,coordsA="data",coordsB="data",arrowstyle="-",shrinkA=0,shrinkB=0,mutation_scale=20)
            arrows.append(connector)
            return arrows
            
