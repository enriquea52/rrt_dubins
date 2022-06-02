#!/usr/bin/env python3


import rrt_star_planner as rrt_s_d

import numpy as np
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt

if __name__ == "__main__":

    planner = rrt_s_d.rrt_star_dubins(-7.5, -7.5, 6.0, 0.0, 20.0, 10, 10, 3000)

    #planner = sbp.rrt_star_dubins(6.0, 0.0, -7.5, -7.5, 5.0, 10, 10, 3000)
    #dubins = planner.test_dubins_trayectory(10)

    
    points = np.asarray(planner.plan()); #plt.scatter(points[:,0], points[:,1], s=0.5, c = "red");print("points shape: ",points.shape)


    
    #tree = np.asarray(planner.tree()); plt.scatter(tree[:,0], tree[:,1], s=0.5, c = "blue"); print("tree shape: ",tree.shape)
    path = np.asarray(planner.path()); #plt.scatter(path[:,0], path[:,1], s=0.5, c = "black"); print("path shape: ",path.shape)
    
    indexes = np.arange(0, path.shape[0], 20)

    path_sub = path[indexes]

    plt.scatter(path_sub[:,0], path_sub[:,1], s=0.5, c = "red")
    
    
    
    #path_points = np.asarray(planner.path_points()); plt.scatter(path_points[:,0], path_points[:,1], s=2, c = "blue"); print("path_points shape: ",path_points.shape)
    #rect = Rectangle((3,-8),2,16,linewidth=1,edgecolor='r',facecolor='none')
    #plt.gca().add_patch(rect)
    
    
    


    plt.xlim((-10, 10))
    plt.ylim((-10, 10))
    plt.show()
