#!/usr/bin/env python
import os
import matplotlib.pyplot as plt

dirname = 'output'

# make a list of all txt files in the directory
text_files = [f for f in os.listdir( dirname ) if f.endswith('.txt')]

# open a file;
for filename in text_files:
    with open("%s/%s" % (dirname, filename), 'r') as f:
        print("opened file %s\n" % (filename))
    #   make an array of all numbers in the first column
    #      and all numbers in the second column
        variable = []
        time     = []
        plotname = filename.split('.')[0]
        for line in f:
            cols = [float(x) for x in line.split()]
            variable.append( cols[0] );
            time.append( cols[1] );
        print("making a plot %s" % ( plotname ))
        fig = plt.figure(figsize=(20,2))
        plt.plot( time, variable )
        fig.suptitle( plotname , fontsize=10)
        plt.xlabel("Time")
        plt.ylabel("%s" % (filename.split('_')[0]))
        fig.savefig( "%s/%s/%s.%s" % ( dirname, "plots", plotname, "png"))
        plt.close()


#   make a plot with elements from second column on the x-axis
#      and elements from the first column on the y-axis
#   label the plot with the name of the file;
#   label the x axis "time".
#   if the name of the file contains "phase", label the y-axis "phase"
#   if the name of the file contains "momentum", label the y-axis "momentum"
#   save the plot as a png file






