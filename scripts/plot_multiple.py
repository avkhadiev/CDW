#!/usr/bin/env python
import os
import matplotlib.pyplot as plt

dirname   = 'output'
plotname = 'ET_eta'

# make a list of all txt files in the directory
text_files = [f for f in os.listdir( dirname ) if f.endswith('field.txt')]

# open a file;
i = 0;
velocity  = [[], [], [], []]
field     = [[], [], [], []]
new_label = []
print("making a plot %s" % ( plotname ))
fig = plt.figure()
plt.xlabel("Electric Field")
plt.ylabel("Average Phase Velocity")
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.000001, 1)
plt.xlim(0.01, 1)
fig.suptitle( "Threshold Field ( Thermal Noise )" , fontsize=20)
ax1 = fig.add_subplot(111);
for filename in text_files:
    velocity.append( [] );
    field.append( [] );
    new_label.append( filename.split('.')[0] )
    print(i);
    with open("%s/%s" % (dirname, filename), 'r') as f:
        for line in f:
            cols = [float(x) for x in line.split()]
            velocity[i].append( cols[0] );
            field[i].append( cols[1] );
    i = i + 1;
ax1.scatter( field[0], velocity[0], c = 'b', label = new_label[0] )
ax1.scatter( field[1], velocity[1], c = 'g', label = new_label[1] )
ax1.scatter( field[2], velocity[2], c = 'm', label = new_label[2] )
plt.legend(loc = 'upper left')
fig.savefig( "%s/%s/%s.%s" % ( 'output', "plots", plotname, "png"))
plt.close()
