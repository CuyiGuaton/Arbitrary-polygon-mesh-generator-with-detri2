import random
import getopt, sys
import numpy as np
from math import sqrt

def write_points(xPoints, yPoints):
    largo =len(xPoints)
    f = open('data.dat', 'w')
    f.write("2 rbox {} D2 z\n".format(largo))
    f.write('{}\n'.format(largo))
    for i in range(0,largo):
        f.write('{0} {1}\n'.format(xPoints[i], yPoints[i]))

    f.close()



def equitalero(nPoints, xPoints, yPoints):
    L = 2
    h = sqrt(3) # = sqrt(3)*2/2
    i = 0
    for i in range(0, nPoints):
        for j in range(0, nPoints):
            xPoints.append(L*j if i%2 == 0 else L*j + L/2)
            yPoints.append(i*h)
        if i%2 == 0:
            xPoints.append(L*nPoints + L/2)
            yPoints.append(i*h)
        else:
            xPoints.append(L*0)
            yPoints.append(i*h)
    

def generate_random_points(nPoints, xPoints, yPoints):

    uwu = nPoints

    for _ in range(0,nPoints-4):
        x = random.randrange(0, uwu)
        y = random.randrange(0, uwu)
        xPoints.append(int(x))
        yPoints.append(int(y))

    xPoints.append(0)
    yPoints.append(0)

    xPoints.append(uwu)
    yPoints.append(uwu)

    xPoints.append(0)
    yPoints.append(uwu)

    xPoints.append(uwu)
    yPoints.append(0)
    move_edges(nPoints, xPoints, yPoints)


def move_edges(nPoints, xPoints, yPoints):
    tolerance =  0.05
    n = nPoints

    for i in range(0, len(xPoints)):
        #eje x
        if xPoints[i] >= nPoints*(1.0-tolerance): 
            xPoints[i] = n
            continue
        if yPoints[i] >= nPoints*(1.0-tolerance): 
            yPoints[i] = n
            continue
        if xPoints[i] <= nPoints*tolerance: 
            xPoints[i] = 0
            continue
        if yPoints[i] <= nPoints*tolerance: 
            yPoints[i] = 0
            continue

    
def generate_uniform_points(nPoints, xPoints, yPoints, dist = 0):
    
    for i in range(int(np.sqrt(nPoints))):
        for j in range(int(np.sqrt(nPoints))):
            xPoints.append(j)
            yPoints.append(i)   

def generate_semi_uniform_points(nPoints, xPoints, yPoints, dist = 0):
    uwu = int(np.sqrt(nPoints))
    for i in range(uwu):
        for j in range(uwu):
            xPoints.append(j + random.uniform(-dist, dist))
            yPoints.append(i + random.uniform(-dist, dist))    
    #xPoints.append(0 - dist)
    #yPoints.append(0 - dist)
#
    #xPoints.append(uwu + dist)
    #yPoints.append(uwu + dist)
#
    #xPoints.append(0 - dist)
    #yPoints.append(uwu + dist)
#
    #xPoints.append(uwu + dist)
    #yPoints.append(0 - dist)


def main():
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]

    if(len(argument_list) == 2):
        nPoints = int(argument_list[0])
        selection = int(argument_list[1])
        dist = 0
    else:
        nPoints = int(argument_list[0])
        selection = int(argument_list[1])
        dist = float(argument_list[2])

    random.seed(138)
    xPoints = []
    yPoints = []

    #dist = nPoints*0.7
    
    if(selection == 0):
        generate_random_points(nPoints, xPoints, yPoints)
    elif selection == 1:
        equitalero(nPoints, xPoints, yPoints)
    elif selection == 2:
        generate_semi_uniform_points(nPoints, xPoints, yPoints, dist)

    write_points(xPoints, yPoints)
    

if __name__ == "__main__":
    main()
