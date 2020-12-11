import random
import getopt, sys
import numpy as np
from math import sqrt

def write_points(nPoints, xPoints, yPoints):
    largo =len(xPoints)
    f = open('autodata.node', 'w')
    f.write("{} 2 0 0\n".format(largo))
    for i in range(0,largo):
        f.write('{0} {1} {2}\n'.format(i, xPoints[i], yPoints[i]))
    f.write('\n')
    f.close()



def equitalero(nPoints, xPoints, yPoints):
    L = 2
    h = sqrt(3) # = sqrt(3)*2/2
    i = 0
    for i in range(0, nPoints):
        if i%2 != 0:
            xPoints.append(L*0)
            yPoints.append(i*h)
        for j in range(0, nPoints):
            xPoints.append(L*j if i%2 == 0 else L*j + L/2)
            yPoints.append(i*h)
        if i%2 == 0:
            xPoints.append(L*(nPoints-1) + L/2)
            yPoints.append(i*h)
        
 

def generate_random_points(nPoints, xPoints, yPoints):
    uwu = nPoints
    s = set([(0, 0), (uwu,uwu), (0, uwu), (uwu,0)])

    while len(s) != uwu:
        x = random.randrange(0, uwu)
        y = random.randrange(0, uwu)
        s.add(move_point(uwu, x, y))

    for _ in range(0,len(s)):
        l = s.pop()
        xPoints.append(l[0])
        yPoints.append(l[1])



def move_point(nPoints, xPoints, yPoints):
    tolerance =  0.005
    #tolerance =  5/(nPoints/10)
    n = nPoints
    r = random.uniform(0, 1)
    if r > 0.5:
        if xPoints >= nPoints*(1.0-tolerance): 
            xPoints = n
        if yPoints >= nPoints*(1.0-tolerance): 
            yPoints = n
        if xPoints <= nPoints*tolerance: 
            xPoints = 0
        if yPoints <= nPoints*tolerance: 
            yPoints = 0
    else:
        if xPoints <= nPoints*tolerance: 
            xPoints = 0            
        if yPoints <= nPoints*tolerance: 
            yPoints = 0
        if xPoints >= nPoints*(1.0-tolerance): 
            xPoints = n
        if yPoints >= nPoints*(1.0-tolerance): 
            yPoints = n
    return (xPoints, yPoints)

    
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


 
 #Input Largo inicial, incremento
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

    write_points(nPoints, xPoints, yPoints)
    

if __name__ == "__main__":
    main()
