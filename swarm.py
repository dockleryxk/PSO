#Particle Swarm Optimizatiion
#Richard Jeffords
#05/2015

import sys 
import math
import random
from random import randint


if sys.argv[1] == '-v':
    verbose = True
else:
    verbose = False

if sys.argv[1] == 'test' or sys.argv[2] == 'test':
    thisIsATest = True
elif len(sys.argv) < 12:
    sys.stderr.write('USAGE: [-v] swarm.py numParticles inertia cognition socialRate localRate worldWidth worldHeight maxVelocity maxEpochs k fname\n')
    sys.exit(0)
else:
    thisIsATest = False

if(not thisIsATest):
    ########### command line args
    numParticles = int(sys.argv[1])
    inertia      = float(sys.argv[2])
    cognition    = float(sys.argv[3])
    socialRate   = float(sys.argv[4])
    localRate    = float(sys.argv[5])
    worldWidth   = float(sys.argv[6])
    worldHeight  = float(sys.argv[7])
    maxVelocity  = float(sys.argv[8])
    maxEpochs    = int(sys.argv[9])
    k            = int(sys.argv[10])
    fname        = str(sys.argv[11])
    fname += '.csv'

else:
    ########### for testing
    numParticles = 20
    inertia      = 0.95
    cognition    = 2.0
    socialRate   = 2.0
    localRate    = 2.0
    worldWidth   = 100.0
    worldHeight  = 100.0
    maxVelocity  = 2.0
    k            = 0
    maxEpochs = 10000

########### data representation
pList = []
class Particle:
    #value, x_pos, y_pos
    gBest     = [0.0, 0, 0]
    bestIndex = 0

    #takes index in pList as constructor argument
    def __init__(self, i):
        #x,y coords, randomly initialized
        self.x          = randint(-worldWidth/2,worldWidth/2)
        self.y          = randint(-worldHeight/2,worldHeight/2)
        #x,y velocity
        self.velocity_x = 0.0
        self.velocity_y = 0.0
        #personal best
        #[fitness value, x coord, y coord]
        self.pBest      = [Q(self.x, self.y), self.x, self.y]
        self.index      = i
        #local best
        self.lBest      = []
        self.lBestIndex = 0
        #array for neighbor indicies
        self.neighbors  = []
    #for printing particle info
    def __str__(self):
        if k > 0:
            return '  i: '+str(self.index)+'\n  x: '+str(self.x)+'\n  y: '+str(self.y)+'\nv_x: '+str(self.velocity_x)+'\nv_y: '+str(self.velocity_y)+'\n  b: '+str(self.pBest[0])+'\n  l: '+str(self.lBest)+'\n'
        else:
            return '  i: '+str(self.index)+'\n  x: '+str(self.x)+'\n  y: '+str(self.y)+'\nv_x: '+str(self.velocity_x)+'\nv_y: '+str(self.velocity_y)+'\n  b: '+str(self.pBest[0])+'\n'

###########

def createParticles():
    global pList
    global numParticles
    global k
    #create particle list
    for i in range(0,numParticles):
        pList.append(Particle(i))

    #fill neighbor lists
    if k > 0:
        for p in pList:
            for x in range(p.index-(k/2),p.index+(k/2)+1):
                if x > numParticles:
                    p.neighbors.append(x%numParticles)
                elif x < 0:
                    p.neighbors.append(numParticles+x)
                elif x == numParticles:
                    p.neighbors.append(0)
                else:
                    p.neighbors.append(x)
        updatelBest()
    
    #initialize global and local bests
    updategBest()


########### distance functions

def mdist():
    global worldWidth
    global worldHeight
    return float(math.sqrt((worldWidth**2.0)+(worldHeight**2.0))/2.0)
########### 
def pdist(p_x,p_y):
    return float(math.sqrt(((p_x-20.0)**2.0) + ((p_y-7.0)**2.0)))
########### 
def ndist(p_x,p_y):
    return float(math.sqrt(((p_x+20.0)**2.0) + ((p_y+7.0)**2.0)))

########### Problem 1
def Q(p_x,p_y):
    return float(100.0 * (1.0 - (pdist(p_x,p_y)/mdist())))

########### Problem 2
'''
def Q(p_x,p_y):
    return float((9.0 * max(0.0, 10.0 - (pdist(p_x,p_y)**2))) + (10.0 * (1.0 - (pdist(p_x,p_y)/mdist()))) + (70.0 * (1.0 - (ndist(p_x,p_y)/mdist()))))
'''

###########

def printParticles():
    global pList
    print '\ngBest: ', Particle.gBest
    print 'index: ', Particle.bestIndex, '\n'
    for p in pList:
        print p

###########

def updateVelocity():
    global pList
    global inertia
    global socialRate
    global cognition
    global maxVelocity
    global localRate
    global k
    rand1 = random.uniform(0.0,1.0)
    rand2 = random.uniform(0.0,1.0)
    rand3 = random.uniform(0.0,1.0)
    v_x    = 0.0
    v_y    = 0.0
    v_x2   = 0.0
    v_y2   = 0.0
    flag_x = False
    flag_y = False

    for p in pList:
        #velocity update with neighbors
        if k > 0:
            v_x = inertia * p.velocity_x + cognition * rand1 * (p.pBest[1] - p.x) + socialRate * rand2 * (Particle.gBest[1] - p.x) + localRate * rand3 * (p.lBest[1] - p.x)
            v_y = inertia * p.velocity_y + cognition * rand1 * (p.pBest[2] - p.y) + socialRate * rand2 * (Particle.gBest[2] - p.y) + localRate * rand3 * (p.lBest[2] - p.y)

        #velocity update without neighbors
        #velocity' = inertia * velocity + c_1 * r_1 * (personal_best_position - position) + c_2 * r_2 * (global_best_position - position) 
        else: 
            v_x = inertia * p.velocity_x + cognition * rand1 * (p.pBest[1] - p.x) + socialRate * rand2 * (Particle.gBest[1] - p.x)
            v_y = inertia * p.velocity_y + cognition * rand1 * (p.pBest[2] - p.y) + socialRate * rand2 * (Particle.gBest[2] - p.y)

        #scale velocity
        #if abs(velocity) > maximum_velocity^2 
        #velocity = (maximum_velocity/sqrt(velocity_x^2 + velocity_y^2)) * velocity 
        if abs(v_x) > maxVelocity:
            v_x2 = (maxVelocity/math.sqrt(v_x**2 + v_y**2)) * v_x
            flag_x = True
        if abs(v_y) > maxVelocity:
            v_y2 = (maxVelocity/math.sqrt(v_y**2 + v_y**2)) * v_y
            flag_y = True

        #use flag to determine which temp variable to use
        #that way, v_x and v_y aren't altered by the scaling
        if flag_x:
            p.velocity_x = v_x2
        else:
            p.velocity_x = v_x
        if flag_y:
            p.velocity_y = v_y2
        else:
            p.velocity_y = v_y

###########

def updatePosition():
    global pList
    for p in pList:
        #position' = position + velocity' 
        p.x += p.velocity_x
        p.y += p.velocity_y

###########

def updatepBest():
    global pList
    for p in pList:
        #if(Q(position) > Q(personal_best_position)) 
        #personal_best_position = position 
        if Q(p.x,p.y) > p.pBest[0]:
            p.pBest = [Q(p.x,p.y), p.x, p.y]

###########

def updategBest():
    global pList
    tmp = Particle.gBest
    tmpIndex = Particle.bestIndex
    for p in pList:
        #if(Q(position) > Q(global_best_position)) 
        #global_best_position = position 
        if Q(p.x,p.y) > tmp[0]:
            tmp      = [Q(p.x,p.y), p.x, p.y]
            tmpIndex = p.index
    Particle.gBest     = tmp
    Particle.bestIndex = tmpIndex

###########

def updatelBest():
    global pList
    tmp = [0.0, 0, 0]
    tmpIndex = 0
    for p in pList:
        for n in p.neighbors:
            #find the local best Q value
            if Q(pList[n].x,pList[n].y) > tmp[0]:
                tmp = [Q(pList[n].x, pList[n].y), pList[n].x, pList[n].y]
                tmpIndex = pList[n].index
        p.lBest      = tmp
        p.lBestIndex = tmpIndex
        #reset tmp
        tmp = [0.0, 0, 0]

###########

def calcError():
    global pList
    global numParticles
    error_x = 0.0
    error_y = 0.0

    #for each particle k:
    #error_x += (position_x[k] - global_best_position_x)^2 
    #error_y += (position_y[k] - global_best_position_y)^2
    for p in pList:
        error_x += (p.x - Particle.gBest[1])**2
        error_y += (p.y - Particle.gBest[2])**2

    #Then
    #error_x = sqrt((1/(2*num_particles))*error_x) 
    #error_y = sqrt((1/(2*num_particles))*error_y) 
    error_x = math.sqrt((1.0/(2.0*numParticles))*error_x)
    error_y = math.sqrt((1.0/(2.0*numParticles))*error_y)

    return [error_x, error_y]

###########

def paramsToCSV():
    global numParticles 
    global inertia     
    global cognition    
    global socialRate   
    global localRate    
    global worldWidth   
    global worldHeight  
    global maxVelocity  
    global maxEpochs
    global k            
    global fname
    f = open(fname, 'a+')
    f.write('numParticles,inertia,cognition,socialRate,localRate,worldWidth,worldHeight,maxVelocity,maxEpochs,k\n'+str(numParticles)+','+str(inertia)+','+str(cognition)+','+str(socialRate)+','+str(localRate)+','+str(worldWidth)+','+str(worldHeight)+','+str(maxVelocity)+','+str(maxEpochs)+','+str(k)+'\nx error,y error\n')
    f.close()


###########

def errorToCSV(e):
    global fname
    f = open(fname, 'a+')
    f.write(str(e[0])+','+str(e[1])+'\n')
    f.close()
        
###########

def plotToCSV():
    global pList
    f = open(fname,'a+')
    f.write('\n\n\nx values,y values')
    for p in pList:
        f.write(str(p.x)+','+str(p.y)+'\n')
    f.close()

########### main

createParticles()
if not thisIsATest:
    paramsToCSV()

epochs = 0
#######
while True:
    ###
    updateVelocity()
    ###
    updatePosition()
    ###
    updatepBest()
    ###
    updategBest()
    if k > 0:
        updatelBest()
    ###
    error = calcError()
    if not thisIsATest:
        errorToCSV(error)
    if verbose:
        print error
    ###
    if error[0] < 0.01 and error[1] < 0.01:
        break
    elif epochs > maxEpochs:
        break
    ###
    epochs += 1
#######
updatepBest()
updategBest()
if k > 0:
    updatelBest()
printParticles()    
if not thisIsATest:
    plotToCSV()
print 'epochs: ', epochs
