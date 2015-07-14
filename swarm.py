#Particle Swarm Optimizatiion
#Copyright Richard Jeffords 2015
#GNU General Public License
"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
     
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
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
    """assign arguments based on command line arguments"""
    ########### command line args
    NP        = int(sys.argv[1])
    I         = float(sys.argv[2])
    C         = float(sys.argv[3])
    SR        = float(sys.argv[4])
    LR        = float(sys.argv[5])
    WW        = float(sys.argv[6])
    WH        = float(sys.argv[7])
    MV        = float(sys.argv[8])
    maxEpochs = int(sys.argv[9])
    K         = int(sys.argv[10])
    FN        = str(sys.argv[11])
    FN        += '.csv'

else:
    """use sample arguments to test the script"""
    ########### for testing
    NP        = 20
    I         = 0.95
    C         = 2.0
    SR        = 2.0
    LR        = 2.0
    WW        = 100.0
    WH        = 100.0
    MV        = 2.0
    K         = 0
    maxEpochs = 10000
    FN        = ''

########### data representation

class ParticleList:
    """ParticleList encapsulates the list of particles and functions used to 
    manipulate their attributes
    """

    def __init__(self, NP, I, C, SR, LR, WW, WH, MV, K, FN):
        """create an array, assign values, and initialize each particle"""
        self.pList        = []
        self.numParticles = NP
        self.inertia      = I
        self.cognition    = C
        self.socialRate   = SR
        self.localRate    = LR
        self.worldWidth   = WW
        self.worldHeight  = WH
        self.maxVelocity  = MV
        self.k            = K
        self.fname        = FN
        self.createParticles()

    def createParticles(self):
        """create a list of particles and then create neighborhoods if it's called for (k > 0)"""
        for i in range(0,self.numParticles):
            self.pList.append(self.Particle(i, self.worldWidth, self.worldHeight, self.k))

        #fill neighbor lists
        if self.k > 0:
            for p in self.pList:
                for x in range(p.index-(self.k/2),p.index+(self.k/2)+1):
                    if x > self.numParticles:
                        p.neighbors.append(x%self.numParticles)
                    elif x < 0:
                        p.neighbors.append(self.numParticles+x)
                    elif x == self.numParticles:
                        p.neighbors.append(0)
                    else:
                        p.neighbors.append(x)
            self.updatelBest()
        
        #initialize global and local bests
        self.updategBest()

    ###########

    class Particle:
        """this class is used for each particle in the list and all of their attributes"""
        #[Q value, x_pos, y_pos]
        gBest     = [0.0, 0, 0]
        bestIndex = 0

        #takes index in pList as constructor argument
        def __init__(self, i, worldWidth, worldHeight, K):
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
            self.k          = K
        #for printing particle info
        def __str__(self):
            """Creates string representation of particle"""
            ret = """  i: {self.index!s}
            x: {self.x!s}
            y: {self.y!s}
            v_x: {self.velocity_x!s}
            v_y: {self.velocity_y!s}
            b: {self.pBest[0]!s}""".format(**locals())
            if self.k > 0:
                return ret+'  l: '+str(self.lBest)+'\n'
            else:
                return ret+'\n'

    ###########

    def printParticles(self):
        """prints out useful info about each particle in the list"""
        print '\ngBest: ', self.Particle.gBest
        print 'index: ', self.Particle.bestIndex, '\n'
        for p in self.pList:
            print p

    ###########

    def updateVelocity(self):
        """at each timestep or epoch, the velocity of each particle is updated
        based on the inertia, current velocity, cognition, social rate, 
        and optionally local rate. Of course, there's some choas too.
        """
        rand1 = random.uniform(0.0,1.0)
        rand2 = random.uniform(0.0,1.0)
        rand3 = random.uniform(0.0,1.0)
        v_x    = 0.0
        v_y    = 0.0
        v_x2   = 0.0
        v_y2   = 0.0
        flag_x = False
        flag_y = False

        for p in self.pList:
            #velocity update with neighbors
            if self.k > 0:
                v_x = inertia * p.velocity_x + cognition * rand1 * (p.pBest[1] - p.x) + socialRate * rand2 * (Particle.gBest[1] - p.x) + localRate * rand3 * (p.lBest[1] - p.x)
                v_y = inertia * p.velocity_y + cognition * rand1 * (p.pBest[2] - p.y) + socialRate * rand2 * (Particle.gBest[2] - p.y) + localRate * rand3 * (p.lBest[2] - p.y)

            #velocity update without neighbors
            #velocity' = inertia * velocity + c_1 * r_1 * (personal_best_position - position) + c_2 * r_2 * (global_best_position - position) 
            else: 
                v_x = self.inertia * p.velocity_x + self.cognition * rand1 * (p.pBest[1] - p.x) + self.socialRate * rand2 * (self.Particle.gBest[1] - p.x)
                v_y = self.inertia * p.velocity_y + self.cognition * rand1 * (p.pBest[2] - p.y) + self.socialRate * rand2 * (self.Particle.gBest[2] - p.y)

            #scale velocity
            #if abs(velocity) > maximum_velocity^2 
            #velocity = (maximum_velocity/sqrt(velocity_x^2 + velocity_y^2)) * velocity 
            if abs(v_x) > self.maxVelocity:
                v_x2 = (self.maxVelocity/math.sqrt(v_x**2 + v_y**2)) * v_x
                flag_x = True
            if abs(v_y) > self.maxVelocity:
                v_y2 = (self.maxVelocity/math.sqrt(v_y**2 + v_y**2)) * v_y
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

    def updatePosition(self):
        """update particle postions based on velocity"""
        for p in self.pList:
            #position' = position + velocity' 
            p.x += p.velocity_x
            p.y += p.velocity_y

    ###########

    def updatepBest(self):
        """at each epoch, check to see if each particle's current position
        is its best (or closest to the solution) yet
        """
        for p in self.pList:
            #if(Q(position) > Q(personal_best_position)) 
            #personal_best_position = position 
            if Q(p.x,p.y) > p.pBest[0]:
                p.pBest = [Q(p.x,p.y), p.x, p.y]

    ###########

    def updategBest(self):
        """find the best position of all the particles in the list"""
        tmp = self.Particle.gBest
        tmpIndex = self.Particle.bestIndex
        for p in self.pList:
            #if(Q(position) > Q(global_best_position)) 
            #global_best_position = position 
            if Q(p.x,p.y) > tmp[0]:
                tmp      = [Q(p.x,p.y), p.x, p.y]
                tmpIndex = p.index
        self.Particle.gBest     = tmp
        self.Particle.bestIndex = tmpIndex

    ###########

    def updatelBest(self):
        """optionally find the best position out of a neighborhood"""
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

    def calcError(self):
        """calculate the error at each epoch"""
        error_x = 0.0
        error_y = 0.0

        #for each particle p:
        #error_x += (position_x[k] - global_best_position_x)^2 
        #error_y += (position_y[k] - global_best_position_y)^2
        for p in self.pList:
            error_x += (p.x - self.Particle.gBest[1])**2
            error_y += (p.y - self.Particle.gBest[2])**2

        #Then
        #error_x = sqrt((1/(2*num_particles))*error_x) 
        #error_y = sqrt((1/(2*num_particles))*error_y) 
        error_x = math.sqrt((1.0/(2.0*self.numParticles))*error_x)
        error_y = math.sqrt((1.0/(2.0*self.numParticles))*error_y)

        return [error_x, error_y]

    ###########

    def paramsToCSV(self):
        """put the parameters at the top of the CSV file"""
        f = open(fname, 'a+')
        f.write('numParticles,inertia,cognition,socialRate,localRate,worldWidth,worldHeight,maxVelocity,maxEpochs,k\n'+str(self.numParticles)+','+str(self.inertia)+','+str(self.cognition)+','+str(self.socialRate)+','+str(self.localRate)+','+str(self.worldWidth)+','+str(self.worldHeight)+','+str(self.maxVelocity)+','+str(self.maxEpochs)+','+str(self.k)+'\nx error,y error\n')
        f.close()

    ###########

    def errorToCSV(self, e):
        """print the error at each epoch to produce an error over time graph"""
        f = open(fname, 'a+')
        f.write(str(e[0])+','+str(e[1])+'\n')
        f.close()
            
    ###########

    def plotToCSV(self):
        """print the points at the end to create a scatter plot, or at each epoch
        to try for a gif animation
        """
        f = open(fname,'a+')
        f.write('\n\n\nx values,y values')
        for p in self.pList:
            f.write(str(p.x)+','+str(p.y)+'\n')
        f.close()

########### distance functions

def mdist():
    global WW
    global WH
    return float(math.sqrt((WW**2.0)+(WH**2.0))/2.0)
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
#initialize particle list

particles = ParticleList(NP, I, C, SR, LR, WW, WH, MV, K, FN)

########### main

if not thisIsATest:
    particles.paramsToCSV()

epochs = 0
#######
while True:
    """each run through this loop represents and epoch"""
    ###
    particles.updateVelocity()
    ###
    particles.updatePosition()
    ###
    particles.updatepBest()
    ###
    particles.updategBest()
    if K > 0:
        particles.updatelBest()
    ###
    error = particles.calcError()
    if not thisIsATest:
        particles.errorToCSV(error)
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
particles.updatepBest()
particles.updategBest()
if K > 0:
    particles.updatelBest()
particles.printParticles()    
if not thisIsATest:
    particles.plotToCSV()
print 'epochs: ', epochs


