#Particle Swarm Optimizatiion
#Copyright Richard Jeffords 2015
__author = "Richard Jeffords"
__email  = "rajeffords@gmail.com"
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
    offset  = 1
else:
    verbose = False
    offset  = 0

if sys.argv[1] == 'test' or sys.argv[2] == 'test':
    a_test = True
elif len(sys.argv) < 12:
    sys.stderr.write('USAGE: [-v] swarm.py num_particles inertia cognition social_rate local_rate world_width world_height max_velocity max_epochs num_neighbors fname\n')
    sys.exit(0)
else:
    a_test = False

if(not a_test):
    """assign arguments based on command line arguments"""
    ########### command line args
    params = ({'num_particles' : int(sys.argv[1 + offset]),
               'inertia'       : float(sys.argv[2 + offset]),
               'cognition'     : float(sys.argv[3 + offset]),
               'social_rate'   : float(sys.argv[4 + offset]),
               'local_rate'    : float(sys.argv[5 + offset]),
               'world_width'   : float(sys.argv[6 + offset]),
               'world_height'  : float(sys.argv[7 + offset]),
               'max_velocity'  : float(sys.argv[8 + offset]),
               'max_epochs'    : int(sys.argv[9 + offset]),
               'num_neighbors' : int(sys.argv[10 + offset]),
               'fname'         : (str(sys.argv[11 + offset]) + '.csv')})

else:
    """use sample arguments to test the script"""
    ########### for testing
    params = ({'num_particles' : 20,
               'inertia'       : 0.95,
               'cognition'     : 2.0,
               'social_rate'   : 2.0,
               'local_rate'    : 2.0,
               'world_width'   : 100.0,
               'world_height'  : 100.0,
               'max_velocity'  : 2.0,
               'max_epochs'    : 10000,
               'num_neighbors' : 0,
               'fname'         : ''})

########### data representation

class Particle_List:
    """Particle_List encapsulates the list of particles and functions used to 
    manipulate their attributes
    """

    def __init__(self, args):
        """create an array, assign values, and initialize each particle"""
        self.p_list        = []
        self.num_particles = args['num_particles']
        self.inertia       = args['inertia']
        self.cognition     = args['cognition']
        self.social_rate   = args['social_rate']
        self.local_rate    = args['local_rate']
        self.world_width   = args['world_width']
        self.world_height  = args['world_height']
        self.max_velocity  = args['max_velocity']
        self.max_epochs    = args['max_epochs']
        self.num_neighbors = args['num_neighbors']
        self.fname         = args['fname']
        self._create_particles()

    def _create_particles(self):
        """create a list of particles and then create neighborhoods if it's called for (k > 0)"""
        for i in range(self.num_particles):
            self.p_list.append(self.Particle(i, self.world_width, self.world_height, self.num_neighbors))

        #fill neighbor lists
        if self.num_neighbors > 0:
            for p in self.p_list:
                begin = p.index - (self.num_neighbors/2)
                end   = p.index + (self.num_neighbors/2) + 1

                for x in range(begin,end):
                    if x > self.num_particles:
                        p.neighbors.append(x%self.num_particles)
                    elif x < 0:
                        p.neighbors.append(self.num_particles+x)
                    elif x == self.num_particles:
                        p.neighbors.append(0)
                    else:
                        p.neighbors.append(x)
            self.update_local_best()
        
        #initialize global and local bests
        self.update_global_best()

    ###########

    class Particle:
        """this class is used for each particle in the list and all of their attributes"""
        #[Q value, x_pos, y_pos]
        global_best = [0.0, 0, 0]
        best_index  = 0

        #takes index in p_list as constructor argument
        def __init__(self, i, world_width, world_height, num_neighbors):
            #x,y coords, randomly initialized
            self.x                = randint(-world_width/2,world_width/2)
            self.y                = randint(-world_height/2,world_height/2)
            #x,y velocity
            self.velocity_x       = 0.0
            self.velocity_y       = 0.0
            #personal best
            #[fitness value, x coord, y coord]
            self.personal_best    = [Q(self.x, self.y), self.x, self.y]
            self.index            = i
            #local best
            self.local_best       = []
            self.local_best_index = 0
            #array for neighbor indicies
            self.neighbors        = []
            self.num_neighbors    = num_neighbors
        #for printing particle info
        def __str__(self):
            """Creates string representation of particle"""
            ret = """
            index: {self.index!s}
            x coordinate: {self.x!s}
            y coordinate: {self.y!s}
            x velocity: {self.velocity_x!s}
            y velocity: {self.velocity_y!s}
            personal best: {self.personal_best[0]!s}
            """.format(**locals())
            if self.num_neighbors > 0:
                return ret+'  local best: '+str(self.local_best)+'\n'
            else:
                return ret+'\n'

    ###########

    def print_particles(self):
        """prints out useful info about each particle in the list"""
        print '\nglobal_best: ', self.Particle.global_best
        print 'index: ', self.Particle.best_index, '\n'
        for p in self.p_list:
            print p

    ###########

    def update_velocity(self):
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

        for p in self.p_list:
            #velocity update with neighbors
            if self.num_neighbors > 0:
                v_x = inertia * p.velocity_x + cognition * rand1 * (p.personal_best[1] - p.x) + social_rate * rand2 * (Particle.global_best[1] - p.x) + local_rate * rand3 * (p.local_best[1] - p.x)
                v_y = inertia * p.velocity_y + cognition * rand1 * (p.personal_best[2] - p.y) + social_rate * rand2 * (Particle.global_best[2] - p.y) + local_rate * rand3 * (p.local_best[2] - p.y)

            #velocity update without neighbors
            #velocity' = inertia * velocity + c_1 * r_1 * (personal_best_position - position) + c_2 * r_2 * (global_best_position - position) 
            else: 
                v_x = self.inertia * p.velocity_x + self.cognition * rand1 * (p.personal_best[1] - p.x) + self.social_rate * rand2 * (self.Particle.global_best[1] - p.x)
                v_y = self.inertia * p.velocity_y + self.cognition * rand1 * (p.personal_best[2] - p.y) + self.social_rate * rand2 * (self.Particle.global_best[2] - p.y)

            #scale velocity
            #if abs(velocity) > maximum_velocity^2 
            #velocity = (maximum_velocity/sqrt(velocity_x^2 + velocity_y^2)) * velocity 
            if abs(v_x) > self.max_velocity:
                v_x2 = (self.max_velocity/math.sqrt(v_x**2 + v_y**2)) * v_x
                flag_x = True
            if abs(v_y) > self.max_velocity:
                v_y2 = (self.max_velocity/math.sqrt(v_y**2 + v_y**2)) * v_y
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

    def update_position(self):
        """update particle postions based on velocity"""
        for p in self.p_list:
            #position' = position + velocity' 
            p.x += p.velocity_x
            p.y += p.velocity_y

    ###########

    def update_personal_best(self):
        """at each epoch, check to see if each particle's current position
        is its best (or closest to the solution) yet
        """
        for p in self.p_list:
            #if(Q(position) > Q(personal_best_position)) 
            #personal_best_position = position 
            if Q(p.x,p.y) > p.personal_best[0]:
                p.personal_best = [Q(p.x,p.y), p.x, p.y]

    ###########

    def update_global_best(self):
        """find the best position of all the particles in the list"""
        tmp       = self.Particle.global_best
        tmp_index = self.Particle.best_index
        for p in self.p_list:
            #if(Q(position) > Q(global_best_position)) 
            #global_best_position = position 
            if Q(p.x,p.y) > tmp[0]:
                tmp       = [Q(p.x,p.y), p.x, p.y]
                tmp_index = p.index
        self.Particle.global_best = tmp
        self.Particle.best_index  = tmp_index

    ###########

    def update_local_best(self):
        """optionally find the best position out of a neighborhood"""
        tmp = [0.0, 0, 0]
        tmp_index = 0
        for p in self.p_list:
            for n in p.neighbors:
                #find the local best Q value
                if Q(self.p_list[n].x, self.p_list[n].y) > tmp[0]:
                    tmp = [Q(self.p_list[n].x, self.p_list[n].y), self.p_list[n].x, self.p_list[n].y]
                    tmp_index = self.p_list[n].index
            p.local_best       = tmp
            p.local_best_index = tmp_index
            #reset tmp
            tmp = [0.0, 0, 0]

    ###########

    def calc_error(self):
        """calculate the error at each epoch"""
        error_x = 0.0
        error_y = 0.0

        #for each particle p:
        #error_x += (position_x[k] - global_best_position_x)^2 
        #error_y += (position_y[k] - global_best_position_y)^2
        for p in self.p_list:
            error_x += (p.x - self.Particle.global_best[1])**2
            error_y += (p.y - self.Particle.global_best[2])**2

        #Then
        #error_x = sqrt((1/(2*num_particles))*error_x) 
        #error_y = sqrt((1/(2*num_particles))*error_y) 
        error_x = math.sqrt((1.0/(2.0*self.num_particles))*error_x)
        error_y = math.sqrt((1.0/(2.0*self.num_particles))*error_y)

        return [error_x, error_y]

    ###########

    def params_to_CSV(self):
        """put the parameters at the top of the CSV file"""
        f = open(self.fname, 'a+')
        f.write(('parameters\n'
        +'num_particles,inertia,cognition,social_rate,local_rate,world_width,world_height,max_velocity,max_epochs,num_neighbors\n'
        +str(self.num_particles)+','
        +str(self.inertia)+','
        +str(self.cognition)+','
        +str(self.social_rate)+','
        +str(self.local_rate)+','
        +str(self.world_width)+','
        +str(self.world_height)+','
        +str(self.max_velocity)+','
        +str(self.max_epochs)+','
        +str(self.num_neighbors)+
        '\n\n\nerror,over,time\n'+
        'x error,y error\n'))
        f.close()

    ###########

    def error_to_CSV(self, e):
        """print the error at each epoch to produce an error over time graph"""
        f = open(self.fname, 'a+')
        f.write(str(e[0])+','+str(e[1])+'\n')
        f.close()
            
    ###########

    def plot_to_CSV(self):
        """print the points at the end to create a scatter plot, or at each epoch
        to try for a gif animation
        """
        f = open(self.fname,'a+')
        f.write('\n\n\nfinal,coordinates\nx values,y values\n')
        for p in self.p_list:
            f.write(str(p.x)+','+str(p.y)+'\n')
        f.close()

########### distance functions

def mdist():
    global params
    return float(math.sqrt((params['world_width']**2.0)+(params['world_height']**2.0))/2.0)
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

particles = Particle_List(params)

########### main

if not a_test:
    particles.params_to_CSV()

epochs = 0
#######
while True:
    """each run through this loop represents and epoch"""
    ###
    particles.update_velocity()
    ###
    particles.update_position()
    ###
    particles.update_personal_best()
    ###
    particles.update_global_best()
    if params['num_neighbors'] > 0:
        particles.update_local_best()
    ###
    error = particles.calc_error()
    if not a_test:
        particles.error_to_CSV(error)
    if verbose:
        print error
    ###
    if error[0] < 0.01 and error[1] < 0.01:
        break
    elif epochs > params['max_epochs']:
        break
    ###
    epochs += 1
#######
particles.update_personal_best()
particles.update_global_best()
if params['num_neighbors'] > 0:
    particles.update_local_best()
particles.print_particles()    
if not a_test:
    particles.plot_to_CSV()
print 'epochs: ', epochs


