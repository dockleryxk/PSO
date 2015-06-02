# PSO
#### Overview
From [Wikipedia](http://en.wikipedia.org/wiki/Particle_swarm_optimization):
```
In computer science, particle swarm optimization (PSO) is a computational method 
that optimizes a problem by iteratively trying to improve a candidate solution 
with regard to a given measure of quality. PSO optimizes a problem by having a 
population of candidate solutions, here dubbed particles, and moving these 
particles around in the search-space according to simple mathematical formulae 
over the particle's position and velocity. Each particle's movement is influenced 
by its local best known position but, is also guided toward the best known 
positions in the search-space, which are updated as better positions are found 
by other particles. This is expected to move the swarm toward the best solutions.
```
#### Parameters
The following is some guidance for the ranges of the parameters. 

#### Number of Particles
A typical range is 20-40. For many problems, as few as 10 particles may be 
enough. For larger or more difficult problems, 100, 200, or more particles 
might be more appropriate. 

#### Inertia
Generally the range is 0-1, typically very close to 1. 

#### Cognition, Social, and Local parameters
Usually they are nearly equal and typically around 2, but they can range 
from 0-4. If the local parameter is set to 0, then the neighbors take no 
part in the optimization.

#### World Width and World Height
The range in the x and y directions, respectively. A good range is 
-50 to 50 in both directions. In that case, both the world width and height 
would be 100. 

#### Maximum Velocity
Limits how much a particle can move for a given iteration. Typical values 
to try are 1, 2, and 10. Some implementations define maximum velocity as the 
world width or height. 

#### Usage
swarm.py [-v] [test] numParticles inertia cognition socialRate localRate 
worldWidth worldHeight maxVelocity maxEpochs k fname

* -v invokes the verbose option, which will print the error as the algorithm
runs in the format [error_x, error_y]
* using the test argument will automatically provide the below parameters, -v works
    * Number of Particles = 20
    * Inertia             = 0.95
    * Cognition           = 2
    * Social Rate         = 2
    * Local Rate          = 2
    * World Width         = 100
    * World Height        = 100
    * Max Velocity        = 2
    * k                   = 0 (no neighbors taken into account)
    * Max Epochs          = 10,000
* no CSV file will be generated, which is normally specified by fname
Note: do not include an extension for fname, .csv is appended automatically

#### Additional Info
The problems solved in this particular program are as show:
![Problems](http://i.imgur.com/G6yKPUZ.png)
These are hard coded in, and only Problem 1 remains uncommented, however
anyone could write their own problem and make it into a Q function.

To give a better idea of what these particles are doing, here are two gifs
of the particles coverging at the global max (20, 7) and the local max
(-20, -7) over time, respectively.

![Global Max](http://i.imgur.com/C6EIyyZ.gif) ![Local Max](http://i.imgur.com/hVmY8DB.gif)

Further, here is a line plot of the average error of the particles over time.
This run was using the test parameters for problem 1:

![error over time](http://i.imgur.com/qgeQQ4k.jpg)

It is observed that the particles quickly near the solution, but "bounce around"
until they have achieved an error of less that 0.01 each. Riveting.
