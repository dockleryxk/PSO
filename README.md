# PSO
## Overview
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
## Parameters
The following is some guidance for the ranges of the parameters. 

## Number of Particles
A typical range is 20-40. For many problems, as few as 10 particles may be 
enough. For larger or more difficult problems, 100, 200, or more particles 
might be more appropriate. 

## Inertia
Generally the range is 0-1, typically very close to 1. 

##Cognition, Social, and Local parameters
Usually they are nearly equal and typically around 2, but they can range 
from 0-4. If the local parameter is set to 0, then the neighbors take no 
part in the optimization.

## World Width and World Height
The range in the x and y directions, respectively. A good range is 
-50 to 50 in both directions. In that case, both the world width and height 
would be 100. 

## Maximum Velocity
Limits how much a particle can move for a given iteration. Typical values 
to try are 1, 2, and 10. Some implementations define maximum velocity as the 
world width or height. 


