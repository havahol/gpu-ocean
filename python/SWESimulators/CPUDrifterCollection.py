# -*- coding: utf-8 -*-

"""
This python class implements a DrifterCollection living on the CPU.

Copyright (C) 2018  SINTEF ICT

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


from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import time

import Common
import BaseDrifterCollection

class CPUDrifterCollection(BaseDrifterCollection.BaseDrifterCollection):
    """
    Class holding the collection of drifters.
    """ 
    def __init__(self, numDrifters, observation_variance=0.1,
                 boundaryConditions=Common.BoundaryConditions(), 
                 domain_size_x=1.0, domain_size_y=1.0):
        """
        Creates a GlobalParticles object for drift trajectory ensemble.

        numDrifters: number of drifters in the collection, not included the observation
        observation_variance: uncertainty of observation position
        boundaryConditions: BoundaryConditions object, relevant during re-initialization of particles.    
        """
        
        # Call parent constructor
        super(CPUDrifterCollection, self).__init__(numDrifters,
                                         observation_variance=observation_variance,
                                         boundaryConditions=boundaryConditions,
                                         domain_size_x=domain_size_x, 
                                         domain_size_y=domain_size_y)
        
        # One position for every particle plus observation
        self.positions = np.zeros((self.numDrifters + 1, 2))
        
    def copy(self):
        """
        Makes an independent indentical copy of the current object
        """
    
        copyOfSelf = CPUDrifterCollection(self.numDrifters,
                                observation_variance = self.observation_variance,
                                boundaryConditions = self.boundaryConditions,
                                domain_size_x = self.domain_size_x, 
                                domain_size_y = self.domain_size_y)
        copyOfSelf.positions = self.positions.copy()
        
        return copyOfSelf
    
    
    
    ### Implementation of abstract GETs
    
    def getDrifterPositions(self):
        return self.positions[:-1,:].copy()
    
    def getObservationPosition(self):
        return self.positions[-1, :].copy()
    
    
    
    ### Implementation of abstract GETs
    
    def setDrifterPositions(self, newDrifterPositions):
        np.copyto(self.positions[:-1,:], newDrifterPositions) 
        # Signature of copyto: np.copyto(dst, src)
    
    def setObservationPosition(self, newObservationPosition):
        np.copyto(self.positions[-1,:], newObservationPosition)
        
    
    def drift(self, eta, hu, hv, H0, nx, ny, dx, dy, dt, \
              x_zero_ref, y_zero_ref, sensitivity=1, doPrint=False):
        """
        Using the eta, hu, hv and H0 fields on the given grid to move
        all particles dt forward in time.
        
        Function copied from notebook where it first was implemented. 
        
        Assumes halo of 2 ghost cells and periodic boundary conditions.
        """
        # Change positions by reference
        positions = self.positions

        totNumDrifters = self.positions.shape[0]
        # Loop over all drifters (drifters + obs)
        for i in range(totNumDrifters):
            if doPrint: print "---------- Particle " + str(i) + " ---------------"
            x0, y0 = positions[i,0], positions[i,1]
            if doPrint: print "(x0, y0): ", (x0,y0)

            # First, find which cell each particle is in
            cell_id_x = int(np.ceil(x0/dx) + x_zero_ref)
            cell_id_y = int(np.ceil(y0/dy) + y_zero_ref)
            

            if (cell_id_x < 0 or cell_id_x > nx + 4 or cell_id_y < 0 or cell_id_y > ny + 4):
                print "ERROR! Cell id " + str((cell_id_x, cell_id_y)) + " is outside of the domain!"
                print "\t\Particle position is: " + str((x0, y0))

            if doPrint: print "cell values in x-direction: ", ((cell_id_x-2-0.5)*dx, (cell_id_x-2+0.5)*dx)
            if doPrint: print "cell values in y-direction: ", ((cell_id_y-2-0.5)*dy, (cell_id_y-2+0.5)*dy)

            h = H0 + eta[cell_id_y, cell_id_x]
            u = hu[cell_id_y, cell_id_x]/h
            v = hv[cell_id_y, cell_id_x]/h

            if doPrint: print "Velocity: ", (u, v)

            x1 = sensitivity*u*dt + x0
            y1 = sensitivity*v*dt + y0
            if doPrint: print "(x1, y1): ", (positions[i,0], positions[i,1])

            positions[i,0] = x1
            positions[i,1] = y1


        # Check what we assume is periodic boundary conditions    
        self.enforceBoundaryConditions()
        #applyPeriodicBoundaryConditionsToParticles(positions, nx, ny, dx, dy)
    
    
    def enforceBoundaryConditions(self):
        """
        Enforces boundary conditions on all particles in the ensemble, and the observation.
        This function should be called whenever particles are moved, to enforce periodic 
        boundary conditions for particles that have left the domain.
        """
        
        if (self.boundaryConditions.isPeriodicNorthSouth() and self.boundaryConditions.isPeriodicEastWest()):
            # Loop over particles
            for i in range(self.getNumDrifters() + 1):
                x, y = self.positions[i,0], self.positions[i,1]

                x, y = self._enforceBoundaryConditionsOnPosition(x,y)

                self.positions[i,0] = x
                self.positions[i,1] = y
        else:
            # TODO: what does this mean in a non-periodic boundary condition world?
            #print "WARNING [GlobalParticle.enforceBoundaryConditions]: Functionality not defined for non-periodic boundary conditions"
            #print "\t\tDoing nothing and continuing..."
            pass
    
    
    
  