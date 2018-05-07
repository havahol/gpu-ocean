import unittest
import time
import numpy as np
import sys
import gc

from testUtils import *

sys.path.insert(0, '../')
from SWESimulators import Common
from SWESimulators.CPUDrifterCollection import *
from driftTests.BaseDrifterTest import BaseDrifterTest

class CPUDrifterTest(BaseDrifterTest):

    def setUp(self):
        super(CPUDrifterTest, self).setUp()

    def tearDown(self):
        super(CPUDrifterTest, self).tearDown()
    

    def create_small_drifter_set(self):
        self.smallDrifterSet = CPUDrifterCollection(self.numDrifters,
                                                     self.observationVariance,
                                                     self.boundaryCondition)

    def create_resampling_drifter_set(self):
        self.resamplingDrifterSet = CPUDrifterCollection(self.resampleNumDrifters)

    def create_large_drifter_set(self, size, domain_x, domain_y):
        return CPUDrifterCollection(size, domain_size_x=domain_x, domain_size_y=domain_y)

    
    def issue_drift(self, drifters):
        eta, hu, hv = self.sim.download()
        Hi = self.sim.downloadBathymetry()[0]
        drifters.drift(eta, hu, hv, np.max(Hi), self.nx, self.ny, self.dx, self.dy, self.dt,
                       self.x_zero_ref, self.y_zero_ref)
    
