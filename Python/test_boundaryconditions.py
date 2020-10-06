from pyfluids.boundaryconditions import *

no_slip_algo = NoSlipBCAlgorithm()
no_slip_bc_adapter = NoSlipBCAdapter()
no_slip_bc_adapter.algorithm = no_slip_algo