from pyfluids.kernel import LBMKernel, KernelType

kernel = LBMKernel(KernelType.BGK)
kernel.use_forcing = True
kernel.forcing_in_x1 = 1
kernel.forcing_in_x2 = 1
kernel.forcing_in_x3 = 1
kernel.set_forcing(1, 2, 3)
print(kernel)
