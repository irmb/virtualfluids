:: don't close cmd on error
if not defined in_subprocess (cmd /k set in_subprocess=y ^& %0 %*) & exit )
:: @ECHO OFF
mpiexec -n 2 C:\Users\Master\Documents\MasterAnna\VirtualFluids_dev\build\bin\Release\MusselOyster.exe