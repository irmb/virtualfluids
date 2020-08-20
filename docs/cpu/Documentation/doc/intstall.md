# Installation    {#installation}
## Step 1: CMake  
Download and install [CMake](http://www.cmake.org). 
If you have Windows OS oder X-Windows for Linux/Unix start *cmake-gui*.
If you work in terminal start *ccmake*.
## Step 2: project configuration
You need set a path for *VirtualFluids* source code. It is allays *source* directory at the end of path.
E.g. *c:/vf/source*. You need also set a path to build binaries. E.g. *c:/vf/bin*. If you click **Configure** button start a project wizard and you can choose a compiler.
Set **Grouped** and **Advanced** check boxes in CMake. In the group *USE* you have follow important options:
* USE_BOND    - using an agent based communication framework for fault tolerance BOND
* USE_METIS   - using a domain decomposition tool METIS    
* USE_MPI       - using MPI for distributed memory parallelization   
* USE_YAML     - using YAML
* USE_ZOLTAN - using domain decomposition tool ZOLTAN  

There are externals library and you need additionally compile them.    
## Step 4: project compilation
Compile you project with suitable compiler