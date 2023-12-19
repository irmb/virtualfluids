# Getting Start with Docker

This page describes how to start using and developing VirtualFluids.

## 1. Requirements

### Platforms
We mainly support two platforms: Windows and Linux. We are aware that there are multiple ways of developing VirtualFluids. 
This tutorial will focus on the usage of Visual Studio Code (VS Code), a terminal and an Ubuntu Docker image. The reason for this is that this workflow is currently possible on both Windows and Linux. Additionally, the use of Docker has the great advantage that the installation of necessary software is reduced to a minimum. The Docker image already contains compiler, cmake, mpi, openmp, clang-format and much more ...

**Important Note:** To use the GPU part of VirtualFluids in Docker, at least **Windows 11** is required.

### Installation of the Required Software
### Git
As mentioned before, [Git](https://git-scm.com/) is necessary to get access to VirtualFluids. Therefore at least some basic knowledge in Git is required. If you are new to Git, this [free book](https://git-scm.com/book/en/v2) is a good starting point.
Git needs to be installed and accessible in your terminal:
```
git --version
   > git version 2.25.1
```

#### Cuda Driver
_You can skip this part if you are planning to only use the CPU version of VirtualFluids._

Download and install the current driver for your corresponding GPU from [here](https://www.nvidia.com/Download/index.aspx?lang=en-us).

**Important Note:** It is possible that on Ubuntu systems, the driver installed with apt-get does not work with Docker. In that case, download and install it manually.

#### Docker
[Docker](https://www.docker.com/products/docker-desktop/) is a free software which helps to virtualise the whole development environment of VirtualFluids. Please download and install [Docker Desktop](https://www.docker.com/products/docker-desktop/) (on Windows) or the Docker CLI on Linux.

Only on Windows:
 - The installation process of Docker Desktop will guide you through enabling WSL2 (Windows Subsystem for Linux) if you haven't already.
 - If you already have it installed, update it to the latest version and enable Settings - General - Use the WSL2 backed engine.
 - To be able to use the Docker CLI from inside WSL2 (not just from PowerShell/cmd), enable the integration in Settings - Resources - WSL INTEGRATION.

#### VS Code
[VS Code](https://code.visualstudio.com/) is a free and open source, very powerful editor. What makes it so powerful are its many [extensions](https://code.visualstudio.com/docs/editor/extension-marketplace).
1. download and install VS Code
2. install the following extension within VS Code: [Remote - Containers](vscode:extension/ms-vscode-remote.remote-containers). More infos can be found in this tutorial: [click](https://code.visualstudio.com/docs/remote/containers-tutorial).


## 2. Get, Build and Run it ...
Now we are done installing the necessary software, and we can start working with VirtualFluids.

1. Open a shell or PowerShell and clone the repository onto your machine:
```
git clone https://git.rz.tu-bs.de/irmb/virtualfluids.git
```
_Please notice that this is the open source version of VirtualFluids.
Also: Currently, it is not possible to access the Git server with an SSH key outside the TU Braunschweig network._

2. Open VS Code and open the folder containing VirtualFluids in VS Code.
   - the Docker information is stored here: .devcontainer/devcontainer.json
     - "image" points to the current Docker image of VirtualFluids, which is stored [here](https://git.rz.tu-bs.de/irmb/virtualfluids/container_registry).
     - **Important**: "runArgs": ["--gpus","all"] is only necessary when the GPU version is used.

3. Click: <kbd>CTRL</kbd> + <kbd>SHIFT</kbd> + <kbd>P</kbd> and type: Remote-Containers: Open Folder in Container ...
   - this process can take a couple of minutes when it is performed for the first time.
   - when everything works fine, you can see this in the bottom left corner: ![devcontainer_vscode](img/vscode/devcontainer_vscode.png)

4. Open a new terminal within VS Code: Menu -> Terminal -> New Terminal
5. Build VirtualFluids in the terminal. The option "-DCMAKE_CUDA_ARCHITECTURE" is only necessary when the GPU version is used (CMAKE_CUDA_ARCHITECTURE should correspond to the [compute capability](https://en.wikipedia.org/wiki/CUDA#GPUs_supported) of your GPU):
```
   mkdir build && cd build
   cmake --preset=all_make -DCMAKE_CUDA_ARCHITECTURE=70 ..
   make -j 8
```
6. Run the examples, e.g. LidDrivenCavity:
```
   ./bin/LidDrivenCavityGPU
```
7. The result files of this simulation are stored in: ./build/output/DrivenCavity

8. The result files of VirtualFluids are mostly in [VTK](https://kitware.github.io/vtk-examples/site/VTKFileFormats/) format. These files can be visualised with the free software [Paraview](https://www.paraview.org/).

## 3. Advanced Topics

### Using configuration files
A developer of VirtualFluids can use a "configuration file" to specify some individual settings for the creation of the software. The files are located in ```CMake/cmake_config_files/```. The naming scheme of those files is ```<hostname>.config.cmake```. If the hostname of your computer is correct, this file is then loaded during the cmake process.

Things which can be defined in such a file could be
- the cuda architecture of the used gpu: e.g. ```set(CMAKE_CUDA_ARCHITECTURES 60)```,
- an application which is under a frequent use or development: e.g. ```list(APPEND USER_APPS "apps/gpu/LBM/WTG_RUB")```.

So far, we are tracking those files in our git repository. We will probably change this in the future, as these files only make sense on the individual computers.


**Configuration file together with docker:**
When a docker container is used, the hostname defaults to be the container id. In this case, the hostname needs to be changed:
- if the container is started manually with ```docker run```, the hostname can be passed as a parameter with ```--hostname <hostname>```
- when docker is used together with VS Code, the hostname can be set in ```.devcontainer.json```

```
 "runArgs": ["--hostname=${localEnv:HOSTNAME}"], // HOSTNAME needs to be known by the VS Code environment. It is probably necessary to add "export HOSTNAME=<hostname>" to the config file of your host machine's bash.
```
In this case it can be necessary to set the hostname in the config file of the bash as well (e.g. .bashrc).
```
export HOSTNAME=<hostname>
```


### Adding mounts to your Dev Container
[Bind mounts](https://docs.docker.com/storage/bind-mounts/) can be used to persist data or to provide additional data into containers. For example, you may want to provide STL files to VirtualFluids or output simulation results into the file system of the host machine.

When you use a bind mount, a file or directory on the host machine is mounted into a container. The file or directory is referenced by its absolute path on the host machine (```source```) and the path where the file or directory is mounted in the container (```target```).
When working with a [Dev Container](https://code.visualstudio.com/remote/advancedcontainers/add-local-file-mount), bind mounts are added to the ```devcontainer.json``` file:

```json
{
    "name": "virtual-fluids-env",
    // ...
    "mounts": [
        "source=/local/source/path/on/host/machine/,target=/target/path/in/container/,type=bind",
    ]
}
```
**Important Note:** Do not use spaces in the ```mount``` property.

You can reference local environment variables or the local path of the workspace. In this example, a local directory in the host's home folder is mounted into the container's workspace folder:
```json
{
    "name": "virtual-fluids-env",
    // ...
    "mounts": [
        "source=${localEnv:HOME}/STLs,target=${containerWorkspaceFolder}/stl,type=bind"
    ]
}
```
A documentation of all variables which can be referenced in ```devcontainer.json``` can be found [here](https://code.visualstudio.com/docs/remote/devcontainerjson-reference#_variables-in-devcontainerjson).

**For Windows users**:
On Windows, ```${localEnv:HOME}``` is located in the WSL Linux file system. In the file explorer, you can find it here:  ```\\wsl.localhost\Ubuntu-20.04\home\<username>\```. It is also possible to bind folders from the Windows file system, but this can have a negative impact on performance  ([Docker WSL 2 Best Practices](https://docs.docker.com/desktop/windows/wsl/#best-practices)).