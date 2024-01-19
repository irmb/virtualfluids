<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Debugging with VSCode

We can debug our VirtualFluids application within the docker container in VS Code. This wiki page describes how to do that.

1. Most important: VirtualFluids default build is __Release__. For debugging, the build type must be changed during the cmake process:
    ```
   cmake --preset=all_make -DCMAKE_BUILD_TYPE=Debug ..
    ```
    The variable ```CMAKE_BUILD_TYPE``` is a native cmake variable and used to specifiy the build type. The available options are: Debug, Release, RelWithDebInfo, MinSizeRel ([cmake-doc](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html)).

2. In VS Code we can now add a debug configuration file to our project. The debugging configuration information is stored in the file ```.vscode/launch.json```. To add the first debug configuration, we need to create this file manually.

    Then one possible configuration for gdb may look like this:
    ![devcontainer_vscode](img/vscode/launch-json.png)

    A file to start from can be found here: [launch.json](img/vscode/templates/launch.json).

3. A very important entry in this file is ```"program": ```.
Here we need to add the path to the executable we want to debug.
4. Now we can start debugging by adding a breakpoint to the specific file. To do this, we click on a line number on the left:

    ![devcontainer_vscode](img/vscode/breakpoint.png)

    The debugging can be started in the debugging tab:

    ![devcontainer_vscode](img/vscode/debugging-start.png)

    Here we can also monitor variables and the call stack:

    ![devcontainer_vscode](img/vscode/debugging.png)


Additionally VS Code provides a nice documentation as well:
- [Debugging in General](https://code.visualstudio.com/docs/editor/debugging)
- [C++ launch.json Reference](https://code.visualstudio.com/docs/cpp/launch-json-reference)
