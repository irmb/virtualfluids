<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Docker Virtualisation

The development environment for VirtualFluids is based on Docker. Docker is a tool that allows you to run applications in containers. A container is a lightweight, standalone, executable package of software that includes everything needed to run an application: code, runtime, system tools, system libraries, and settings.

The description of the VirtualFluids Docker image can be found in the folder `Containers` in the root folder of the repository. The newest Docker image is based on the `ubuntu:22.04` image and contains all dependencies to run VirtualFluids. 

Every user of VirtualFluids can use the Dockerfile to create a new Docker image. However, a prebuild Docker image is available in the Gitlab CI Docker Registry here: [https://git.rz.tu-bs.de/irmb/VirtualFluids/container_registry](https://git.rz.tu-bs.de/irmb/VirtualFluids/container_registry).

The devcontainer for Visual Studio Code and most oft the GitLab-CI jobs are based on the prebuild Docker image.

## Update the prebuild image
To update the prebuild image one need to have at least Maintainer rights in the Gitlab project. Additionally, at this point in time it is only possible to push to the registry from our Linux Gitlab Runner!!

First login to the Gitlab registry with the following command:
```bash
docker login git.rz.tu-bs.de:4567
```
To authenticate you need a personal access token. The token can be created in the Gitlab settings under `Access Tokens`. The token needs the `read_registry` and `write_registry` scope.

A new image with the tag name `vf_base' can be created from the Containers folder with the following command:
```bash
docker build -t vf_base -f Ubuntu22_04.Dockerfile .
```

We recommend to use the `vf_base` tag name for the first prebuild image and to create two new tags for the new image. To push the tags afterwards, the URL of the registry needs to part of the tag. Additionally, the first tag should be contain new version and the second tag contain the suffix `latest`. All current versions can be found [here](https://git.rz.tu-bs.de/irmb/VirtualFluids/container_registry/116). The following commands can be used to create the new tags:
```bash
docker image tag vf_base git.rz.tu-bs.de:4567/irmb/virtualfluids/ubuntu22_04:latest
docker image tag vf_base git.rz.tu-bs.de:4567/irmb/virtualfluids/ubuntu22_04:1.3
```

and to push the tags to the registry:
```bash
docker image push --all-tags git.rz.tu-bs.de:4567/irmb/virtualfluids/ubuntu22_04
```

Afterwards one can check if the image is available [here](https://git.rz.tu-bs.de/irmb/VirtualFluids/container_registry/116).
Then the new image can be used in the Gitlab CI jobs and the devcontainer for Visual Studio Code.