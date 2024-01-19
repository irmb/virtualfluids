<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Release

## Release Policy

MAJOR.MINOR.PATCH
- Patch: Bugfixes, no new features
- Minor: New features
- Major: will be considered when VirtualFluids achieves a stable state


Release of VirtulFluids can be done by the core developers of VirtualFluids. A Release requires the following steps:

## 0. Check Issues and Merge Requests

- [ ] Check if there are open issues or merge requests that are assigned to the current Release
    - [ ] If there are open issues or merge requests, check if they are still relevant for the current Release
    - [ ] If they are still relevant, assign them to the next Release
    - [ ] If they are not relevant anymore, close them

## 1. Version Number and Authors

- [ ] Update the version number in the authors.cff
- [ ] Update the version number in the CMakeLists.txt
- [ ] Update the version number in the pyproject.toml
- [ ] Update the version number in the sonar-project.properties
- [ ] check authors in authors.cff
    - authors of the current Release were directly involved in the development of the current Release
    - authors of the current Release are listed in the authors.cff
    - people who are not involved in the current Release are listed in AUTHORS.md

## 2. Update the Changelog

- [ ] Update the Changelog.md

## 3. Prepare Release

1. [ ] Merge the develop branch into main branch
2. [ ] Create a tag for the current Release with the version number
3. Tag and Main Branch are automatically mirrored to https://github.com/irmb/virtualfluids
4. [ ] Create a new Release on gitlab and github based on the tag
5. When Zenodo sees the new release on github, it automatically creates a new version on Zenodo (this requires that the Zenodo account is linked to the github account and the repository is enabled on Zenodo)

## Repositories

- Main Repository: https://git.rz.tu-bs.de/irmb/virtualfluids
- Mirror: https://github.com/irmb/virtualfluids
- Zenodo: https://zenodo.org/records/10283049 (DOI: https://doi.org/doi/10.5281/zenodo.10283048)
