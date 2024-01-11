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

- [ ] Check if there are open issues or merge requests which are assigned to the current Release

## 1. Version Number and Authors

- [ ] Update the version number in the authors.cff
- [ ] Update the version number in the CMakeLists.txt
- [ ] Update the version number in the pyproject.toml
- [ ] Update the version number in the sonar-project.properties
- [ ] check authors in authors.cff
    - authors of the current Release were directly involved in the development of the current Release
    - authors of the current Release are listed in the authors.cff
    - people  which are not involved in the current Release are listed in AUTHORS.md

## 2. Update the Changelog

- [ ] Update the Changelog.md

## 3. Prepare Release

1. [ ] Merge the develop branch into main
2. [ ] Create a tag for the current Release with the version number
3. Tag and Main Branch are automatically mirrored to https://github.com/irmb/virtualfluids
4. When Zenodo sees the new tag on github, it automatically creates a new version on Zenodo 

## Repositories

- Main Repository: https://git.rz.tu-bs.de/irmb/virtualfluids
- Mirror: https://github.com/irmb/virtualfluids
- Zenodo: https://zenodo.org/records/10283049 (DOI: https://doi.org/doi/10.5281/zenodo.10283048)
