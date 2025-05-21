<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Release

## Release Policy

MAJOR.MINOR.PATCH
- Patch: Bugfixes, no new features
- Minor: New features
- Major: will be considered when VirtualFluids achieves a stable state


Release of VirtualFluids can be done by the core developers of VirtualFluids. A Release requires the following steps:

## 0. Check Issues and Merge Requests

- Check if there are open issues or merge requests that are assigned to the current Release
    - If there are open issues or merge requests, check if they are still relevant for the current Release
    - If they are still relevant, assign them to the next Release
    - If they are not relevant anymore, close them

## 1. Version Number and Authors

- Update the version number and the release date in the CITATION.cff
- Update the version number in the CMakeLists.txt
- Update the version number in the pyproject.toml
- Update the version number in the sonar-project.properties
- check authors in CITATION.cff
    - authors who were directly involved in the development of the current Release are listed in the CITATION.cff
    - people who are not involved in the current Release are listed in AUTHORS.md

## 2. Update the Changelog

- Update the Changelog.md

## 3. Prepare Release

1. Merge the develop branch into main branch on https://git.rz.tu-bs.de/irmb/VirtualFluids
    - this can be done by creating a Merge Request from develop to main [here](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/merge_requests/new?merge_request%5Bsource_branch%5D=develop&merge_request%5Btarget_branch%5D=main)
2. Create a tag on [gitlab](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/tags/new) from the main branch. Name of the tag can be the current version number.
3. Tag and Main Branch are automatically mirrored to https://github.com/irmb/virtualfluids
4. Create a new Release on [gitlab](https://git.rz.tu-bs.de/irmb/VirtualFluids/-/releases/new) and [github](https://github.com/irmb/virtualfluids/releases/new) based on the new tag
5. When Zenodo sees the new release on github, it automatically creates a new version on Zenodo. If everything went well, you can find it [here](https://zenodo.org/account/settings/github/repository/irmb/virtualfluids).
   - This step requires that the Zenodo account is linked to the github account and the repository is enabled on Zenodo. Learn how to link an account [here](https://help.zenodo.org/docs/profile/linking-accounts/).
    - You can find more information about integrating zenodo with github [here](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content).


## Repositories

- Main Repository: https://git.rz.tu-bs.de/irmb/virtualfluids
- Mirror: https://github.com/irmb/virtualfluids
- Zenodo: https://zenodo.org/records/10283049 (DOI: https://doi.org/doi/10.5281/zenodo.10283048)
