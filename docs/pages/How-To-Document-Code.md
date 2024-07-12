<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# How To Document Code

To build sustainable research software, it is mandatory to document code. 
Even if it turns out that some developments are not continued, documentation is important to help future scientists to learn from the own experiences.

The documentation of the source code takes place…

- in commit messages  
  - Keeping commit messages concise, clear, and at a higher level of abstraction is a good practice for effective communication within a version control system like Git.
  - The commit message should briefly summarize what the commit does to the code. 

- in source code
  - VirtualFluids is using Doxygen to generate Documentation from within the source code
   - In most of the cases comment should describe ***why*** something was implemented and not ***how***.
   - if available add formulas, references to scripts, paper, and all information you got linked to the implemented code

- adding <!-- DOXYGEN_MAKE_REF -->[unit](Unit-Tests.md) and <!-- DOXYGEN_MAKE_REF -->[regression](Regression-Tests.md) tests.


## Doxygen
The [Doxygen](https://www.doxygen.nl/) documentation is automatically generated and deployed using a Gitlab-CI job. The pages are hosted on Gitlab page here: [VirtualFluids Documentation](https://irmb.gitlab-pages.rz.tu-bs.de/VirtualFluids/index.html).
The main configuration of the Doxygen Pages is done in the file `docs/doxygen/Doxyfile`.
For styling the pages are based on the [doxygen-awesome-css](https://jothepro.github.io/doxygen-awesome-css/) theme. 

The Documentation is build from two main parts:
### The source code

This part is sorted into groups using different special comments in the source code. For example the the file `src/lbm/collision/CollisionParameter.h` is sorted into the group `lbm` and the subgroup `collision`. This is done by the following comment in the source code:
```
//! \addtogroup collision
//! \ingroup lbm
//! \{

    // ... functions and variables and everything else ...

//! \}
 ```

### The markdown files in the docs/pages folder
The markdown files are automatically translated into html files. To add the files also in the sidebar, the file `docs/doxygen/DoxygenLayout.xml` needs to be updated.