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

- adding <!-- DOXYGEN_MAKE_REF -->[unit](Unit-Tests.md) and <!-- DOXYGEN_MAKE_DOXYGEN_MAKE_REFABSOLUTE -->[regression](Regression-Tests.md) tests.
