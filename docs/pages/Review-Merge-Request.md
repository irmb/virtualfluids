<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Review Merge Request

This document describes the process of reviewing a merge request in VirtualFluids. The review process is an important part of the development process. The review process ensures that the code changes are correct, documented, and tested. We follow a 4-eye principle, which means that every merge request needs to be reviewed by at least one other developer. Thereby the knowledge of the code is spread among the developers and the code quality is improved.

## Review Process

The reviewer is responsible for checking the code changes in the merge request. Thereby it is important to keep in mind that the review process is not only about finding bugs but also about improving the code quality. The review-process should be seen as a learning process for both the author and the reviewer. The reviewer should provide constructive feedback to the author. While we have some recommendations for the review process, the review process is not a strict process. The reviewer can adapt the process to the specific merge request.

## Recommendations for the Review Process
 
- The code changes are correct and do not break the existing code
- The code changes are documented in the source code (doxygen in-code)
- The code changes are documented in the documentation (markdown pages)
- The code changes are tested with unit tests (existing or new tests)
- The code changes are tested with regression tests (existing or new tests)
- The code changes fullfill our coding guidelines (see  <!-- DOXYGEN_MAKE_REF -->[Coding Guidelines](Coding-Guidelines.md))
    - e.g. the variable names are meaningful and the code is readable in general
