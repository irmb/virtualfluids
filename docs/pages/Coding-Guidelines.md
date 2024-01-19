<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Coding Guidelines

The current VirtualFluids Coding Guidelines are available as PDF here: [PDF Coding Guidelines](https://git.rz.tu-bs.de/irmb/codingguidelines/-/blob/master/main.pdf)

## Missing changes in Coding Guidelines:

### Filenames

Usual filenames must not contain the word "test", "tests" or "mock". Otherwise, the file is put into test target and not into production code by cmake.

### Names

We stick to camel case naming

| Type | Spelling | Example |
|------|----------|---------|
| Subfolder | CamelCase | BoundaryCondition |
| File | CamelCase | BoundaryCondition.h |
| Class | CamelCase | Grid |
| Variable | camelCase | velocity |
| Namespace | camelCase | collisionKernel |

### Name of Configuration Files

appname_zusatz.cfg

### Includes

We want to include self implemented files always starting with the library name. e.g.: `#include <lbm/constants/D3Q27.h>`

We distinguish between quotation marks ( "..." ) for self-implemented libraries and angle brackets ( \<...\> ) for external libraries.

Includes have to be grouped according to their abstraction (standard library / external libraries / internal libraries / same library)

Includes per group have to be in alphabetical order.

e.g.:

```cpp
// standard library
#include <iostream>
#include <string>

// external libraries
#include <cuda.h>
#include <cuda-runtime.h>

// internal libraries
#include <lbm/ChimeraTransformation.h>
#include <lbm/constants/D3Q27.h>

// same library
#include "cpu/core/BoundaryConditions/BCArray.h"
#include "cpu/core/BoundaryConditions/BCFunction.h"
```

```plaintext
```

---

```cpp
int i = 0;
```

gzadfghsdfh `variable` dfdfdfd