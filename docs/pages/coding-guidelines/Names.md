<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->
# Names

In general, names should be self-explanatory and as short as possible. We are using camelCase for naming variables and functions. The following table gives an overview of the naming conventions used in the code.

| Type | Spelling | Example |
|------|----------|---------|
| Subfolder | PascalCase | BoundaryCondition |
| File | PascalCase | BoundaryCondition.h |
| Class | PascalCase | GridBuilder |
| Variable | camelCase | inflowVelocity |
| Namespace | snake_case | collision_kernel |

Below are some more detailed guidelines for naming.

## Names representing types must be in mixed case starting with upper case (PascalCase).

Example: Line, SavingsAccount

Common practice in the C++ development community.

## Variable names and attribute names must be in mixed case starting with lowercase (camelCase).

Example:
```cpp
line, inflowVelocity
```

Common practice in the C++ development community. Makes variables easy to distinguish from types, and effectively resolves potential naming collision, as in the declaration `Line line`;

## Names representing methods or functions must be verbs and written in mixed case starting with lower case (camelCase).

Example: 
```cpp
getName(), computeTotalWidth()
```

Common practice in the C++ development community. This is identical to variable names, but functions in C++ are already distinguishable from variables by their specific form.

## Names representing namespaces should be snake_case.

Example: 
```cpp
model::analyzer, io::iomanager, common::math::geometry, vf::compare_nups
```

Common practice in the C++ development community

## All names must be written in English.

Example: 
```cpp
fileName; // NOT: dateiName
```

English is the preferred language for international development.

## Variables with a large scope should have long names, variables with a small scope can have short names.

Example: NA

Scratch variables used for temporary storage or indices are best kept short. A programmer reading such variables should be able to assume that their value is not used outside of a few lines of code. Common scratch variables for integers are i, j, k, m, n and for characters c and d.

## The name of the object is implicit and should be avoided in a method name.

Example: 
```cpp
line.getLength(); // NOT: line.getLineLength();
```

The latter seems natural in the class declaration, but proves superfluous in use, as shown in the example.

## The plural form should be used on names representing a collection of objects.

```cpp
vector<Point> points;
int neighbors[10];
```

Enhances readability since the name gives the user an immediate clue of the type of the variable and the operations that can be performed on its elements.

## Abbreviations in names should be avoided

Example:
```cpp
computeAverage(); //NOT: compAvg();
```

There are two types of words to consider. First are the common words listed in a language dictionary. These must never be abbreviated. Write:

* command instead of cmd
* copy instead of cp
* point instead of pt
* compute instead of cp
* initialize instead of init

etc.

Then there are domain-specific phrases that are more naturally known through their abbreviations/acronyms. These phrases should be kept abbreviated. Write:

* html instead of HypertextMarkupLanguage
* cpu instead of CentralProcessingUnit 

etc.

## Naming pointers specifically should be avoided.

Example:
```cpp
Line* line; // NOT: Line* pLine; // NOT: Line* linePtr;
```

Many variables in a C/C++ environment are pointers, so a convention like this is almost impossible to follow. Also, objects in C++ are often oblique types where the specific implementation should be ignored by the programmer. Only when the actual type of object is of special significance, the name should emphasize the type.

## Negated boolean variable names must be avoided.

Example:
```cpp
bool isError; // NOT: isNoError
bool isFound; // NOT: isNotFound
```

The problem arises when such a name is used in conjunction with the logical negation operator as this results in a double negative. It is not immediately apparent what !isNotFound means.

