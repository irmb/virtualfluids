# Unit Tests

This page describes how to add unit tests to VirtualFluids. VirtualFluids uses the C++ testing and mocking framework [GoogleTest](http://google.github.io/googletest/).

## 0. Test Structure in VirtualFluids
VirtualFluids is build upon multiple libraries `<library>` (e.g. [basics](https://git.rz.tu-bs.de/irmb/virtualfluids/-/tree/main/src/basics)). Every library can have a corresponding test executable. The test executable is called `<library>Test` and is created automatically by the CMake build system if the following two conditions are met:
1. The CMakeLists.txt of the libraries contains: `vf_add_tests()` (e.g. the basics library: [CMakeLists.txt](https://git.rz.tu-bs.de/irmb/virtualfluids/-/blob/main/src/basics/CMakeLists.txt))
2. The library contains a file following the naming convention: `<fileName>Test.cpp` (e.g: [StringUtilTest.cpp](https://git.rz.tu-bs.de/irmb/virtualfluids/-/blob/main/src/basics/StringUtilities/StringUtilTest.cpp))


## 1. Building and Running Tests
The tests can be built by running the following commands:
```
cmake .. -DBUILD_VF_UNIT_TESTS=ON
make
```
or by using the presets:
```
cmake .. --preset=all_make
make
```


The tests can then be run by executing the specific test executable:
```
./bin/<library>Test
```
or all test executables can be run with:
```
ctest
```

Additionally all test executables are automatically executed by our [continuous integration pipeline](https://git.rz.tu-bs.de/irmb/virtualfluids/-/pipelines) during each push to the repository.

## 2. Adding Tests

We will show you a simple example on how to add tests to VirtualFluids:
To add tests for the class ```ExampleClass``` declared in ```ExampleClass.h``` you need to create a file named ```ExampleClassTest.cpp```.

The following code block shows a simple test case:
```cpp
#include <gmock/gmock.h>
#include "ExampleClass.h"


auto RealEq = [](auto value) {
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value);
#else
    return testing::FloatEq(value);
#endif
};


TEST(ExampleClassTest, multiplyIntByReal)
{
    //arrange
    real number = 1.1;
    int multiplicator = 2;
    ExampleClass sut();

    // act
    result = sut.multiplyIntByReal(multiplicator, number);

    // assert
    EXPECT_THAT(result, RealEq(2.2));
}
```

The signature of the test contains the test suite's name along with the name of the test: `` TEST(TestSuiteName, TestName)``.

The first step inside the ``TEST()`` function is to **arrange** the variables and objects needed for the test.

The next step ist to **act** on the target behavior. The act step should cover the main thing to be tested. For example this could be calling a member function which should be tested.

The third and final step is to **assert** the expected outcome. The result or response of the act step is checked. The assertion(s) determine(s) whether the test passes or fails.

### Assertions with googletest
For the assert step googleTest provides two options: [EXPECT and ASSERT](http://google.github.io/googletest/reference/assertions.html). Upon failure, `EXPECT_` macros generate nonfatal failures and allow the current function to continue running, while `ASSERT_` macros generate fatal failures and abort the current function. The above example uses the ``EXPECT_THAT`` macro. With ``EXPECT_THAT()``  you can use Google Test's [predefined matchers](http://google.github.io/googletest/reference/matchers.html) from the testing namespace (for example ``testing::IsTrue()`` or ``testing::IsEmpty()`` or your own custom matchers. The example above uses the custom matcher ``RealEQ()`` for testing the equality of two real numbers (float or double).



## 3. Common Problems
When you test a class which depends on CUDA you may get a build error like this:
```
fatal error: cuda_runtime.h: No such file or directory
11 | #include <cuda_runtime.h>
```

To fix this problem, you need to specify CUDA as language in CMake by adding the following code to ``CMakeList.txt``:
```
if(BUILD_VF_UNIT_TESTS)
    set_source_files_properties(ExampleClassTest.cpp PROPERTIES LANGUAGE CUDA)
endif()

```
For VirtualFluidsGPU you can find the CMakeList here: ``VirtualFluids/src/gpu/VirtualFluids_GPU/CMakeLists.txt``.

<!-- Gleiches Problm wie in Punkt 2. Running the Tests -->

## 4. Further Information
 You can find further information on how to write tests in [GoogleTest User’s Guide](http://google.github.io/googletest/).