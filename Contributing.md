# Contributing
If you want to contribute to VirtualFluids, your help is very welcome.
To contribute use a merge request as follows:

## How to make a clean merge request

- Create a personal fork of VirtualFluids.
- Clone the fork on your local machine. Your remote repo on gitlab is called `origin`.
- Add the original repository as a remote called `upstream`.
- If you created your fork a while ago be sure to pull upstream changes into your local repository.
- Create a new branch to work on! Branch from `develop`
- Implement/fix your feature, comment your code.
- Follow the code style of the project, including indentation. We provide files to help you with automatic formatting for clang-format and cmake-format.
- Run the <!-- DOXYGEN_MAKE_REF -->[tests](docs/pages/unit-tests).
- Write or adapt tests (<!-- DOXYGEN_MAKE_REF -->[Unit-Tests](docs/pages/unit-tests) and <!-- DOXYGEN_MAKE_REF -->[Regression-Tests](docs/pages/regression-tests)) as needed.
- Add or change the <!-- DOXYGEN_MAKE_REF -->[documentation](docs/pages/how-to-document-code) as needed.
- Push your branch to your fork on gitlab, the remote `origin`.
- From your fork open a merge request in the correct branch. Target the project's `develop` branch
- …
- If we requests further changes just push them to your branch. The MR will be updated automatically.
- Once the merge request is approved and merged you can pull the changes from `upstream` to your local repo and delete
your extra branch(es).

