# Contributing

If you want to contribute to VirtualFluids, your help is very welcome.
To contribute use a pull request as follows:

### How to make a clean pull request

- Create a personal fork of VirtualFluids.
- Clone the fork on your local machine. Your remote repo on gitea is called `origin`.
- Add the original repository as a remote called `upstream`.
- If you created your fork a while ago be sure to pull upstream changes into your local repository.
- Create a new branch to work on! Branch from `develop` or `open_source`.
- Implement/fix your feature, comment your code.
- Follow the code style of the project, including indentation.
- Run the tests.
- Write or adapt tests as needed.
- Add or change the documentation as needed.
- Push your branch to your fork on gitea, the remote `origin`.
- From your fork open a pull request in the correct branch. Target the project's `develop` or `open_source` branch
- …
- If we requests further changes just push them to your branch. The PR will be updated automatically.
- Once the pull request is approved and merged you can pull the changes from `upstream` to your local repo and delete
your extra branch(es).

And last but not least: Always write your commit messages in the present tense. Your commit message should describe what the commit, when applied, does to the code – not what you did to the code.

## Documentation

To build sustainable research software, it is mandatory to document code. 
Even if it turns out that some developments are not continued, documentation is important to help future scientists to learn from the own experiences.  

The documentation of the source code takes place…

- in commit messages  
  - As it is possible to put all the information into the commit messages, we want to keep the messages short and on a higher level of abstraction.
  - The commit message should briefly summarize what the commit does to the code. 

- in source code
  - VirtualFluids is using Doxygen to generate Documentation from within the source code
   - In most of the cases comment should describe ***why*** something was implemented and not ***how***.
   - if available add formulars, references to scripts, paper, and all information you got linked to the implemented code
