# TTK-Clemson

### Installing TTK and ParaView

These instructions are intended for Linux or Mac users. We strongly discourage to develop under windows.

The first step requires you to clone the group repository from GitHub. If you are already familiar with git or GitHub you can use the method that you prefer. If you are not you can follow these simple steps:

- register a GitHub account
- ask an administrator of the TTK-Clemson repository to be added to it
- install [GitKraken](https://www.gitkraken.com) on your computer (if you are using MacOS you can also use [GitHub Desktop](https://desktop.github.com)
- in GitKraken select `Clone repository` to create a local copy of the ttk-clemson repository. We will come back to the use of GitKraken later.


Now, follow the instructions [here](https://topology-tool-kit.github.io/installation.html) for installing Paraview and ttk. Remember to select the instructions for your operating system, and for version 0.9.7 of TTK along with the version 5.6.0 of ParaView. Use the Paraview and TTK folders found in the GitHub repository instead of downloading them from scratch.


### Rules for GitHub

At this point you should have a functioning version of Paraview and TTK. That means that if you open Paraview, you are able to see and load TTK modules. If this is the case you are almost ready to start creating your first module.

One last thing to learn is how to correctly use the GitHub repository. The repository that you have downloaded now is called `main`. Ideally, this part of the repository should always contain a fully functioning and stable version of TTK (i.e., a version where all the modules are stables and compilable). Of course this is rarely the case when we are developing some new module. To avoid issues please follow these steps:

- before creating a new module, create a branch in GitKraken and name it as your new module. For example, if you want to create a module `ScalarFieldSmoother` go in GitKraken (or GitHub desktop) and create a new branch there with the name `develop-ScalarFieldSmoother`.
- from this point you will be able to develop the new module without affecting the work of the other people.
- once you are done with developing the module and everything is fully functioning you can use `merge` to merge your work back into the main repository.


### Creating a new module in TTK-Clemson and using GitHub for keeping track of your work

- create a new branch in GitKraken with the name `develop-nameModule`
- using the terminal get into the ttk-0.9.7 folder and create a new module with the command `scripts/createTTKmodule.sh nameModule`
- look at GitKraken now. You will notice that it is keeping track of all the new files that have been generated. This is the main reason for using git. Git provides a service for keeping track of your work and all the features in your code. Every time you finish to update a functionality of your module you should upload your updates committing them to the repository. For a quick overview of the main idea behind git you can look at [this guide](https://support.gitkraken.com/start-here/guide/) are at many other resources on the web
