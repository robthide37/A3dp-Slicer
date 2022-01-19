# Step by Step Visual Studio 2019 Instructions

### Install the tools

Install Visual Studio Community 2019 from [visualstudio.microsoft.com/vs/](https://visualstudio.microsoft.com/vs/). Older versions are not supported as PrusaSlicer requires support for C++17.
Select all workload options for C++ and make sure to launch Visual Studio after install (to ensure that the full setup completes).

Install git for Windows from [gitforwindows.org](https://gitforwindows.org/)
If you haven't, I advise to use a gui like tortoisegit.
Download and run the exe accepting all defaults

### Download sources

Clone the respository.  To place it in C:\src\REPO_NAME (use Slic3r, SuperSlicer or PrusaSlicer as REPO_NAME), run:
```
c:> mkdir src
c:> cd src
c:\src> git clone https://github.com/supermerill/SuperSlicer.git
```

### Run the automatic build script

The script `build_win.bat` will automatically find the default Visual Studio installation, set up the build environment, and then run both CMake and MSBuild to generate the dependencies and application as needed. If you'd rather do these steps manually, you can skip to the [Manual Build Instructions](#manual-build-instructions) in the next section. Otherwise, just run the following command to get everything going with the default configs:

```
c:\src>cd c:\src\REPO_NAME
c:\src\REPO_NAME>build_win.bat -d=..\REPO_NAME-deps -r=console
```

The build script will run for a while (over an hour, depending on your machine) and automatically perform the following steps:
1. Configure and build [deps](#compile-the-dependencies) as RelWithDebInfo with `c:\src\REPO_NAME-deps` as the destination directory
2. Configure and build all [application targets](#compile-slic3r) as RelWithDebInfo
3. Launch the resulting `superslicer-console.exe` binary

You can change the above command line options to do things like:
* Change the destination for the dependencies by pointing `-d` to a different directory such as: `build_win.bat -d=s:\REPO_NAMEDeps`
* Open the solution in Visual Studio after the build completes by changing the `-r` switch to `-r=ide`
* Generate a release build without debug info by adding `-c=Release` or a full debug build with `-c=Debug`
* Perform an incremental application build (the default) with: `build_win.bat -s=app-dirty`
* Clean and rebuild the application: `build_win.bat -s=app`
* Clean and rebuild the dependencies: `build_win.bat -s=deps`
* Clean and rebuild everything (app and deps): `build_win.bat -s=all`
* _The full list of build script options can be listed by running:_ `build_win.bat -?`

### Troubleshooting

You're best off initiating builds from within Visual Studio for day-to-day development. However, the `build_win.bat` script can be very helpful if you run into build failures after updating your source tree. Here are some tips to keep in mind:
* The last several lines of output from `build_win.bat` will usually have the most helpful error messages.
* If CMake complains about missing binaries or paths (e.g. after updating Visual Studio), building with `build_win.bat` will force CMake to regenerate its cache on an error.
* After a deps change, you may just need to rebuild everything with the `-s=all` switch.
* Reading through the instructions in the next section may help diagnose more complex issues.

* For PSAPI_LIB NOTFOUND: find your psapi.dll, which is usually on `C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\um\x64`. You can set your env variables `WindowsSdkDir` to `C:/Program Files (x86)/Windows Kits/10` and `WindowsSDKVersion` to an available version (the name of a directory inside C:/Program Files (x86)/Windows Kits/10/include). You can also just add the directory path of your psapi.dll in your PATH env variable.
* If you have "STL fixing by the Netfabb service will not be compiled", then it's the same problem with `WindowsSdkDir` and `WindowsSDKVersion`.

# Manual Build Instructions

_Follow the steps below if you want to understand how to perform a manual build, or if you're troubleshooting issues with the automatic build script._

### Compile the dependencies.
Dependencies are updated seldomly, thus they are compiled out of the Slic3r source tree.
Go to the Windows Start Menu and Click on "Visual Studio 2019" folder, then select the ->"x64 Native Tools Command Prompt" to open a command window and run the following:
```
cd c:\src\REPO_NAME\deps
mkdir build
cd build
cmake .. -G "Visual Studio 16 2019" -DDESTDIR="c:\src\REPO_NAME-deps"

msbuild /m ALL_BUILD.vcxproj // This took 13.5 minutes on my machine: core I7-7700K @ 4.2Ghz with 32GB main memory and 20min on a average laptop
```

### Generate Visual Studio project file for Slic3r, referencing the precompiled dependencies.
Go to the Windows Start Menu and Click on "Visual Studio 2019" folder, then select the ->"x64 Native Tools Command Prompt" to open a command window and run the following:
```
cd c:\src\REPO_NAME\
mkdir build
cd build
cmake .. -G "Visual Studio 16 2019" -DCMAKE_PREFIX_PATH="c:\src\REPO_NAME-deps\usr\local"
```

Note that `CMAKE_PREFIX_PATH` must be absolute path. A relative path like "..\..\REPO_NAME-deps\usr\local" does not work.

### Compile Slic3r. 

Double-click c:\src\REPO_NAME\build\REPO_NAME.sln to open in Visual Studio 2019.
OR
Open Visual Studio for C++ development (VS asks this the first time you start it).

Select REPO_NAME_app_gui as your startup project (right-click->Set as Startup Project).

Run Build->Rebuild Solution once to populate all required dependency modules.  This is NOT done automatically when you build/run.  If you run both Debug and Release variants, you will need to do this once for each.

Debug->Start Debugging or press F5

Slic3r should start. You're up and running!

note: Thanks to @douggorgen for the original guide, as an answer for a issue 


