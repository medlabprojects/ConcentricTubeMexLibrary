# ConcentricTubeMexLibrary

This library contains the kinematics for concentric tube robots with mex files for implementation in MatLab. 
This so far only supports robots with two through four tubes. For each of these, we have solutions for the general case, 
at S, and with Psi. 

The general case is (e.g. ThreeTube) solve for the kinematics directly given a rotation and translation at the base. 

AtS (e.g. ThreeTubeAtS) something where arclength is evaluated at specified point(s). 

With Psi (e.g. ThreeTubeWithPsi) takes into account torsion that would occur along the straight sections of the tubes. 
Each folder has a test file that gives a simple robot configuration and plots to test if the mex built correctly. 


### Note: this is code written by Hunter years ago, found on Margaret's network drive, then adapted and uploaded by Jesse to this repo
Some things might need tweaking later down the line but the main cpp functions are not changed. 

## How to use?
Each folder should have a test that makes a robot and runs the kinematic code, this is a good example to follow. 
In general, each will take in the same general set up. Let's use TwoTubeMex as an example. 
We use an object called a struct to hold the information for the robot; its similar to a class but without any methods. 
The mex function will take for the first n inputs (where n is the number of tubes), a struct with the information for each tube. 
The tube structure is the following:
```
OD: Outer diameter
ID: Inner diameter
k:  Tube precurvature
L:  Total tube length
Lt: Transmission length
E:  Young's Modulus
G:  Shear Modulus
```

The next two inputs are arrays with the base rotation (alpha or psi) and base translation (beta) for each tube.
 
Alpha: rotation at the base of the tubes

(Psi: rotation from torsion at each section)

Beta: distance between base of tube and the imaginary front plate (**must be negative**)



The output will be another struct with the following fields:
```
p_tip [3x1] - position of tip in base frame
q_tip [4x1] - quaternion orientation of tip in base frame
J_tip [6xnTubes] - jacobian of tip in base frame
s     [Nx1] - sampled backbone positions
p     [Nx3] - positions of each link
q     [Nx4] - orientations of each link
J     [Nx1] - cell containing 6x4 jacobians 
n (optional)- number of sampled points
```

Note: WithPsi will also contain the field
```
psi [Nx3]   - psi values at each link
```


Here is an example: 
```
tube1.OD = 1.0e-3;
tube1.ID = 0.6e-3;
tube1.k = 10;
tube1.Lt = 150e-3;
tube1.L = tube1.L + 33e-3;

tube1.E = 60e9;
tube1.G = tube1.E/2/(1.4);

tube2.OD = 1.0e-3;
tube2.ID = 0.6e-3;
tube2.k = 10;
tube2.Lt = 100e-3;
tube2.L = tube2.Lt + 60e-3;
tube2.E = 60e9;
tube2.G = tube2.E/2/(1.4);

alpha = [.2*pi, .3*pi];
beta = [-15e-3, -7e-3];
kin = TWoTubeMex(tube1, tube2, alpha, beta);
```



## Mex files

### What is a mex file?
A mex file is basically a way for matlab to take Cpp code and translate it into a function that you can use in a matlab script. 
Why would you want to do this? Cpp functions can generally be written to be much faster and extremely efficient. Also if all the 
code is already written in Cpp and you want to test it out in a matlab script, you can just turn it into a mex file instead of 
going through and re-writting everything by hand.

Macs and Windows both use different mex file extensions for some reason. Macs use the '.mexmaci64' and Windows use '.mexw64'
The respective one will be generated based on the machine running the code (i.e macs will automatically generate the mac mex files)

### Building the mex file
Mex files will contain all of the dependencies (i.e. any additional libraries or includes) that function needs to run, 
so they themselves can run independently, but you need all of these dependencies to **build** the mex file in the first place.

The TubeMex.cpp files need the ReferenceCannulaKinematics(AtS) folders to run. Those in turn rely on two Cpp libraries, Eigen and boost. 
These are both standard libraries and can be downloaded here for [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page#Download) and 
here for [Boost](https://www.boost.org). Put these in some folder on your drive. For example, I made a folder called 'CppLibraries' with
a subfolder called 'include' that now contains both Eigen and Boost. 

There are 2 ways to build the mex files. The first is in some folders, there is script called 'BuildMex.m'. Make sure the paths 
to the Eigen and Boost libraries are correct. Running this will generate the mex file. This might be using an old version of Cpp to compile though
so not recommended. 

The second way is with the MexMaker.m script. This is a bit more manual so we'll go through it line by line

```
mex -v CXXFLAGS='$CXXFLAGS --std=c++17 -fno-omit-frame-pointer'...
```

The first line calls mex and sets some flags (options)..

"-v" makes it verbose so you get all the warnings and outputs. 
"--std=c++17" defines which cpp version to use. Originally it was c++11 but c++17 is the latest version. Sometimes with these later versions
they change things slightly and there might be some syntax errors. It'll usually tell you how to fix them.
"-fno-omit-frame-pointer" I'm not really sure but it seems important 

```
 TwoTubeMex/WithPsi/TwoTubeMexWithPsi.cpp ...
```
This next line points to the make cpp file that you want to turn into a mex. 

--CHANGE THIS TO THE FILE YOU WANT TO COMPILE--

```
ReferenceCannulaKinematics/Rk8_Coeffs.cpp...
ReferenceCannulaKinematics/TaggedInterval.cpp...
```
These I'm not sure just leave as is.

```
-I./ReferenceCannulaKinematics...
-I/Users/jessedalmeida/CppLibraries/include/Eigen...
-I/Users/jessedalmeida/CppLibraries/include/boost
```
These point those dependencies I talked about earlier. ReferenceCannulaKinematics won't have to change since it should be in this same folder.
However, if you're doing an '__TubeMexAtS' file then you need to change this reference to be 'ReferenceCannulaKinematicsAtS' instead. 
**The last two should be the absolute path to where you saved the Eigen and Boost libraries**

Once you make these changes, run it and you should see a new mex file pop up!

## Notes
There will be a couple of warnings when you build the mex file. That's ok don't worry. 
Or maybe try to fix them if you're feeling ~spicy~ helpful. 

## Common Errors
### No C++ complier
You should have a C++ compiler somewhere on your device. The mex file will automcatically look for one. 
For instance though, I needed to download the XCode app (all 12GBs) from the App Store when I first got my mac and it included
the compiler Clang++. 










