mex -v CXXFLAGS='$CXXFLAGS --std=c++17 -fno-omit-frame-pointer'...
    TwoTubeMexWithPsi/TwoTubeMexWithPsi.cpp ...
    ReferenceCannulaKinematics/Rk8_Coeffs.cpp...
    ReferenceCannulaKinematics/TaggedInterval.cpp...
    -I./ReferenceCannulaKinematics...
    -I/Users/jessedalmeida/CppLibraries/include/Eigen...
    -I/Users/jessedalmeida/CppLibraries/include/boost


%% Test
robot1.tube1.OD = 1.18e-3;
robot1.tube1.ID = 0.889e-3;
robot1.tube1.k = 1/(58.7e-3);
robot1.tube1.Lt = 54.6e-3;
robot1.tube1.L = robot1.tube1.Lt + 33.6-3;
robot1.tube1.E = 60e9;
robot1.tube1.G = robot1.tube1.E/2/(1.33);

robot1.tube2.OD = 1.8e-3;
robot1.tube2.ID = 1.61e-3;
robot1.tube2.k = 0;
robot1.tube2.Lt = 41.8e-3;
robot1.tube2.L = robot1.tube1.Lt + 33.3-3;
robot1.tube2.E = 60e9;
robot1.tube2.G = robot1.tube1.E/2/(1.33);

robot1.base = eye(4);

alphas = [0,0];
betas = [-2 -1];
kin = TwoTubeMexWithPsi(robot1.tube1, robot1.tube2, alphas, betas)