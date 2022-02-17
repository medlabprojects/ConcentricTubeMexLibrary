mex -v CXXFLAGS='$CXXFLAGS -std=c++11 -fno-omit-frame-pointer' ThreeTubeMex.cpp...
    ../ReferenceCannulaKinematics/Rk8_Coeffs.cpp...
    ../ReferenceCannulaKinematics/TaggedInterval.cpp...
    -I../ReferenceCannulaKinematics...
    -I"C:/Libraries/include"...
    -I"C:/Libraries/include/Eigen"
%%
tube1.OD = 1.0e-3;
tube1.ID = 0.6e-3;
tube1.k = 10;
tube1.L = 150e-3;
tube1.Lt = 100e-3;
tube1.E = 60e9;
tube1.G = tube1.E/2/(1.4);

tube2.OD = 1.0e-3;
tube2.ID = 0.6e-3;
tube2.k = 10;
tube2.L = 100e-3;
tube2.Lt = 60e-3;
tube2.E = 60e9;
tube2.G = tube2.E/2/(1.4);

tube3.OD = 1.0e-3;
tube3.ID = 0.6e-3;
tube3.k = 10;
tube3.L = 50e-3;
tube3.Lt = 25e-3;
tube3.E = 60e9;
tube3.G = tube3.E/2/(1.4);

psiL = [0,0,0];
beta = [-120e-3,-80e-3,-40e-3];

tic;
for i=1:1000
[R] = ThreeTubeMex(tube1, tube2, tube3, psiL, beta);
end
toc;

R.p_tip
R.q_tip
R.J_tip