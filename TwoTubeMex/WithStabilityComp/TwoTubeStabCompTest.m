%% Test stability and compliance mex file
tube1.OD = 1.0e-3;
tube1.ID = 0.6e-3;
tube1.k = 25;
tube1.L = 150e-3;
tube1.Lt = 100e-3;
tube1.E = 50e9;
tube1.G = tube1.E/2/(1.4);

tube2.OD = 1.0e-3;
tube2.ID = 0.6e-3;
tube2.k = 25;
tube2.L = 100e-3;
tube2.Lt = 60e-3;
tube2.E = 50e9;
tube2.G = tube2.E/2/(1.4);

alpha = [0 0];
beta = [-1e-3 -1e-3];

[kin] = TwoTubeMexWithStabilityComp(tube1, tube2, alpha, beta)