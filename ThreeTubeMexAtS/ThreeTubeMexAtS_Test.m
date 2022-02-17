tube1.OD = 1.0e-3;
tube1.ID = 0.6e-3;
tube1.k = 10;
tube1.L = 150e-3;
tube1.Lt = 100e-3;
tube1.E = 50e9;
tube1.G = tube1.E/2/(1.4);

tube2.OD = 1.0e-3;
tube2.ID = 0.6e-3;
tube2.k = 10;
tube2.L = 100e-3;
tube2.Lt = 60e-3;
tube2.E = 50e9;
tube2.G = tube2.E/2/(1.4);

tube3.OD = 1.0e-3;
tube3.ID = 0.6e-3;
tube3.k = 10;
tube3.L = 50e-3;
tube3.Lt = 25e-3;
tube3.E = 50e9;
tube3.G = tube3.E/2/(1.4);

psiL = [0,0,0];
beta = [-100e-3,-60e-3,-25e-3];

tic;
time = 0;
sEval = 0.000123456;
N = 42;

R = ThreeTubeMexAtS(tube1, tube2, tube3, psiL, beta, N, sEval);
time = time + toc;
s_max = max(R.s);
s_min = min(R.s);
s_interp = linspace(s_min, s_max, 100);
p_interp = interp1(flipud(R.s), flipud(R.p), s_interp,'pchip');

plot3(p_interp(:,1), p_interp(:,2), p_interp(:,3))
axis tight;
daspect([1 1 1]);
drawnow;

toc;

R.p_tip;
R.q_tip;
R.J_tip;

disp(['Kinematics time: ' num2str(time)])

