%% TwoTubeMexWithPsi Test Script
tube1.OD = 1.0e-3;
tube1.ID = 0.6e-3;
tube1.k1 = 25;
tube1.k2 = 60;
tube1.L = 150e-3;
tube1.Lc_transition = 25e-3;
tube1.Lt = 100e-3;
tube1.E = 50e9;
tube1.G = tube1.E/2/(1.4);

tube2.OD = 1.0e-3;
tube2.ID = 0.6e-3;
tube2.k1 = 25;
tube2.k2 = 50;
tube2.L = 100e-3;
tube2.Lc_transition = 20e-3;
tube2.Lt = 60e-3;
tube2.E = 50e9;
tube2.G = tube2.E/2/(1.4);


beta = [-100e-3,-60e-3];
[az,el]=view
time = 0;
figure;
for i=1:50
psiL = [2*pi*i/100,6*pi*i/100]
tic;
[R] = TwoTubeMexWithMultiCurve(tube1, tube2, psiL, beta, 20);
time = time + toc;
s_max = max(R.s);
s_min = min(R.s);
s_interp = linspace(s_min, s_max, 100);
p_interp = interp1(flipud(R.s), flipud(R.p), s_interp,'pchip');
plot3(p_interp(:,1), p_interp(:,2), p_interp(:,3))
daspect([1 1 1]);
xlim([-25 25]*1e-3)
ylim([-25 25]*1e-3)
zlim([-50 50]*1e-3)
% view(az,el)
grid on;
drawnow
end

disp(['Kinematics time: ' num2str(time)])

