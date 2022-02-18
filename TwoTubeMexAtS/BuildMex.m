
%BEGIN BUILD

%pwd should return ......\TwoTubeMexAtS
%make sure this is the case

%The Eigen library location, so that the
% compiler can find the include <Eigen/Dense>
EIGEN_LOCATION='/Users/jessedalmeida/CppLibraries/include/Eigen';
% EIGEN_LOCATION='C:\Libraries\C++\Eigen';


%The Boost library location
BOOST_LOCATION='/Users/jessedalmeida/CppLibraries/include/boost'; 
% BOOST_LOCATION='C:\Libraries\C++\boost_1_64_0';


% Shouldn't need to change anything below here

ipath_kin = ['-I"' fullfile(pwd,'..','ReferenceCannulaKinematicsAtS') '"'];
ipath_eigen = ['-I"' EIGEN_LOCATION '"'];
ipath_boost = ['-I"' BOOST_LOCATION '"'];

%Build the program with 'mex'
%(by default optimizations are on)
mex('-v', ...
     ipath_kin, ...
     ipath_eigen, ...
     ipath_boost, ...
     'TwoTubeMexAtS.cpp', ...
     fullfile(pwd,'..','ReferenceCannulaKinematicsAtS','Rk8_Coeffs.cpp'), ...
     fullfile(pwd,'..','ReferenceCannulaKinematicsAtS','TaggedInterval.cpp')...
     );
     
%END BUILD
%
%%

%
% Example use
%
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


beta = [-100e-3,-60e-3];
% [az,el]=view
time = 0;
for i=1:1
psiL = [2*pi*i/100,6*pi*i/100];
tic;
[R] = TwoTubeMexAtS(tube1, tube2, psiL, beta, 0, .011111111111);
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
drawnow
end

disp(['Kinematics time: ' num2str(time)])

