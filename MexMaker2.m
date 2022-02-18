%% MexMaker2

% ----------CAN CHANGE----------
% The Eigen library location, so that the compiler can find the include <Eigen/Dense>
EIGEN_LOCATION='/Users/jessedalmeida/CppLibraries/include/Eigen';

% The Boost library location
BOOST_LOCATION='/Users/jessedalmeida/CppLibraries/include/boost';

% Target file
file = 'TwoTubeMex';

% withPsi?
withPsi = false;

% AtS? 
atS = false;

% -------- do not change please (unless good reason then go ahead)---------
if withPsi && atS
    disp('Cannot have withPsi and AtS, please choose one or the other')
end

% Shouldn't need to change anything below here
ipath_kin = ['-I"' fullfile(pwd,'..','ReferenceCannulaKinematics') '"'];
ipath_eigen = ['-I"' EIGEN_LOCATION '"'];
ipath_boost = ['-I"' BOOST_LOCATION '"'];

psi = '';
s = '';
if withPsi
    psi = 'WithPsi';
end
if atS
    s = 'AtS';
end

fname = [file 'Mex'  psi s ];

targetpath = [fname '' fname];
cxxflag = 'CXXFLAGS="\$CXXFLAGS" --std=c++17 -fno-omit-frame-pointer';

%Build the program with 'mex' (by default optimizations are on)
mex('-v', cxxflag, ...
     ipath_kin, ...
     ipath_eigen, ...
     ipath_boost, ...
     'TwoTubeMexWithPsi.cpp', ...
     '../ReferenceCannulaKinematics/Rk8_Coeffs.cpp', ...
     '../ReferenceCannulaKinematics/TaggedInterval.cpp');

%END BUILD
%