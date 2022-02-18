
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
