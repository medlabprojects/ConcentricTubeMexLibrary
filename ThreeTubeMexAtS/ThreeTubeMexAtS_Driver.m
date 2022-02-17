%
%BEGIN BUILD
%
%The Eigen library location, so that the
% compiler can find the include <Eigen/Dense>
% EIGEN_LOCATION='E:\Libraries\include';
EIGEN_LOCATION='C:\Libraries\C++\Eigen';


%The Boost library location
% BOOST_LOCATION='E:\Libraries\include'; 
BOOST_LOCATION='C:\Libraries\C++\boost_1_64_0';

%
% Shouldn't need to change anything below here
%
ipath_kin = ['-I"' fullfile(pwd,'..','ReferenceCannulaKinematicsAtS') '"'];
ipath_eigen = ['-I"' EIGEN_LOCATION '"'];
ipath_boost = ['-I"' BOOST_LOCATION '"'];

%Build the program with 'mex'
%(by default optimizations are on)
mex('-v', ...
     ipath_kin, ...
     ipath_eigen, ...
     ipath_boost, ...
     'ThreeTubeMexAtS.cpp', ...
     '../ReferenceCannulaKinematicsAtS/Rk8_Coeffs.cpp', ...
     '../ReferenceCannulaKinematicsAtS/TaggedInterval.cpp');
%
%END BUILD
%
%%

