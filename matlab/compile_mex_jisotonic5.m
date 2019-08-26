path0=fileparts(mfilename('fullpath')); %directory of this script
mex([path0,'/jisotonic5_mex.cpp'],'-output',[path0,'/jisotonic5_mex']);
