addpath('./mdaio');
addpath('./JRCLUST');

temp_path = 'D:\JeremyMagland\K15\jrc3';
raw_fname = 'D:\JeremyMagland\K15\tet_K=15_1.mda';
prb_fname = '.\JRCLUST\prb\tetrode.prb';
param_fname = 'D:\JeremyMagland\K15\tet_K=15_1.meta';
settings_fname = 'D:\JeremyMagland\K15\settings.prm';
firings_out_fname = 'firings_out.mda';

% Clear previous 
try rmdir(temp_path, 's'); catch, end

p_jrclust(temp_path,raw_fname,prb_fname,param_fname,settings_fname,firings_out_fname);

