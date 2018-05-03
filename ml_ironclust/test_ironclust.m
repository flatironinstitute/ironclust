% addpath('./mdaio');
% addpath('./JRCLUST');
% 
% temp_path = 'D:\JeremyMagland\K15\jrc3';
% raw_fname = 'D:\JeremyMagland\K15\tet_K=15_1.mda';
% prb_fname = '.\JRCLUST\prb\tetrode.prb';
% param_fname = 'D:\JeremyMagland\K15\tet_K=15_1.meta';
% settings_fname = 'D:\JeremyMagland\K15\settings.prm';
% firings_out_fname = 'firings_out.mda';
% 
% % Clear previous 
% try rmdir(temp_path, 's'); catch, end
% 
% p_ironclust(temp_path,raw_fname,prb_fname,param_fname,settings_fname,firings_out_fname);
% 

%matlab -nodisplay -r 
% addpath('/home/jamesjun/Dropbox/flatiron/ironclust/ml_ironclust/'); 
addpath('./JRCLUST'); 
addpath('./mdaio'); 
try 
    if ispc()
        p_ironclust('c:\temp', 'D:\mountainsort\K15\raw.mda', 'D:\mountainsort\K15\geom.csv', './JRCLUST/tetrode_template.prm', 'firings.mda', 'D:\mountainsort\K15\params.txt');
    else
        p_ironclust('/tmp/mountainlab-tmp/tempdir__DPdLUY', ...
            '/tmp/mountainlab-tmp/output_da379cf4ea75964c8155cf72fab47cee7b12d48f_timeseries_out.mda', ...
            '/home/jamesjun/mlscript/data/geom.csv','','/tmp/mountainlab-tmp/output__firings_out.mda', ...
            '/tmp/mountainlab-tmp/tempdir__DPdLUY/argfile.txt');
    end
%     p_ironclust( ...
%         '/tmp/mountainlab-tmp/tempdir_5566f41e19_pQojAF', ...
%         '/tmp/mountainlab-tmp/output_f509f39b380e8484d5a14ebd2744754cb1ee88f3_timeseries_out', ...
%         '/tmp/mountainlab-tmp/output_98828a223d010792aef9f1c7da1e6df7b2f3f279_geometry_out', ...
%         '/home/jamesjun/mountainsort_examples/examples/001_sort_datasets/datasets/synth_K10/firings.mda', ...
%         '/tmp/mountainlab-tmp/tempdir_5566f41e19_pQojAF/argfile.txt'); 
%     p_ironclust('~/Downloads/K15a/tmp/', ...
%         '~/Downloads/K15a/raw.mda', ...
%         '~/Downloads/K15a/geom.csv', ...
%         'tetrodes.prm', ...
%         'firings.mda', ...
%         '');
%     p_ironclust('c:\temp', 'D:\mountainsort\K15\raw.mda', 'D:\mountainsort\K15\geom.csv', './JRCLUST/tetrode_template.prm', 'firings.mda', 'D:\mountainsort\K15\params.txt');
%     p_ironclust(temp_path, raw_fname, geom_fname, prm_fname, arg_fname)

catch
    fprintf(2,lasterr());
end
%quit;
