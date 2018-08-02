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
catch
    fprintf(2,lasterr());
end
%quit;
