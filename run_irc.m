function run_irc(vcDir_in, vcDir_out)
% arguments
% -----
% dir_in: input directory
% dir_out: output directory

if ~isdeployed() 
    addpath('./matlab'); 
    addpath('./mdaio');
    addpath('./jsonlab-1.5');
end

cSep = filesep();

p_ironclust(vcDir_out, ...
  [vcDir_in, cSep, 'raw.mda'], ...
  [vcDir_in, cSep, 'geom.csv'], ...
  'tetrode_template.prm', ...
  [vcDir_in, cSep, 'firings_true.mda'], ...
  [vcDir_in, cSep, 'firings_out.mda'], ...
  [vcDir_in, cSep, 'params.json']);

exit();
end %func
