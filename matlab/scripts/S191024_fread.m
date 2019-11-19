%% fread test

 vcFile = 'D:\recordings\BuzsakiLab\DanEnglish\amplifier.dat';
%vcFile = 'C:\tmp\groundtruth\hybrid_synth\drift_siprobe\rec_64c_1200s_11\raw.mda';
% fid = fopen(vcFile, 'r');
S_dir = dir(vcFile);
nBytes = S_dir.bytes;
nRead = 100;
nSkip = floor(nBytes/nRead);
vrT_read = zeros(1, nRead);
fid = [];
parfor iRead = 1:nRead
    t1=tic;
    switch 3
        case 3
            fid = fopen(vcFile, 'r');
            fseek(fid, (iRead-1)*nSkip, 'bof');
            fclose(fid); fid = [];
        case 2, a=fread(fid, nSkip, '*int8');
        case 1, fseek(fid, (iRead-1)*nSkip, 'bof');
    end
    t1=toc(t1);
    vrT_read(iRead) = t1;
%     disp(t1);
end
if ~isempty(fid)
    fclose(fid);
end
disp(mean(vrT_read));
hold on; plot(sort(vrT_read));