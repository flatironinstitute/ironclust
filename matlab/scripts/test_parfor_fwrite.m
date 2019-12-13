vcFile = 'test_parfor.bin';
fid_w = fopen(vcFile, 'w');
n_repeat = 1000;
tic
for i=1:10
    fseek(fid_w, (i-1)*8*n_repeat, 'bof');
    fwrite(fid_w, repmat(i, [n_repeat, 1]), 'double');
end
toc
fclose(fid_w);

%%
fid_r = fopen(vcFile, 'r');
a=fread(fid_r, inf, '*double');
fclose(fid_r);