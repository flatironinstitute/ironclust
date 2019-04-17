function cellstr2file(cs, vcFile)
fid=fopen(vcFile,'w'); cellfun(@(x)fprintf(fid,'%s\n',x), cs); fclose(fid);
end %func