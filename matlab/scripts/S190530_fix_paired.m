%% subsample monotrode data, determine the SNR of the recordings
% put scale_factor
% convert data type for kampff which is originally in int16
% remove snr below 1


vcDir_in = '/mnt/ceph/users/jjun/groundtruth_original/paired_recordings';
vcDir_out = '/mnt/ceph/users/jjun/groundtruth/paired_recordings';
%csDir = {'boyden32c', 'crcns', 'kampff', 'mea64c'};
csDir = {'mea64c'};
files_per_study = 100;
min_snr = 1; % one GT per paired
scale_factor = [.195, 0.305176, 2.34, .1042];
csFiles_copy = {'firings_true.mda' ,'geom.csv', 'params.json', 'raw.mda'};
readmda_ = @(x)irc('call','readmda', {x});
writemda_ = @(x,y)irc('call','writemda', {x,y});
for iDir = 1:numel(csDir)
    scale_factor1 = scale_factor(iDir);
    vcDir_in1 = fullfile(vcDir_in, csDir{iDir});
    vcDir_out1 = fullfile(vcDir_out, csDir{iDir});
    csFiles_rec1 = irc('call', 'dir_', [vcDir_in1, '/**/raw.mda']);
    
    % select `files_per_study` number of recordings based on SNR
    if numel(csFiles_rec1) > files_per_study
        vrSnr_rec1 = zeros(size(csFiles_rec1));
        for iFile = 1:numel(csFiles_rec1)
            vcFile_gt_mda11 = strrep(csFiles_rec1{iFile}, 'raw.mda', 'firings_true.mda');
            S_gt1 = irc('call', 'mda2gt1_', {vcFile_gt_mda11});
            vrSnr_rec1(iFile) = S_gt1.vrSnr_min_clu(1);
        end %for
        vi_select_rec1 = find(vrSnr_rec1 > min_snr);    
        if numel(vi_select_rec1) > files_per_study
            [~, viSrt2] = sort(vrSnr_rec1(vi_select_rec1), 'ascend');
            vii1 = irc('call', 'subsample_vr', {viSrt2, files_per_study});
            vi_select_rec2 = vi_select_rec1(vii1);
        else
            vi_select_rec2 = vi_select_rec1;
        end
        csDir_rec2 = cellfun(@(x)fileparts(x), csFiles_rec1(vi_select_rec2), 'UniformOutput',0);
    else
        csDir_rec2 = cellfun(@(x)fileparts(x), csFiles_rec1, 'UniformOutput',0);
    end    
    
    for iDir2 = 1:numel(csDir_rec2)
        vcDir_in_rec2 = csDir_rec2{iDir2};
        vcDir_out_rec2 = strrep(vcDir_in_rec2, vcDir_in, vcDir_out); 
        irc('call', 'mkdir_', {vcDir_out_rec2});
        for iCopy = 1:numel(csFiles_copy)
            vcFile_in_ = fullfile(vcDir_in_rec2, csFiles_copy{iCopy});
            vcFile_out_ = fullfile(vcDir_out_rec2, csFiles_copy{iCopy});
            irc('call', 'delete_', {vcFile_out_});
            switch csFiles_copy{iCopy}
                case {'firings_true.mda', 'geom.csv'}
                    copyfile(vcFile_in_, vcFile_out_);
                case 'raw.mda'
                    switch csDir{iDir}
                        case {'kampff', 'crcns'} % int16, zero-centered step of 1
                            mr_ = int16(readmda_(vcFile_in_));
                            writemda_(mr_, vcFile_out_);
                        case 'mea64c' % uint16, 2^15 centered
                            mr_ = int16(int32(readmda_(vcFile_in_)) - 2^15);
                            writemda_(mr_, vcFile_out_);
                        case 'boyden32c' % single, zero-centered step of .195
                            mr_ = readmda_(vcFile_in_);
                            mr_ = int16((mr_-median(mr_,2)) / .195);
                            writemda_(mr_, vcFile_out_);
                        otherwise
                            copyfile(vcFile_in_, vcFile_out_);
                    end
                case 'params.json'
                    S_json_ = loadjson(vcFile_in_);
                    if ~isfield(S_json_, 'scale_factor')
                        S_json_.scale_factor = scale_factor1;
                        irc('call', 'struct2json_', {S_json_, vcFile_out_});
                    else
                        copyfile(vcFile_in_, vcFile_out_);
                    end                    
            end %switch
        end % for iCopy
    end % for iFile2
    fprintf('Saved %s: %d/%d recordings\n', csDir{iDir}, numel(csDir_rec2), numel(csFiles_rec1));
end %forreadmda