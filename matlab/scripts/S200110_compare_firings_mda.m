% sRateHz = 30000;
% load_firings_ = @(x)irc2('call', 'load_firings_', {x, sRateHz});
% S_irc = load_firings_('D:\Globus\recordings\kf19_nt27\irc2\firings.mda');


% S_ms4 = load_firings_('D:\Globus\recordings\kf19_nt27\kf19_nt27.firings_orig.mda');
% describe firings.mda

%% validate S_irc using S_ms4

% preselect from S_irc, remove durations that does not overlap
% skip clusters if there is no overlap what so ever

irc2 auto D:\Globus\recordings\kf19_nt27\irc2\raw_geom.prm

vcFile_gt_mda = 'D:\Globus\recordings\kf19_nt27\kf19_nt27.firings_orig.mda';
vcFile_clu_mda = 'D:\Globus\recordings\kf19_nt27\irc2\firings.mda';
irc2('compare-mda', vcFile_gt_mda, vcFile_clu_mda); % 3 hours

irc2('describe-mda', vcFile_gt_mda);
irc2('describe-mda', vcFile_clu_mda);