function [tbl_acc, tbl_prc, tbl_rcl, tbl_cpu, csSorter, vnUnits_acc1] = sf_json2tbl_study(S_json, snr_thresh)

cS_sar = S_json.StudyAnalysisResults;

[tbl_acc, tbl_prc, tbl_rcl, tbl_cpu, vnUnits_acc1] = deal([]);
csSorter = {};
nStudy = numel(cS_sar);
for iStudy = 1:nStudy
    cS_sr1 = cS_sar{iStudy}.sortingResults;
    vrSnr1 = cS_sar{iStudy}.trueSnrs;
    vl1 = vrSnr1 >= snr_thresh;
    vnUnits_acc1(iStudy) = sum(vl1);
    for iSorter = 1:numel(cS_sr1)        
        S_sorter11 = cS_sr1{iSorter};
        csSorter{iSorter} = S_sorter11.sorterName;
        if ~isempty(S_sorter11.accuracies)
            tbl_acc(iStudy, iSorter) = mean_(S_sorter11.accuracies(vl1));
            tbl_prc(iStudy, iSorter) = mean_(S_sorter11.precisions(vl1));
            tbl_rcl(iStudy, iSorter) = mean_(S_sorter11.recalls(vl1));
            tbl_cpu(iStudy, iSorter) = mean_(S_sorter11.cpuTimesSec);
        else
            tbl_acc(iStudy, iSorter) = nan;
            tbl_prc(iStudy, iSorter) = nan;
            tbl_rcl(iStudy, iSorter) = nan;
            tbl_cpu(iStudy, iSorter) = nan;
        end
    end
end %for

end %func

function a = mean_(b)

if isempty(b), a=nan; return; end
if iscell(b), b = cell2mat(b); end
a = nanmean(b);

end %func