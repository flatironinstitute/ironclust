%% density clustering

P = struct('mu1', 0, 'sd1', 1, 'mu2', 4, 'sd2', 1, 'n1', 100, 'n2', 100, 'knn', 30, 'fontSize', 16);
rng(0);
vrX1 = randn([P.n1,1])*P.sd1 + P.mu1;
vrX2 = randn([P.n2,1])*P.sd2 + P.mu2;
figure_ = @()figure('Unit','Normalized','Color','w','OuterPosition',[0 0 .5 .5],'Toolbar','none','Menubar','figure');
xlabel_ = @(x)xlabel(x,'FontSize',P.fontSize, 'Interpreter', 'tex');
ylabel_ = @(x)ylabel(x,'FontSize',P.fontSize, 'Interpreter', 'tex');

%%
vrEdge = P.mu1-5:.5:P.mu2+5;
%figure; hold on; plot(vrX1,'.'); plot(vrX2,'.');

hFig = figure_(); hold on; %histogram(vrX1, vrEdge); histogram(vrX2, vrEdge);
fh_ksdensity = @(x)ksdensity(x,'bandwidth',.3);
set(gcf,'Color','w');
[y_,x_]=fh_ksdensity(vrX1); area(x_,y_,'FaceAlpha',.5,'FaceColor','b');
[y_,x_]=fh_ksdensity(vrX2); area(x_,y_,'FaceAlpha',.5,'FaceColor','r');

xlabel_('x');
ylabel_('obs. freq.'); grid on;

%% pairwise density plot

viClu_gt = [repmat(2,P.n1,1); repmat(1,P.n2,1)];
vrX = [vrX1; vrX2];

% calc Rho
[mrD_knn, miD_knn] = pdist2(vrX, vrX, 'euclidean', 'Smallest', P.knn);
vrRho = 1./mrD_knn(end,:)';

% calc Delta
mrD = pdist2(vrX, vrX);
mrD(vrRho <= vrRho') = nan;
[vrDelta, viDelta] = min(mrD);
vrDelta = vrDelta' .* vrRho;
[vrDelta_srt, viDelta_srt] = sort(vrDelta, 'descend');
[imax1, imax2] = deal(viDelta_srt(1), viDelta_srt(2));
vrDelta(imax1) = max(vrDelta)*1.1;

hFig=figure_(); hold on; plot(vrRho, vrDelta, 'k.'); 
plot(vrRho(imax1), vrDelta(imax1), 'r*');
plot(vrRho(imax2), vrDelta(imax2), 'b*');
plot(get(gca,'XLim'), ones(1,2), 'r-');
grid on;
xlabel_('\rho := 1/d\_knn'); ylabel_('\delta := d\_min* / d\_knn');

%% detect knn collision
% viCl_clu = find(vrDelta > 1);
% figure; hold on; grid on;
% plot(vrRho, vrDelta, 'k.');
% plot(vrRho(viCl_clu(1)), vrDelta(viCl_clu(1)), 'b*');
% plot(vrRho(viCl_clu(2)), vrDelta(viCl_clu(2)), 'r*');


%% remove merged
miD_knn_clu = miD_knn(2:end,viCl_clu);
vrDelta_clu = vrDelta(viCl_clu);
vl_keep_clu = true(numel(viCl_clu), 1);
for iiClu=1:numel(viCl_clu)
    iClu = viCl_clu(iiClu);
    viClu_knn1 = miD_knn_clu(:,iiClu);
    viCl_clu1 = viCl_clu(vrDelta_clu >= vrDelta_clu(iiClu));
    if any(ismember(viClu_knn1, viCl_clu1))
        vl_keep_clu(iiClu) = false;
    end
end %for


%% plot nearesat links
vi_nn = viDelta(:);
mrX_nn = [vrX(:), vrX(vi_nn)];
mrY_nn = [vrRho(:), vrRho(vi_nn)];
mrX_nn(:,end+1) = nan;
mrY_nn(:,end+1) = nan;

hFig = figure_(); hold on; grid on; 
if 0
    plot(mrX_nn(:), mrY_nn(:), '-', 'Color', repmat(.75,1,3));
end
plot(vrX, vrRho, 'k.');
plot(vrX(imax1), vrRho(imax1), 'r*');
plot(vrX(imax2), vrRho(imax2), 'b*');
xlabel_('x');
ylabel_('\rho := 1/d\_knn');

% animate density line
[vrRho_srt, viRho_srt] = sort(vrRho, 'descend');
hLine = plot(get(gca,'XLim'), [nan, nan], 'k-');
hClu1 = plot(nan,nan,'ro','MarkerSize',3);
hClu2 = plot(nan,nan,'bo','MarkerSize',3);
hLink = plot(nan,nan, '-', 'Color', repmat(.75,1,3));
viClu = zeros(size(vrX));
viClu(imax1) = 1;
viClu(imax2) = 2;
[vrX_clu1, vrX_clu2, vrY_clu1, vrY_clu2, mrX_edge, mrY_edge] = deal([]);
hTitle = title('Accuracy: ', 'FontSize', P.fontSize);
vrT_wait = logspace(log10(.01), log10(.25), numel(vrX));
vrT_wait = vrT_wait(end:-1:1);
for iiSpk = 2:numel(vrX)
    [iSpk1, rho1] = deal(viRho_srt(iiSpk), vrRho_srt(iiSpk));
    x1 = vrX(iSpk1);
    hLine.YData = repmat(rho1, 1, 2);
    iSpk1_parent = viDelta(iSpk1);
    if viClu(iSpk1) == 0
        iClu1_parent = viClu(iSpk1_parent);
        switch iClu1_parent
            case 1
                [vrX_clu1(end+1), vrY_clu1(end+1)] = deal(x1, rho1);
                set(hClu1, 'XData', vrX_clu1, 'YData', vrY_clu1);
            case 2  
                [vrX_clu2(end+1), vrY_clu2(end+1)] = deal(x1, rho1);
                set(hClu2, 'XData', vrX_clu2, 'YData', vrY_clu2);
            otherwise
                continue;
        end
        viClu(iSpk1) = iClu1_parent;
        mrX_edge(:,end+1) = [vrX(iSpk1), vrX(iSpk1_parent), nan];
        mrY_edge(:,end+1) = [vrRho(iSpk1), vrRho(iSpk1_parent), nan];
        set(hLink, 'XData', mrX_edge(:), 'YData', mrY_edge(:));
        vl_ = viClu > 0;
        accuracy1 = nanmean(viClu(vl_) == viClu_gt(vl_));
        hTitle.String = sprintf('Accuracy: %0.2f (%0.0f%% done)', accuracy1, mean(vl_)*100);
    end
    pause(vrT_wait(iiSpk));
    drawnow limitrate;
end
delete(hLine);
% 96.5% accuracy
% compare with density dip

%% 2D demo