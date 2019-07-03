%% density clustering

P = struct('mu1', 0, 'sd1', 1, 'mu2', 4, 'sd2', 1, 'n1', 1000, 'n2', 1000, 'knn', 30, 'fontSize', 16, 'amp_drift', 1, 'nTime', 10);
rng(0);
vrY1 = linspace(0,2*pi,P.n1)';
vrY2 = linspace(0,2*pi,P.n2)';
vrX1 = randn([P.n1,1])*P.sd1 + P.mu1 + sin(vrY1) * P.amp_drift;
vrX2 = randn([P.n2,1])*P.sd2 + P.mu2 + sin(vrY2) * P.amp_drift;
figure_ = @()figure('Unit','Normalized','Color','w','OuterPosition',[0 0 .5 .5],'Toolbar','none','Menubar','figure');
figure1_ = @()figure('Unit','Normalized','Color','w','OuterPosition',[0 0 .5 1],'Toolbar','none','Menubar','figure');
xlabel_ = @(x)xlabel(x,'FontSize',P.fontSize, 'Interpreter', 'tex');
ylabel_ = @(x)ylabel(x,'FontSize',P.fontSize, 'Interpreter', 'tex');
zlabel_ = @(x)zlabel(x,'FontSize',P.fontSize, 'Interpreter', 'tex');

%% assign time using poisson distribution
figure_();  
plot(vrX1, vrY1, 'b.', vrX2, vrY2, 'r.');
xlabel_('x');
ylabel_('y');
grid on; set(gca,'YLim', [0,2*pi]);

%% combine two time series
vrY = [vrY1(:); vrY2(:)];
vrX = [vrX1(:); vrX2(:)];
viClu_gt = [repmat(2,P.n1,1); repmat(1,P.n2,1)];

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

%% slice by time, 10 slices



%% pairwise density plot
speed_plot = 30;
vrT_wait = logspace(log10(.01), log10(.25), numel(vrX)) / speed_plot;

% calc Rho
mrXY = [vrX, vrY];
[mrD_knn, miD_knn] = pdist2(mrXY, mrXY, 'euclidean', 'Smallest', P.knn);
vrRho = 1./mrD_knn(end,:)';

% calc Delta
mrD = pdist2(mrXY, mrXY);
mrD(vrRho <= vrRho') = nan;
[vrDelta, viDelta] = min(mrD);
vrDelta = vrDelta' .* vrRho;
[~, viDelta_srt] = sort(vrDelta, 'descend');
[imax1, imax2] = deal(viDelta_srt(1), viDelta_srt(2));
vrDelta(imax1) = max(vrDelta)*1.1;

hFig=figure1_(); 
hAx1=subplot(2,1,1); hold on; grid on; 
plot(vrRho, vrDelta, 'k.'); 
plot(vrRho(imax1), vrDelta(imax1), 'r*');
plot(vrRho(imax2), vrDelta(imax2), 'b*');
hLine1 = plot(vrRho(imax1)*[1, 1], get(gca,'YLim'), 'k-');
grid on;
xlabel_('\rho'); ylabel_('\delta');

vi_nn = viDelta(:);
mrX_nn = [vrX(:), vrX(vi_nn)];
mrY_nn = [vrRho(:), vrRho(vi_nn)];
mrX_nn(:,end+1) = nan;
mrY_nn(:,end+1) = nan;

hAx2=subplot(2,1,2); hold on; grid on; 
plot(vrX, vrY, 'k.');
plot(vrX(imax1), vrY(imax1), 'b*');
plot(vrX(imax2), vrY(imax2), 'r*');
xlabel_('x');
ylabel_('y');
set(gca,'YLim', [0,2*pi]);

pause;

% animate density line
[~, viRho_srt] = sort(vrRho, 'descend');
ch1_clu = {plot(hAx1,nan,nan,'ro','MarkerSize',3), ...
    plot(hAx1,nan,nan,'bo','MarkerSize',3)};
ch2_clu = {plot(hAx2,nan,nan,'ro','MarkerSize',3), ...
    plot(hAx2,nan,nan,'bo','MarkerSize',3)};
ch1_edge_clu = {plot(hAx1,nan,nan, '-', 'Color', [.5,0,0]), ...
        plot(hAx1,nan,nan, '-', 'Color', [0,0,.5])};
ch2_edge_clu = {plot(hAx2,nan,nan, '-', 'Color', [.5,0,0]), ...
        plot(hAx2,nan,nan, '-', 'Color', [0,0,.5])};
    
viClu = zeros(size(vrX));
viClu(imax1) = 2;
viClu(imax2) = 1;
[cmrXY1_clu, cmrXY1_edge_clu, cmrXY2_clu, cmrXY2_edge_clu] = deal({[],[]});
hTitle = title('Accuracy: ', 'FontSize', P.fontSize);
vrT_wait = vrT_wait(end:-1:1);
for iiSpk = 2:numel(vrX)
    iSpk1 = viRho_srt(iiSpk);
    [rho1, delta1] = deal(vrRho(iSpk1), vrDelta(iSpk1));
    [x1, y1] = deal(vrX(iSpk1), vrY(iSpk1));
    iSpk1_parent = viDelta(iSpk1);
    if viClu(iSpk1) == 0
        iClu_ = viClu(iSpk1_parent);
        
        % update hAx1 
        hLine1.XData = [rho1, rho1];
        mr1_ = [cmrXY1_clu{iClu_}; rho1, delta1];
        cmrXY1_clu{iClu_} = mr1_;
        set(ch1_clu{iClu_}, 'XData', mr1_(:,1), 'YData', mr1_(:,2));
        
        mr_edge1_ = [cmrXY1_edge_clu{iClu_};
            vrRho(iSpk1), vrDelta(iSpk1);
            vrRho(iSpk1_parent), vrDelta(iSpk1_parent);
            nan, nan];
        cmrXY1_edge_clu{iClu_} = mr_edge1_;
        set(ch1_edge_clu{iClu_}, 'XData', mr_edge1_(:,1), 'YData', mr_edge1_(:,2));

        % update hAx2
        mr2_ = [cmrXY2_clu{iClu_}; x1, y1];
        cmrXY2_clu{iClu_} = mr2_;
        set(ch2_clu{iClu_}, 'XData', mr2_(:,1), 'YData', mr2_(:,2));

        mr_edge2_ = [cmrXY2_edge_clu{iClu_};
            vrX(iSpk1), vrY(iSpk1);
            vrX(iSpk1_parent), vrY(iSpk1_parent);
            nan, nan];
        cmrXY2_edge_clu{iClu_} = mr_edge2_;
        set(ch2_edge_clu{iClu_}, 'XData', mr_edge2_(:,1), 'YData', mr_edge2_(:,2));

        viClu(iSpk1) = iClu_;
        vl_ = viClu > 0;
        accuracy1 = nanmean(viClu(vl_) == viClu_gt(vl_));
        hTitle.String = sprintf('Accuracy: %0.2f (%0.0f%% done)', accuracy1, mean(vl_)*100);
    end
    pause(vrT_wait(iiSpk));
    drawnow limitrate;
end
delete(hLine1);
% 96.5% accuracy
% compare with density dip

%% 2D demo