% simuate emg time trace
% generate arbitrary waveforms and spatially weigh
% distribute nerves evenly and see how it moves

%% in 3d space distribute lines. bundle of muscle fibers. say 32 of them.

P = struct('r_s', .8, 'r_e', 1, 'n_e', 16, 'n_s', 8, 'r_s1', .7, ...
    'xy_ratio', 1, 'phi_e', 2*pi/32, 'n_r', 1000);

%vrR_s = P.r_s * mod([1:P.n_s],2) + P.r_s1 * mod([0:P.n_s-1],2);
vrX_s = P.r_s .* cos((1:P.n_s)./P.n_s*(2*pi));
vrY_s = P.r_s .* sin((1:P.n_s)./P.n_s*(2*pi)) * P.xy_ratio;

get_pos_ = @(r, n, p)r .* [cos((1:n)'./n*(2*pi)+p), sin((1:n)'./n*(2*pi)+p)];
mrPos_e = get_pos_(P.r_e, P.n_e, P.phi_e);
mrPos_r = get_pos_(P.r_e, P.n_r, 0);

%% 2D geometry model
figure_ = @(x)figure('Color','w','Name',x,'NumberTitle','off');
figure_('electrode:red, nerve:black'); hold on;
title('2D model');
axis equal; grid on;

plot(vrX_s, vrY_s, 'ko');
plot(mrPos_r(:,1), mrPos_r(:,2), 'k:');
plot(mrPos_e(:,1), mrPos_e(:,2), 'ro');


%% 3D geometry model
figure_('electrode:red, nerve:black'); hold on;
title('3D model');

% plot arm
model = cylinderModel([-2,0,0,2,0,0,1]);
h=plot(model);
h.EdgeAlpha=0;
h.FaceColor='k';
h.FaceAlpha=.1;

% plot nerve
vrZ_s = [zeros(1, P.n_s), ones(1, P.n_s)];
for i_s = 1:P.n_s
    plot3([-2,2], repmat(vrX_s(i_s),[1,2]), repmat(vrY_s(i_s),[1,2]), 'k-');
end

% plot elec
for i_e = 1:P.n_e
    plot3([-.1,0,.1], repmat(mrPos_e(i_e,1),[1,3]), repmat(mrPos_e(i_e,2),[1,3]), 'r.');
end
view(3);


%% drift 2D movie

P_drift = struct('ang_mod', 2*pi/16/4, 'n_period', 2, 'nt', 1000);

figure_ = @(x)figure('Color','w','Name',x,'NumberTitle','off');
figure_('electrode:red, nerve:black'); hold on;
title('2D model')

plot(vrX_s, vrY_s, 'ko');
plot(mrPos_r(:,1), mrPos_r(:,2), 'k:');
axis equal; grid on; 
vrPhase_t = sin(linspace(0, 2*pi*P_drift.n_period, P_drift.nt)) * P_drift.ang_mod;
plot(mrPos_e(:,1), mrPos_e(:,2), 'mo');
hAx = gca;
h_e1 = [];
uiwait(msgbox('press OK to start movie'));
v = VideoWriter('elec_2d_drift.avi'); open(v);
for it=1:P_drift.nt
    mrPos_e1 = get_pos_(P.r_e, P.n_e, P.phi_e+vrPhase_t(it));
    if ~isempty(h_e1)
        delete(h_e1);
    end
    h_e1 = plot(hAx,mrPos_e1(:,1), mrPos_e1(:,2), 'ro');
%     drawnow limitrate;
    writeVideo(v,getframe(gcf));
end
close(v); disp(['wrote to: ', v.Filename]);

%% drift 3D movie

figure_('electrode:red, nerve:black'); hold on;
title('3D model');

% plot arm
model = cylinderModel([-2,0,0,2,0,0,1]);
h=plot(model);
h.EdgeAlpha=0;
h.FaceColor='k';
h.FaceAlpha=.1;

% plot nerve
vrZ_s = [zeros(1, P.n_s), ones(1, P.n_s)];
for i_s = 1:P.n_s
    plot3([-2,2], repmat(vrX_s(i_s),[1,2]), repmat(vrY_s(i_s),[1,2]), 'k-');
end

% plot elecclc
for i_e = 1:P.n_e
    plot3([-.1,0,.1], repmat(mrPos_e(i_e,1),[1,3]), repmat(mrPos_e(i_e,2),[1,3]), 'm.');
end
view(3);
hAx=gca;
% movie loop
h_e1 = [];
uiwait(msgbox('press OK to start movie'));
F(P_drift.nt) = struct('cdata',[],'colormap',[]);
v = VideoWriter('test1.avi'); open(v);
for it=1:P_drift.nt
    mrPos_e1 = get_pos_(P.r_e, P.n_e, P.phi_e+vrPhase_t(it));
    if ~isempty(h_e1)
        delete(h_e1);
    end
    vrX1 = repmat(mrPos_e1(:,1),[1,3]); 
    vrY1 = repmat(mrPos_e1(:,2),[1,3]);
    vrZ1 = repmat([-.1,0,.1], [P.n_e,1]);
    h_e1 = plot3(hAx, vrZ1(:), vrX1(:), vrY1(:), 'r.');
%     drawnow limitrate;
    writeVideo(v,getframe(gcf));
end
close(v); disp(['wrote to: ', v.Filename]);


%% action potential of the nerve, make it travel at different rate
% import waveform and discritize space, propagate waveform on a line
vcFile_template = '~/ceph/groundtruth/mearec_synth/tetrode/templates_tetrode.h5';
S_h5 = HDF2Struct(vcFile_template);
trWav_clu = S_h5.templates;
% start from zero
trWav_clu = trWav_clu - trWav_clu(1,:,:);
figure; plot(trWav_clu(:,:,1));

[vrMin_clu, viChan_min_clu] = min(min(trWav_clu,[],1),[],2); 

% single channel by selecting the largest
dimm_wav = size(trWav_clu);
mrWav_clu = reshape(trWav_clu, dimm_wav(1), []);
mrWav_clu = mrWav_clu(:, sub2ind(dimm_wav(2:3), viChan_min_clu(:)', 1:dimm_wav(3)));

% sort the neurons by size
[~, viClu_sort] = sort(vrMin_clu, 'ascend');

%% pick eight waveforms
viClu_s = 1:10:80;
mrWav_s = mrWav_clu(:, viClu_sort(viClu_s)); 
figure; plot(mrWav_s); title('Selected waveforms')

%% plot travelling wave, speed of pix/step
dx = 1/8;
P_wav = struct('nt', 1000, 'nx', 1000, 'x1_s', -10, 'x2_s', 10, 'x_p', dx/2, 'x_n', -dx/2, 'dy', .1, 'icell', 1);

% figure;

% pick wave and balance sum by boosting positive phase
vrV_s1 = mrWav_s(:, P_wav.icell) .* hanning(dimm_wav(1));
vl1 = vrV_s1>0; vrV_s1(vl1) = vrV_s1(vl1) * abs(sum(vrV_s1(~vl1)) / sum(vrV_s1(vl1)));
vrV_s1 = vrV_s1 / max(abs(vrV_s1));

vrVx = zeros(P_wav.nx, 1, 'single');
vrVx(1:dimm_wav(1)) = flipud(vrV_s1);

[vrVp, vrVn, vrVe] = deal(nan(P_wav.nt, 1, 'single'));
vrX_s = linspace(P_wav.x1_s, P_wav.x2_s, P_wav.nx)';
vrR_p = sqrt((P_wav.x_p - vrX_s).^2 + P_wav.dy^2);
vrR_n = sqrt((P_wav.x_n - vrX_s).^2 + P_wav.dy^2);
mrVxt_s = zeros(P_wav.nx, P_wav.nt, 'single');
%calc_Ve_ = @(x,r)sum(x./r);
% precompute
calc_Ve_ = @(x,r)sum(x.*log(1./r));
for it = 1:P_wav.nt
    %it1 = mod((it-1), size(mrWav_s,1))+1;
    it1 = min(it, dimm_wav(1));
    vrVx(2:end) = vrVx(1:end-1);    
    mrVxt_s(:,it) = vrVx;
    vrVp(it) = calc_Ve_(vrVx, vrR_p);
    vrVn(it) = calc_Ve_(vrVx, vrR_n);
    vrVe(it) = vrVp(it) - vrVn(it);
end
%%
figure_('nerve conduction');  uiwait(msgbox('press OK to start movie'));
ax=[];
ax(1) = subplot(3,1,1); h_Vx = plot(vrVx,'k'); hold on; plot(500+dx/2*100,0,'ro'); plot(500-dx/2*100,0,'bo'); xlabel('x pos'); grid on; hold on; ylabel('V');
ax(2) = subplot(3,1,2); hold on; h_Vp = plot(vrVp,'r'); h_Vn = plot(vrVn,'b'); xlabel('time'); ylabel('V (norm)'); grid on;
ax(3) = subplot(3,1,3); h_Vd = plot(vrVp-vrVn,'g'); xlabel('time'); ylabel('Vp-Vn (norm)'); grid on;
arrayfun(@(x)axis(x, [200 800 -1 1]), ax);

v = VideoWriter('travel.avi'); open(v);
for it=1:P_wav.nt
    vrVx = mrVxt_s(:,it);
    h_Vx.YData = vrVx;
    h_Vp.YData = vrVp(1:it)/max(abs(vrVp));
    h_Vn.YData = vrVn(1:it)/max(abs(vrVn));
    h_Vd.YData = vrVe(1:it)/max(abs(vrVe));
    writeVideo(v,getframe(gcf));
end
close(v); disp(['wrote to: ', v.Filename]);
%%
% plot narrowed waveform
norm_ = @(x)x/max(abs(x));
% figure; imagesc(mrVxt_s); xlabel('time'); ylabel('space');
figure_(''); hold on; plot(norm_(vrVp),'r'); plot(norm_(vrVn),'b'); plot(norm_(vrVe),'g'); ylabel('V(t)'); xlabel('time');
width_ = @(vl)abs(find(vl,1,'last') - find(vl,1,'first'));
fwhm_ = @(x)width_(abs(x)>max(abs(x))/2);
title(sprintf('dx=%0.3f, FWHM: %d, amp: %0.3f', dx, fwhm_(vrVe), min(vrVe)/min(vrVp)));
xlim([0 700]); grid on;

%% plot overlapping spikes
n_overlap = 70;
plot_ = @(x,c,n)plot(x, c,'LineWidth', n);
[vrV1, vrV2, vrV3] = deal(norm_(vrVe));
vrV2(n_overlap:end) = vrV2(1:end-n_overlap+1);
vrV3(1:end-n_overlap+1) = vrV3(n_overlap:end);
figure_(sprintf('dx=%0.2f',dx)); hold on; 
plot_(vrV2+vrV3,'k-',2); %plot_(vrV2,'r-',1); plot_(vrV3, 'b-',1);  
xlim([0 700]); grid on;

