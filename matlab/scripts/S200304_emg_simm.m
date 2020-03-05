% simuate emg time trace
% generate arbitrary waveforms and spatially weigh
% distribute nerves evenly and see how it moves

%% in 3d space distribute lines. bundle of muscle fibers. say 32 of them.

P = struct('r_s', .9, 'r_e', 1, 'n_e', 16, 'n_s', 16, 'r_s1', .7, ...
    'xy_ratio', 1, 'phi_e', 2*pi/32, 'n_r', 1000);

vrR_s = P.r_s * mod([1:P.n_s],2) + P.r_s1 * mod([0:P.n_s-1],2);
vrX_s = vrR_s .* cos((1:P.n_s)./P.n_s*(2*pi));
vrY_s = vrR_s .* sin((1:P.n_s)./P.n_s*(2*pi)) * P.xy_ratio;

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
close(v);

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
close(v);


%% action potential of the nerve, make it travel at different rate
% discretize the space

x = linspace(

