%% Low-rank Hankel Tensor Completion for Traffic State Estimation
% NGSIM dataset demo
% 2021-May, @Xudong Wang

clear,clc
close all

load NGSIM                  % seed = 3000

% colormap
cm_jet= flipud(jet);
cm = flipud(jet);
cm_jet(1,:) = 1;            % speed 0 = white


iter = 1;
res = [];
final = [];
mr = 1 - sum(q(:))/(size(V,1)*size(V,2)); % missing rate
seed = 3000;

%% Main

tau_s = 20:10:50;
tau_t = 20:10:50;

[N,T] = size(veh);

for itau = 1:length(tau_s)
    for jtau = 1:length(tau_t)

        tau = [tau_s(itau), tau_t(jtau)];
        sizeh = [tau N-tau(1)+1 T-tau(2)+1];    % original size of Hankel tensor

        % STH-LRTC parameters
        hal.rho = 1e-6;
        hal.max_rho = 1;
        hal.max_iter = 200;
        hal.beta = 1.1;
        hal.tol = 0.001;
        hal.plotf = 1;
        hal.sizeh = sizeh;
        hal.seed = seed;
        hal.theta = 6;

        tic
        [mat_hat, rmse, mae] = STH_LRTC(veh, V, q, tau, hal);
        toc


        res(iter,:) = [tau, hal.theta, toc, rmse, mae];
        iter = iter + 1; 
    end

end
temp = res(res(:,1)==seed,[6,7]);
[minRMSE, minID] = min(temp(:,1));

final = [final; mr minRMSE temp(minID,2)];


