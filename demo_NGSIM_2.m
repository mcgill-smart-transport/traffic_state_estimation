clear,clc

load ngsim_us101_lane2_deltat_5.0_deltad_10.0.mat
load ngsim_us101_lane2_traj.mat

% colormap
cm_jet= flipud(jet);
cm = flipud(jet);
cm_jet(1,:) = 1;            % speed 0 = white

% cut entire speed at beginning/end 
V(isnan(V)) = 0;
V = V(11:140,21:end-40);

% complete ground truth 
while find(V==0)
    zeroid = find(V==0);
    V(zeroid) = V(zeroid-1);
end

clearvars K Q zeroid            % delete unused variables


%% 
iter = 1;
res = [];


for seed = 3000                  % random seed
    
delta = 0.02;                    

[veh,q] = genData(raw_data, delta, seed); % generate training data

mr = 1 - sum(q(:))/(size(V,1)*size(V,2)); % missing rate

imagesc(veh)
%%

% tau_s = 30:10:50;      % spatial embedding length set
% tau_t = 70:10:90;      % temporal embedding length set

tau_s = 50;
tau_t = 110;

order_flag = 0;         % order_flag=1: balanced matrix;  0: [1 4 2 3]

    for itau = 1:length(tau_s)
        for jtau = 1:length(tau_t)

            tic
            tau = [tau_s(itau), tau_t(jtau)];

            [N,T] = size(veh);
            sizeh = [tau N-tau(1)+1 T-tau(2)+1];    % original size of Hankel tensor

            % determine the order to build a balanced matrix
            if order_flag == 1
                [~, maxid] = max(sizeh);
                [~, minid] = min(sizeh);

                mid = [1,2,3,4];
                mid([maxid,minid]) = [];

                order = [minid maxid mid];
            else
                
                order = [1,4,2,3];
            end

            % halrtc parameters
            hal.rho = 1e-6;
            hal.max_rho = 1;
            hal.max_iter = 200;
            hal.beta = 1.1;
            hal.tol = 0.0001;
            hal.order = order;
            hal.plotf = 1;
            hal.sizeh = sizeh;
            hal.seed = seed;

            for theta = 10

                [mat_hat, rmse, mae] = STH_LRTC(veh, V, q, tau, theta, hal);

                res(iter,:) = [seed,tau,theta,toc,order_flag, rmse, mae];
                iter = iter + 1;
            end

        end

    end
end
