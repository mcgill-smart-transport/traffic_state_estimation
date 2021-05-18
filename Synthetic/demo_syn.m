clc,clear

load syn                    % synthetic data 

cm_jet= flipud(jet);
cm = flipud(jet);
cm_jet(1,:) = 1;            % speed 0 = white
set(0, 'DefaultFigureColormap', cm)



veh = floating2(end:-1:1,:);
veh(isnan(veh)) = 0;

if find(truth==0)
    id = find(truth==0);
    truth(id) = (truth(id-1)+truth(id+1))/2;
end

q = logical(veh);           % mask matrix
mr = 1-sum(q,'all')/(size(veh,1)*size(veh,2));  % missing rate


% convert m/s to ft/s
coef = 3.28084;
% coef = 1;
truth = truth(end:-1:1,:) * coef;
veh = veh * coef;
Trainrmse = sqrt(norm(veh(q(:)) - truth(q(:)), 'fro' )^2/(sum(q(:))));    
    

iter = 1;
res = [];

%% STH-LRTC

tau_s = [2 5 10 20];
tau_t = [5 10 20 40];

tic
temp = [];
[N,T] = size(veh);

for itau = 1:length(tau_s)
    for jtau = 1:length(tau_t)
    
        tau = [tau_s(itau), tau_t(jtau)];
        sizeh = [tau N-tau(1)+1 T-tau(2)+1];    % original size of Hankel tensor
       


    % STH-LRTC parameters
    hal.rho = 1e-6/coef;
    hal.max_rho = 100;
    hal.max_iter = 200;
    hal.beta = 1.1;
    hal.tol = 0.001;
    hal.plotf = 1;
    hal.sizeh = sizeh;
    hal.theta = 3;
       
    [mat_hat, rmse, mae] = STH_LRTC(veh, truth, q, tau, hal);

    res(iter,:) = [tau,hal.theta,rmse, mae];
    iter = iter + 1;

    end
end

[minRMSE, minID] = min(res(:,5));
minMAE = res(minID,6);
final_res = [final_res; mr Trainrmse minRMSE minMAE];


