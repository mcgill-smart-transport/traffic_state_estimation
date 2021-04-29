clc,clear
load syn.mat

cm_jet= flipud(jet);
cm = flipud(jet);
cm_jet(1,:) = 1;            % speed 0 = white
set(0, 'DefaultFigureColormap', cm)


veh = floating5;
veh(isnan(veh)) = 0;


subplot(211)
imagesc(truth)
colorbar
colormap(cm_jet)
caxis([0 40])
title('Ground truth')
subplot(212)
imagesc(veh)
colorbar
caxis([0 40])
title('5% trajectory')


q = logical(veh);           % mask matrix
mr = 1-sum(q,'all')/(size(veh,1)*size(veh,2));  % missing rate


%%

tau_s = [1,2,5,10,20,40];   % spatial embedding length set
tau_t = [5,10,20,40,60];    % temporal embedding length set


iter = 1;
res = [];
order_flag = 1;     % order_flag=1: balanced matrix; 0: [1 2 3 4]

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
        order = [1,2,3,4];
    end
    


    % halrtc parameters
    hal.rho = 1e-6;
    hal.max_rho = 1;
    hal.max_iter = 500;
    hal.beta = 1.1;
    hal.tol = 0.0001;
    hal.order = order;
    hal.plotf = 1;
    hal.sizeh = sizeh;

 
    
    for theta = 5
       
        [mat_hat, rmse, mae] = STH_LRTC(veh, truth, q, tau, theta, hal);
        
        res(iter,:) = [tau,theta,toc,order_flag, rmse, mae];
        iter = iter + 1;
    end

    end

end
