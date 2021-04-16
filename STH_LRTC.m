function [veh_hat, rmse_final, mae] = STH_LRTC(veh, truth, q, tau, theta, opt)
%% STH_LRTC: spatiotemporal Hankel tensor with low-rank tensor completion

% Based on Tensor completion for estimating missing values in visual data (HaLRTC) @Ji Liu

% Input:
% veh: spatiotemporal matrix with missing values (H = mat2hankel(veh, tau, order))
% truth: ground truth matrix (veh = V.*q)
% q: index matrix (1: training data position, 0: testing data position)
% tau: delay embedding lengths (tau = [tau_s, tau_t])
% theta:




% Output:
% veh_mat: estimated full spaiotemporal matrix




rho = opt.rho;
plotf = opt.plotf;
max_iter = opt.max_iter;
beta = opt.beta;
max_rho = opt.max_rho;
tol = opt.tol;
order = opt.order;
sizeh = opt.sizeh;


% matrix to Hankel tensor
obs_raw = mat2hankel(veh, tau, order); 
pos_miss = mat2hankel(~q, tau, order);

dim = size(obs_raw);
[vn, vt] = size(truth);


idx = ~q(:);        % testing index

normX = norm(obs_raw(:));

Y = obs_raw;
X = obs_raw;



for iter = 1 : max_iter
    
    % update M
    z = reshape(X+Y/rho,[dim(1)*dim(2),dim(3)*dim(4)]); % tensor to a balanced matrix
    trun_mat   = SVT_TNN(z, 1/rho,theta);   % shrinkage singular value decomposition
    W = reshape(trun_mat,dim);              % reshape balanced matrix to tensor
    
    
    M_mat_hat = hankel2mat(W,sizeh);        % add Hankel constraint
    M = mat2hankel(M_mat_hat,tau,order);
    
    
    % update X
    Xold = X;
    Xest = -Y/rho + M;
    
    X = obs_raw + Xest.* pos_miss;
    
    veh_hat = hankel2mat(X,sizeh);          % add Hankel constraint
    X =  mat2hankel(veh_hat,tau,order);
    
    
    % update Y
    Y = Y - rho * (M - X);
    

    % update rho
    rho = min (beta * rho, max_rho);
    
    % check convergence
    res(iter) = norm( X(:) - Xold(:) )/ normX ;
    
    
    
    % error on testing data   
    rmse(iter) = sqrt(norm(veh_hat(idx) - truth(idx), 'fro' )^2/(sum(idx)));
    
    if iter == 1 || mod(iter, 1) == 0
        disp(['tau= ' num2str(tau) ' theta= ' num2str(theta) ' iter= ' num2str(iter) ...
            ', rho=' num2str(rho) ', rmse_test=' num2str(rmse(iter)) ', res=' num2str(res(iter))]);
    end
    
   
    
    if plotf && mod(iter,1) == 0
        subplot(311)
        imagesc(truth)
        caxis([0 50])
        subplot(312)
        imagesc(veh)
        caxis([0 50])
        subplot(313)
        imagesc(veh_hat)
        caxis([0 50])
        drawnow()
    end
    
    if res(iter) < tol
        break
    end
   
end

rmse_final = rmse(iter);
mae = sum(abs(veh_hat(idx) - truth(idx)),'all')/sum(idx);

