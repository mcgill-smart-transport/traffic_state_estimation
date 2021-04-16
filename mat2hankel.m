function [H] = mat2hankel(X, tau, order)
%% delay embedding (matrix to tensor) 
% original order: spatial delay * temporal delay * spatial left * temporal left
% X: matrix
% tau: delay embedding length
% order: reshape the tensor to get a balanced matrix later

[N,T] = size(X);

% for 3 dim
if length(tau) == 1              % only spatial delay
    H = zeros(tau,N-tau+1,T);    % incomplete tensor
    for t=1:T
        for k = 1:N-tau+1
            H(:,k,t) = X(k:k+tau-1,t);
        end
    end
    
       
else
    
    % for 4 dim
    
    stau = tau(1);
    ttau = tau(2);

    H = zeros(stau, ttau, N-stau+1,T-ttau+1);   % original hankel tensor


    for t=1:T-ttau+1
        for k = 1:N-tau+1
            H(:,:,k,t) = X(k:k+stau-1,t:t+ttau-1);
        end
    end

end


H = permute(H,order);





  
  
  
  
  
  
  
  
  