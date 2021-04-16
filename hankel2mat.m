function matrix = hankel2mat(tensor, sizeh)

dim = size(tensor);
order = [];

% obtain the original dim order: [stau, ttau, N, T]
if length(unique(dim)) < ndims(tensor)  % if exist same tau_s and tau_t
      for i = 2:ndims(tensor)  
            order = [order find(sizeh(i) == dim)];
      end
else
      for i = 1:ndims(tensor)  
            order = [order find(sizeh(i) == dim)];
      end
end


tensor = permute(tensor,order);

if ndims(tensor) == 3
    [tau, N, T] = size(tensor);
    
    temp = zeros(N+tau-1,T);
    
    for t = 1:T
        for i = 1:N
            idx = i:i+tau-1;
            temp(idx,t) = temp(idx,t) + tensor(:,i,t);
        end
    end
    

    mv = ones(1, N+tau-1)*tau;
    mv(1:tau) = 1:tau;
    mv(end:-1:end-tau+1) = 1:tau;
    
    mv = repmat(mv,[T,1])';

    matrix = temp./mv;


        
elseif ndims(tensor) == 4
    [stau, ttau, N, T] = size(tensor);
    
    temp = zeros(N+stau-1, T+ttau-1);
    
    for t = 1:T
        for i = 1:N
            idx1 = i:i+stau-1;
            idx2 = t:t+ttau-1;            
            temp(idx1,idx2) = temp(idx1, idx2) + (tensor(:,:,i,t));           
        end
        
    end
    
    mv = ones(1, N+stau-1)*stau;
    mv(1:stau) = 1:stau;
    mv(end:-1:end-stau+1) = 1:stau;
    
    mv = repmat(mv,[T+ttau-1,1])';
    
    
    mv2 = ones(1, T+ttau-1)*ttau;
    mv2(1:ttau) = 1:ttau;
    mv2(end:-1:end-ttau+1) = 1:ttau;
    
    mv2 = repmat(mv2,[N+stau-1,1]);
    
    matrix = temp./(mv.*mv2);
     
end




