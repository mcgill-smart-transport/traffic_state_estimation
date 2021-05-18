clear,clc
close all

load syn

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


V = truth(end:-1:1,:);

% convert m/s to ft/s
coef = 3.28084;


V = V * coef;
veh = veh * coef;


q = logical(veh);           % mask matrix
mr = 1-sum(q,'all')/(size(veh,1)*size(veh,2));  % missing rate
Trainrmse = sqrt(norm(veh(q(:)) - V(q(:)), 'fro' )^2/(sum(q(:))));

iter = 1;
res = [];
final_res = [];


[N,T] = size(veh);

deltax = 15*coef;                     % spatial resolution 
deltat = 10;                     % temporal resolution 

idx = deltax * repmat(1:N,T,1)';      % location 
idt = deltat * repmat(1:T,N,1);       % time 


trainId = find(q==1);

vi = veh(trainId);
xi = idx(trainId);
ti = idt(trainId);


%%

for sigma = 295
     for tau = 6

    sumPhi_free = zeros(N,T);
    sumPhi_cong = zeros(N,T);
    phi_free = zeros(N,T);
    phi_cong = zeros(N,T);

    sumVf = zeros(N,T);
    sumVc = zeros(N,T);

    free = -20*coef;             
    cong = 3*coef;

    vthr = 15*coef;
    deltav = 4*coef;
 

    for i = 1:length(vi)
            for x = 1:N                 
                for t = 1:T
                    phi_free(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs((t*deltat - ti(i) - (x*deltax - xi(i))/free))/tau));
                    phi_cong(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs((t*deltat - ti(i) - (x*deltax - xi(i))/cong))/tau));
                end
            end

        
            sumPhi_free =  sumPhi_free + phi_free;
            sumPhi_cong = sumPhi_cong + phi_cong;               % eq(2)
            sumVf = sumVf + phi_free * vi(i);
            sumVc = sumVc + phi_cong * vi(i);
    end
    

    Vfree = sumVf./sumPhi_free;         % eq(4)
    Vcong = sumVc./sumPhi_cong;         % eq(5)


    w = (0.5 * (1 + tanh( (vthr - min(Vfree, Vcong))/deltav)));    % eq(7)


    vhat = w.*Vcong + (1-w).*Vfree;     % eq(6)

    diff = V.*~q - vhat.*~q;
    ntest = sum(~q,'all');
    rmse = sqrt(norm(diff, 'fro' )^2/ntest);
    mae = sum(abs(diff),'all')/ntest;

    res = [res; sigma tau rmse mae];
    
    figure
    subplot(311)
    imagesc(V)
    title('V')
    subplot(312)
    imagesc(vhat)
    colormap(cm)
    title('vhat')
    subplot(313)
    imagesc(veh)
    colormap(cm)
    title('trajectory')
    
    disp(['sigma ' num2str(sigma) ', tau=' num2str(tau) ', rmse_test=' num2str(rmse)]);
    

    end
end

[minRMSE, minID] = min(res(:,3));
minMAE = res(minID,4);

final_res = [final_res; mr Trainrmse minRMSE minMAE];



