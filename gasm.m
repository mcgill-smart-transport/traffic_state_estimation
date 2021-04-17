clear,clc
close all

load NGSIM
cm_jet= flipud(jet);
cm = flipud(jet);
cm_jet(1,:) = 1;            % speed 0 = white
set(0, 'DefaultFigureColormap', cm_jet)

%%

[N,T] = size(veh);

deltax = 10;                    % spatial resolution feet
deltat = 5;                     % temporal resolution second

idx = deltax * repmat(1:N,T,1)';      % location 
idt = deltat * repmat(1:T,N,1);       % time 

imagesc(veh)
%%

trainId = find(q==1);

vi = veh(trainId);
xi = idx(trainId);
ti = idt(trainId);


%%
res = [];
% para
for sigma = 10:5:50
     for tau = 10:5:50
         
    tic
    
    
    sumPhi = zeros(N,T);

    phi = zeros(N,T);
    phi_free = zeros(N,T);
    phi_cong = zeros(N,T);

    sumVf = zeros(N,T);
    sumVc = zeros(N,T);

    
    unit = 0.91;        % km/h to ft/s


    free = 70*unit;
    cong = -15*unit;

    vthr = 60*unit;
    deltav = 20*unit;


    for i = 1:length(vi)
            for x = 1:N                 
                for t = 1:T
                    phi(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs(t*deltat - ti(i))/tau));
                    phi_free(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs((t*deltat - ti(i) - (x*deltax - xi(i))/free))/tau));
                    phi_cong(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs((t*deltat - ti(i) - (x*deltax - xi(i))/cong))/tau));
                end
            end

        
            sumPhi = sumPhi + phi;          % eq(2)
            sumVf = sumVf + phi_free * vi(i);
            sumVc = sumVc + phi_cong * vi(i);
    end
    

    Vfree = sumVf./sumPhi;         % eq(4)
    Vcong = sumVc./sumPhi;         % eq(5)


    w = round(0.5 * (1 + tanh( (vthr - min(Vfree, Vcong))/deltav)));    % eq(7)


    vhat = w.*Vcong + (1-w).*Vfree;     % eq(6)

    diff = V.*~q - vhat.*~q;
    ntest = sum(~q,'all');
    rmse = sqrt(norm(diff, 'fro' )^2/ntest);
    mae = sum(abs(diff),'all')/ntest;

    res = [res; sigma tau rmse mae];
    figure
    imagesc(vhat)
    colormap(cm_jet)
    caxis([0 80])
    drawnow()
    
    disp(['sigma ' num2str(sigma) ', tau=' num2str(tau) ', rmse_test=' num2str(rmse)]);
    toc

    end
end

asm = vhat;
% [a b] = min(res(:,3))

% subplot(511)
% imagesc(Vfree)
% colorbar
% subplot(512)
% imagesc(Vcong)
% colorbar
% subplot(513)
% imagesc(vhat)
% colorbar
% subplot(514)
% imagesc(V)
% colorbar
% subplot(515)
% imagesc(veh)




