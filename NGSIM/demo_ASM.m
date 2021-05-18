
clear,clc

load NGSIM                  % seed = 3000

% colormap 
cm_jet= flipud(jet);
cm = flipud(jet);
cm_jet(1,:) = 1;            % speed 0 = white

set(0, 'DefaultFigureColormap', cm)
set(0, 'DefaultFigureColormap', cm_jet)

%% Main

iter = 1;
res = [];
mr = 1 - sum(q(:))/(130*480);
Trainrmse = sqrt(norm(veh.*q - V.*q, 'fro' )^2/(sum(q(:))));

[N,T] = size(veh);

deltax = 10;                    % spatial resolution feet
deltat = 5;                     % temporal resolution second

idx = deltax * repmat(1:N,T,1)';      % location
idt = deltat * repmat(1:T,N,1);       % time

trainId = find(q==1);

vi = veh(trainId);
xi = idx(trainId);
ti = idt(trainId);



% para
for sigma = 200
    for tau = 10
        tic
                
        sumPhi_free = zeros(N,T);
        sumPhi_cong = zeros(N,T);
        phi = zeros(N,T);
        phi_free = zeros(N,T);
        phi_cong = zeros(N,T);
        
        sumVf = zeros(N,T);
        sumVc = zeros(N,T);
        
        
        unit = 0.91;        % km/h to ft/s   
        free = -60*unit;
        cong = 10;
        
        
        for i = 1:length(vi)
            for x = 1:N
                for t = 1:T
                    phi_free(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs((t*deltat - ti(i) - (x*deltax - xi(i))/free))/tau));
                    phi_cong(x,t) = exp(-( abs(x*deltax - xi(i))/sigma + abs((t*deltat - ti(i) - (x*deltax - xi(i))/cong))/tau));
                end
            end
            
            
            sumPhi_free = sumPhi_free + phi_free;          % eq(2)
            sumPhi_cong = sumPhi_cong + phi_cong;
            sumVf = sumVf + phi_free * vi(i);
            sumVc = sumVc + phi_cong * vi(i);
        end
        
        
        Vfree = sumVf./sumPhi_free;         % eq(4)
        Vcong = sumVc./sumPhi_cong;         % eq(5)
        
        
    end
end


%%

deltav = 10*unit;
vthr = 40;
w = (0.5 * (1 + tanh( (vthr - min(Vfree, Vcong))/deltav)));    % eq(7)


vhat = w.*Vcong + (1-w).*Vfree;     % eq(6)

diff = V.*~q - vhat.*~q;
ntest = sum(~q,'all');
rmse = sqrt(norm(diff, 'fro' )^2/ntest);
mae = sum(abs(diff),'all')/ntest;

res = [res; sigma tau mr Trainrmse rmse mae];

subplot(3,1,1);
imagesc(vhat)
subplot(3,1,2);
imagesc(V)
subplot(3,1,3);
imagesc(veh);


disp(['sigma ' num2str(sigma) ', tau=' num2str(tau) ', rmse_test=' num2str(rmse)]);
toc

asm = vhat;


