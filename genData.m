%% generate training/testing data
function [veh,q] = genData(raw_data, mr, seed)


dt = 5;             % s
dx = 10;            % feet

id = unique(raw_data(:,1));

nid = length(id);
    
rng(seed);
trainId = randperm(nid,ceil(mr*nid));
veh = zeros(150,540);
nn = zeros(150,540);
q = zeros(150,540);


for i = trainId
    
    tr = raw_data(raw_data(:,1)==id(i),2:4);
    tr(:,1) = ceil(tr(:,1)/dt);
    tr(:,2) = max(ceil(tr(:,2)/dx),1);

    time = tr(:,1);
    space = tr(:,2);
    speed = tr(:,3);

    for k = 1:size(tr,1)
        veh(space(k),time(k)) = veh(space(k),time(k)) + speed(k);
        nn(space(k),time(k)) = nn(space(k),time(k)) + 1;
        q(space(k),time(k)) = 1;
    end
    
end
veh = veh./nn;
veh(isnan(veh)) = 0;


% final matrix (remove zeros): size(130*480)
veh = veh(11:140,21:end-40);
q = q(11:140,21:end-40);            % mask matrix
end