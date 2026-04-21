% Analysing the disk's dynamical evolution

clear
%
% File Information
%--------------------------------------------------------------------------
Sim_path = 'C:\Research\diskEvolution\Data\';
Sim_name = 'disk_006md_40d_2';

filepath = fullfile(Sim_path, Sim_name);
cluster_name = [Sim_name '.mat'];
cluster_path = fullfile(Sim_path,'result_clusterInfo/',cluster_name);
%
% Simulation Setting - Dynamics
%--------------------------------------------------------------------------
rho = 3500; % particle density
M_p = 5.5e11; % primary mass
R_p = 380; % primary radius
R_eq = 425; % primary equator radius
SAT = 0; % 0 for no sat, 1 for one sat
R_s = 100; % satellite radius
M_disk = 0.06; % disk mass (M_disk/M_p)
%--------------------------------------------------------------------------
aFRL = 2.533;
G = 6.67384e-11;
alpha = 1.0;
R_disk = alpha*R_p*aFRL;
%--------------------------------------------------------------------------
%
load(cluster_path)
fileno_ini = SERIES(1).fileno;
fileno_end = SERIES(end).fileno(end);
time_stp = SERIES(2).time - SERIES(1).time;

tic;
%
n = fileno_end-fileno_ini+1;
for k = 1:n
    % all the cluster information is recorded in the struct array
    % CLUSTER_SERIES
    fileno = fileno_ini-1+k;
    time(k) = time_stp*(fileno-1000);
end
%
r = 0:0.1:1; 
d = length(r);
a = R_eq + (aFRL*R_p-R_eq).*r;
S = zeros(1,d-1);
for j = 1:d-1
    S(j) = pi*(a(j+1)^2-a(j)^2);
end
%
rhoD = zeros(n,d-1);
MD = zeros(n,1);
L = zeros(n,3);
%
% MyPool = parpool;
for k = 1:n
    fileno = fileno_ini-1+k;
    demfile = [filepath,'\', num2str(fileno), 'X.csv'];
    M = csvread(demfile,1,0);
    %
    pos = M(:,1:3); % read all position
    vel = M(:,8:10); % read all velocity
    rs = M(:,7); %all radius
    np = length(rs);

    pos = pos(1:np-1,:);
    vel = vel(1:np-1,:);
    rs = rs(1:np-1);

    R = sqrt(pos(:,1).^2+pos(:,2).^2+pos(:,3).^2);
    %
    
    diskID = find(R<=R_disk);
    diskID = setdiff(diskID,SERIES(k).PRI.Index,"stable");
    pos = pos(diskID,:);
    vel = vel(diskID,:);
    rs = rs(diskID);
    R = R(diskID);
    %
    np = length(R);
    m = (4/3*pi*rho).*rs.^3;
    %
    % mass of the disk
    MD(k) = sum(m);
    %
    for i = 1:np
        for j = 1:(d-1)
            if R(i)<=a(j+1)
                rhoD(k,j) = rhoD(k,j) + m(i);
                break;
            end
        end
        L(k,:) = L(k,:) + m(i).*cross(pos(i,:),vel(i,:));
    end
    rhoD(k,:) = rhoD(k,:)./S;
end
% delete(MyPool);
toc

%% plot
figure
plot(time./(24*60*60), MD./(M_disk*M_p))
