%--------------------------------------------------------------------------
% clusterInfoFinder_pile_sat.m
%--------------------------------------------------------------------------
% Reconstructing the cluster information from points file .csv.
% Customized for the simulation with a rubble pile satellite in orbit.
% A particle is set to be the representation of the satellite, it should be
% listed in the end in the input .csv file.
% Submodule findCluster_sat.m used.
%
% Wen-Yue Dai, Tsinghua Univ, 2024
%--------------------------------------------------------------------------
clear
%--------------------------------------------------------------------------
% File Information
%--------------------------------------------------------------------------
filepath_0 = 'Data\';
file_name = 'binary_004md_002sat_3500a_30d';

fileno_ini = 1000;
fileno_end = 1432;
%--------------------------------------------------------------------------
% Simulation Setting - Dynamics
%--------------------------------------------------------------------------
rho = 3500;     % particle density
M_p = 5.5e11;   % primary mass
%--------------------------------------------------------------------------
% Searching Details
%--------------------------------------------------------------------------
time_stp = 5000;    % timestep between files
% searching region boundary
Rlim = [0,5000];    % R
Zlim = [-500,500];  % Z
N_msh = 9;          % N_msh = mesh size / particle maximum radius
numPrs = 1e-2;      % numerical precesion of position (m)
N_agg_min = 20;     % minimum cluster size
%
%--------------------------------------------------------------------------
%
tic;
%
n = fileno_end-fileno_ini+1;
for k = 1:n
    % all the cluster information is recorded in the struct array
    % CLUSTER_SERIES
    fileno = fileno_ini-1+k;
    SERIES(k).fileno = fileno;
    SERIES(k).time = time_stp*(fileno-1000);
end
%
E = [1 0 0; 0 1 0; 0 0 1]; % unit matrix
%
filepath = [filepath_0, file_name,'\'];
MyPool = parpool;
parfor k = 1:n
    fileno = fileno_ini-1+k;
    demfile = [filepath, num2str(fileno), 'X.csv'];
    M = csvread(demfile,1,0);
    [CLUSTER_pri, CLUSTER_sat, CLUSTER] = findCluster_sat(M,Rlim,Zlim,N_msh,numPrs);
    %
    cluSize=[];
    cluMass=[];
    cluPos=[];
    cluVel=[];
    cluW=[];
    cluAM=[];
    a=[]; b=[]; c=[];
    satMass = [];
    satPos = [];
    satVel = [];
    satAM = [];
    satW = [];
    %
    N_p = length(M(:,1)); % number of particles (including pri and sat)
    %
    for i = 1:length(CLUSTER)
        clu_tmp = cell2mat(CLUSTER(i));
        mass_tmp = 4/3*pi*rho*M(clu_tmp,7).^3;
        % cluster mass and size
        cluMass(i) = sum(mass_tmp);
        cluSize(i) = length(mass_tmp);
    end
    % filter out large clusters (m>M_agg_min), and sort all clusters in
    % descend order of mass
    id_agg = find(cluSize>=N_agg_min);
    if ~isempty(id_agg)
        cluMass = cluMass(id_agg);
        CLUSTER = CLUSTER(id_agg);
        cluSize = cluSize(id_agg);
        [cluMass, id_agg] = sort(cluMass,'descend');
        CLUSTER = CLUSTER(id_agg);
        cluSize = cluSize(id_agg);
        %
        for i = 1:length(CLUSTER)
            clu_tmp = cell2mat(CLUSTER(i));
            mass_tmp= 4/3*pi*rho*M(clu_tmp,7).^3;
            % cluster size
            cluSize(i) = length(clu_tmp);
            % cluster mass
            cluMass(i) = sum(mass_tmp);
            % cluster mass center (MC) position
            cluPos(i,:) = mass_tmp'*M(clu_tmp,1:3)./cluMass(i);
            % cluster MC velocity
            cluVel(i,:) = mass_tmp'*M(clu_tmp,8:10)./cluMass(i);
            % cluster characteristic angular momentum
            cluAM(i,:) = mass_tmp'*cross(M(clu_tmp,1:3),M(clu_tmp,8:10),2);
            cluAM(i,:) = cluAM(i,:)-cluMass(i)*cross(cluPos(i,:),cluVel(i,:));
            % inertia tensor
            J = zeros(3);
            for j = 1:length(clu_tmp)
                r = M(clu_tmp(j),1:3)-cluPos(i,:);
                J_tmp = (r*r').*E-r'*r;
                J_tmp = mass_tmp(j).*J_tmp;
                J = J + J_tmp;
                
            end 
            J = J + (0.4*4/3*pi*rho*sum(M(clu_tmp,7).^5)).*E;
            % principle moments
            PM = sort(eig(J),'ascend'); % PM = [J1,J2,J3]
            % equivalent axises calculating
            sumPM = sum(PM);
            a(i) = sqrt(5*(0.5*sumPM-PM(1))/cluMass(i));
            b(i) = sqrt(5*(0.5*sumPM-PM(2))/cluMass(i));
            c(i) = sqrt(5*(0.5*sumPM-PM(3))/cluMass(i));
            % spin rate
            cluW(i,:) = cluAM(i,:)*inv(J)';
        end
        %
        SERIES(k).CLUSTER.Index = CLUSTER;
        SERIES(k).CLUSTER.Size = cluSize;
        SERIES(k).CLUSTER.Mass = cluMass;
        SERIES(k).CLUSTER.Size_min = N_agg_min;
        SERIES(k).CLUSTER.Pos = cluPos;
        SERIES(k).CLUSTER.Vel = cluVel;
        SERIES(k).CLUSTER.W = cluW;
        SERIES(k).CLUSTER.AM = cluAM;
        SERIES(k).CLUSTER.Axis = [a' b' c'];
    else
        SERIES(k).CLUSTER.Index = {};
        SERIES(k).CLUSTER.Size = [];
        SERIES(k).CLUSTER.Mass = [];
        SERIES(k).CLUSTER.Size_min = N_agg_min;
        SERIES(k).CLUSTER.Pos = [];
        SERIES(k).CLUSTER.Vel = [];
        SERIES(k).CLUSTER.W = [];
        SERIES(k).CLUSTER.AM = [];
        SERIES(k).CLUSTER.Axis = [];
    end
    %
    SERIES(k).PRI.Index = CLUSTER_pri;
    SERIES(k).SAT.Index = CLUSTER_sat;
    %
    % calculate dynamical properties of satellite
    idSat = find(CLUSTER_sat==N_p-1);
    mass_tmp = (4/3*pi*rho).*M(CLUSTER_sat,7).^3;
    % mass_tmp(idSat) = M_s;
    satMass = sum(mass_tmp);
    % satellite mass center (MC) position
    satPos = mass_tmp'*M(CLUSTER_sat,1:3)./satMass;
    % satellite MC velocity
    satVel = mass_tmp'*M(CLUSTER_sat,8:10)./satMass;
    % satellite characteristic angular momentum
    satAM = mass_tmp'*cross(M(CLUSTER_sat,1:3),M(CLUSTER_sat,8:10),2);
    satAM = satAM-satMass*cross(satPos,satVel);
    % satellite inertial tensor
    J = zeros(3);
    for j = 1:length(CLUSTER_sat)
        r = M(CLUSTER_sat(j),1:3)-satPos;
        J_tmp = (r*r').*E-r'*r;
        J_tmp = mass_tmp(j).*J_tmp;
        J = J + J_tmp;
    end
    J = J + (0.4*4/3*pi*rho*sum(M(CLUSTER_sat,7).^5)).*E;
    % satellite spin rate
    satW = satAM*inv(J);
    % principal moments and equivalent axises
    PM = sort(eig(J),'ascend');
    sumPM = sum(PM);
    a_sat = sqrt(5*(0.5*sumPM-PM(1))/satMass);
    b_sat = sqrt(5*(0.5*sumPM-PM(2))/satMass);
    c_sat = sqrt(5*(0.5*sumPM-PM(3))/satMass);
    %
    SERIES(k).SAT.Mass = satMass;
    SERIES(k).SAT.Pos = satPos;
    SERIES(k).SAT.Vel = satVel;
    SERIES(k).SAT.W = satW;
    SERIES(k).SAT.AM = satAM;
    SERIES(k).SAT.Axis = [a_sat b_sat c_sat];
    %
    % Mass of primary
    idPri = find(CLUSTER_pri==N_p);
    mass_tmp = (4/3*pi*rho).*M(CLUSTER_pri,7).^3;
    mass_tmp(idPri) = M_p;
    SERIES(k).PRI.Mass = sum(mass_tmp);
end
delete(MyPool);
toc

resultname = [filepath_0, 'result_clusterInfo\', file_name, '.mat'];
save(resultname, 'SERIES');