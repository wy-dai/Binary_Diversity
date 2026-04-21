%--------------------------------------------------------------------------
% clusterAnalysis_sat.m
%--------------------------------------------------------------------------
% Calculating the mass and orbital elements of clusters
% INPUT: the output of module clusterInfoFinder_no_sat.m
% 
% Wen-Yue Dai, Tsinghua Univ, 2024
%--------------------------------------------------------------------------
clear
% Parameters and File Infomation
%--------------------------------------------------------------------------
M_p = 5.5e11;
M_d = 0.06;
G = 6.67384e-11;
%
load('disk_006md_40d_2.mat');
%--------------------------------------------------------------------------
%
%% read cluster information
for i = 1:length(SERIES)
    cluNum(i) = length(SERIES(i).CLUSTER.Size);
    cluMass{i} = SERIES(i).CLUSTER.Mass;
    cluAxis{i} = SERIES(i).CLUSTER.Axis;
    cluPos{i} = SERIES(i).CLUSTER.Pos;
    cluVel{i} = SERIES(i).CLUSTER.Vel;
    cluW{i} = SERIES(i).CLUSTER.W;
    if isempty(cluMass{i})
        maxMass(i) = 0;
    else
        maxMass(i) = cluMass{i}(1);
    end
    time(i) = SERIES(i).time;
end
time = time./(24*3600);
%% plot cluster information
figure
subplot(2,1,1)
plot(time,maxMass./M_p)
ylabel('max aggregate mass')
xlabel('t / day')

subplot(2,1,2)
plot(time,cluNum)
ylabel('cluster number')
xlabel('t / day')

%% cluster orbital elements calculateing
%
a = cell([length(SERIES),1]);
i = cell([length(SERIES),1]);
e = cell([length(SERIES),1]);
Omega = cell([length(SERIES),1]);
omega = cell([length(SERIES),1]);
f = cell([length(SERIES),1]);
for j = 1:length(SERIES)
    if isempty(cluMass{j})
        a{j} = 0;
        i{j} = 0;
        e{j} = 0;
        Omega{j} = 0;
        omega{j} = 0;
        f{j} = 0;
    else
        a{j} = [];
        i{j} = [];
        e{j} = [];
        Omega{j} = [];
        omega{j} = [];
        f{j} = [];
        for k = 1:length(cluMass{j})
            pos = cluPos{j}(k,:);
            vel = cluVel{j}(k,:);
            h = cross(pos,vel);
            i{j}(k) = acos(h(3)/norm(h));
            %
            r = norm(pos);
            a{j}(k) = 1/(2/r-vel*vel'/(G*M_p));
            %
            E = cross(vel, h)./(G*M_p) - pos./r;
            e{j}(k) = norm(E);
            %
            N = cross([0 0 1], h);
            Om = N(1)./norm(N);
            if N(2)>0
                Omega{j}(k) = acos(Om);
            else
                Omega{j}(k) = -acos(Om);
            end
            %
            om = N*E'./(norm(N)*norm(E));
            if E(3)>0
                omega{j}(k) = acos(om);
            else
                omega{j}(k) = -acos(om);
            end
            %
            F = E*pos'./(norm(E)*norm(pos));
            if pos*vel'>0
                f{j}(k) = acos(F);
            else
                f{j}(k) = -acos(F);
            end
        end
    end
end
%% max cluster OEle
for j = 1:length(SERIES)
    maxClu_a(j) = a{j}(1);
    maxClu_i(j) = i{j}(1);
    maxClu_e(j) = e{j}(1);
end

%% clusters OEle
% showing the largest n clusters' orbit cluster
n = 3;
a_n = zeros(length(SERIES),n);
i_n = zeros(length(SERIES),n);
e_n = zeros(length(SERIES),n);
Omega_n = zeros(length(SERIES),n);
omega_n = zeros(length(SERIES),n);
f_n = zeros(length(SERIES),n);
m_n = zeros(length(SERIES),n);
for j = 1:length(SERIES)
    if cluNum(j)>=n
        a_n(j,:) = a{j}(1:n);
        i_n(j,:) = i{j}(1:n);
        e_n(j,:) = e{j}(1:n);
        Omega_n(j,:) = Omega{j}(1:n);
        omega_n(j,:) = omega{j}(1:n);
        f_n(j,:) = f{j}(1:n);
        m_n(j,:) = cluMass{j}(1:n);
    else
        cN = cluNum(j);
        if cN>0
            a_n(j,1:cN) = a{j};
            i_n(j,1:cN) = i{j};
            e_n(j,1:cN) = e{j};
            Omega_n(j,1:cN) = Omega{j};
            omega_n(j,1:cN) = omega{j};
            f_n(j,1:cN) = f{j};
            m_n(j,1:cN) = cluMass{j};
        else
            continue;
        end
    end
end

%% clusters and satellite OEle display
figure
subplot(3,1,1)
plot(time,m_n(:,1)./M_p,'red','LineWidth',1.0)
ylabel('m / M_p')
xlabel('t (day)')
hold on
plot(time,m_n(:,2)./M_p,'blue','LineWidth',1.0)
plot(time,m_n(:,3)./M_p,'green','LineWidth',1.0)
legend('1st largest aggr.', '2nd largest aggr.', '3rd largest aggr.')

% subplot(4,1,2)
% plot(time,i_n(:,1).*(180/pi),'red','LineWidth',1.0)
% ylabel('incline')
% xlabel('t (day)')
% hold on
% plot(time,i_n(:,2).*(180/pi),'blue','LineWidth',1.0)
% plot(time,i_n(:,3).*(180/pi),'green','LineWidth',1.0)
% legend('1st largest aggr.', '2nd largest aggr.', '3rd largest aggr.')

subplot(3,1,2)
plot(time,a_n(:,1),'red','LineWidth',1.0)
ylabel('semi-major axis')
xlabel('t (day)')
hold on
plot(time,a_n(:,2),'blue','LineWidth',1.0)
plot(time,a_n(:,3),'green','LineWidth',1.0)
legend('1st largest aggr.', '2nd largest aggr.', '3rd largest aggr.')

subplot(3,1,3)
plot(time,e_n(:,1),'red','LineWidth',1.0)
ylabel('eccentricity')
xlabel('t (day)')
hold on
plot(time,e_n(:,2),'blue','LineWidth',1.0)
plot(time,e_n(:,3),'green','LineWidth',1.0)
legend('1st largest aggr.', '2nd largest aggr.', '3rd largest aggr.')

