%--------------------------------------------------------------------------
% clusterAnalysis_sat.m
%--------------------------------------------------------------------------
% Calculating the mass and orbital elements of clusters and the satellite
% INPUT: the output of module clusterInfoFinder_pile_sat.m or
% clusterInfoFinder_sat.m
% 
% Wen-Yue Dai, Tsinghua Univ, 2024
%--------------------------------------------------------------------------
%
clear
%--------------------------------------------------------------------------
M_p = 5.5e11;
M_d = 0.04;
G = 6.67384e-11;
filepath = 'binary_004md_002sat_3000a_30d.mat';

%% read cluster information
%
M_c = M_p*(1+M_d);
%
load(filepath)
%
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

%% read satellite information
for i = 1:length(SERIES)
    satMass(i) = SERIES(i).SAT.Mass;
    satAxis(i,:) = SERIES(i).SAT.Axis;
end

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
            r = norm(pos,2);
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
%
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

%% satellite orbital elements
for j = 1:length(SERIES)
    satPos = SERIES(j).SAT.Pos;
    satVel = SERIES(j).SAT.Vel;
    sath = cross(satPos,satVel);
    sath2 = sath*sath';
    sati(j) = acos(sath(3)/sqrt(sath2));
    %
    satv2 = satVel*satVel';
    satr = norm(satPos,2);
    sata(j) = 1/(2/satr-satv2/(G*M_c));
    %
    satE = cross(satVel, sath)./(G*M_c) - satPos./satr;
    sate(j) = norm(satE);
    %
    satN = cross([0 0 1], sath);
    satOm = satN(1)./norm(satN);
    if satN(2)>0
        satOmega(j) = acos(satOm);
    else
        satOmega(j) = 2*pi-acos(satOm);
    end
    %
    satom = satN*satE'./(norm(satN)*norm(satE));
    if satE(3)>0
        satomega(j) = acos(satom);
    else
        satomega(j) = 2*pi-acos(satom);
    end
    %
    satF = satE*satPos'./(norm(satE)*norm(satPos));
    if satPos*satVel'>0
        satf(j) = acos(satF);
    else
        satf(j) = 2*pi-acos(satF);
    end
end

%% Output the final OEle for rebound simulation
INPUT1 = sprintf('%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e', m_n(end,1), cluAxis{end}(1,3), a_n(end,1), e_n(end,1), i_n(end,1), Omega_n(end,1), omega_n(end,1), f_n(end,1));
INPUT2 = sprintf('%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e', satMass(end), satAxis(end,3), sata(end), sate(end), sati(end), satOmega(end), satomega(end), satf(end));

%% new style plot
figure(1)
subplot(2,1,1)
plot(time,m_n(:,1)./M_p,'red','LineWidth',1.0)
ylabel('m / M_p')
xlabel('t (day)')
hold on
% plot(time,m_n(:,2)./M_p,'blue','LineWidth',1.0)
% plot(time,m_n(:,3)./M_p,'green','LineWidth',1.0)
plot(time,satMass./M_p,'black','LineWidth',1.0)
% legend('1st largest aggr.', '2nd largest aggr.', 'satellite')
legend('1st largest aggr.', 'satellite')

subplot(2,1,2)
aph = a_n(:,1) .* (1 + e_n(:,1))./380;
per = a_n(:,1) .* (1 - e_n(:,1))./380;
a = [per; aph(end:-1:1)];
t = [time, time(end:-1:1)];
p = fill(t,a,'red');
p.FaceAlpha = 0.2;      
p.EdgeColor = 'none';
hold on

plot(time,a_n(:,1)./380,'red','LineWidth',1.0)

aph = sata .* (1 + sate)./380;
per = sata .* (1 - sate)./380;
a = [per, aph(end:-1:1)];
t = [time, time(end:-1:1)];
p = fill(t,a,'blue');
p.FaceAlpha = 0.2;      
p.EdgeColor = 'none';

plot(time,sata./380,'black','LineWidth',1.0)

aFRL = 2.533;
plot(time,aFRL.*ones(size(time)),'--','Color','black','LineWidth',1.0)

xlabel('t / day')
