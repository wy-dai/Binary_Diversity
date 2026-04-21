% Particle genaration for DEMBody Caculation
% ----------
% The model following Kokubo & Ida. 2000 and Hyodo et al. 2015.
% A debris disk surrounding a spherical asteroid (the primary)
% \Sigma ~ r^b (Power-law distr. for the disk density)
%
% Dai Wen-Yue @ Tsinghua U 2024

clear
N = 80000; % Maxmum number of particles
MAXITER = 900000; % max num iteration

Ifcohesion = true;

R_L = 3; % Lower bound of radius
R_U = 10; % Upper bound of radius
k = -3.0; % Power index for particle size disribution
b = -3.3; % Power index for the disk density distribution
rho = 3500; % particle density kg/m^3
G = 6.67384e-11; % gravitational constant
M_p = 5.5e11; % mass of the primary
M_d = 0.08; % expected disk mass (M_p)
R_p = 380; % Eff radius of the primary (m)
R_eq = 425; % Eff radius of the primary's equator (m)
aFRL = 2.533; % FRL (\rho_particle = 3500, \phi = 0.6)
% omega_p = 2*pi/(60*60*T_p); % spin rate of the primary

a_U = aFRL*R_p; % The outer edge id the disk
a_L = 0.46*aFRL*R_p; % inner edge
H = 0.4*R_p; % disk thickness

%% Granular length distribution - Power law
R = rand([N 1]); % radius
R = (R_L^k - (R_L^k - R_U^k) .* R).^(1/k);

%% Generating particles in a thin disk (surrounding the primary a_L<=r<=a_U)

c = 1; % Number of iteration
n = 1; % Counting the particles
X = []; % X = [x y z a vx vy wz]
B=b+2;
while c< MAXITER
    % tempr = a_L+(a_U-a_L)*rand();
    tempr = (a_L^B - (a_L^B - a_U^B) * rand())^(1/B);
    tempphi = 2*pi*rand();
    tempx = tempr*cos(tempphi);
    tempy = tempr*sin(tempphi);
    tempz = H*(1-2*rand())/2;

    if a_L+R(n)<=tempr && tempr<=a_U % in the disk or not
        if isempty(X)==1
            X(n,1) = tempx;
            X(n,2) = tempy;
            X(n,3) = tempz;
            X(n,4) = tempr;
            n = n+1;
        else
            Xtemp = X(1:n-1,:);
            for i = 1:n-1
                dist(i) = sqrt((tempx-X(i,1))^2+(tempy-X(i,2))^2+(tempz-X(i,3))^2);
                dist(i) = dist(i)-R(i)-R(n); % distance(particle n, particle i)
            end

            if isempty(find(dist<0,1))==1 % overlap or not
                X(n,1) = tempx;
                X(n,2) = tempy;
                X(n,3) = tempz;
                X(n,4) = tempr;
                n = n+1;
            end
        end
    end
    Md = (4/3*pi*rho).*R(1:n).^3;
    Md = sum(Md)/M_p;
    if Md>=M_d
        break
    end
    % if n>N
    %     break;
    % end
    c = c+1;
end
R = R(1:n-1);
M = 4./3.*rho.*pi.*R.^3; % particle mass

%% Velocity of particles (relative to the rotating primary)
for i = 1:n-1
    a = X(i,4);
    interiorID = find(X(:,4)<a); % find particles that r<r(i)
    interiorR = R(interiorID);
    interiorM = 4/3*rho*pi*sum(interiorR.^3);
    totalM = interiorM+M_p;
    omegaSq = G.*totalM./(a^3);
    omega = sqrt(omegaSq);
%     omega_rel = omega - omega_p; % relative rotation (in fixed frame on primary)
    
    X(i,5) = -X(i,2)*omega; % vx
    X(i,6) = X(i,1)*omega; % vy, and notice that vz=0
    X(i,7) = omega; % omega_z

%     % add norm noise
%     sigma = 0.01;
%     vnoise = normrnd(0,sigma)*a*omega;
%     wnoise = normrnd(0,sigma)*omega;
%     X(i,7) = X(i,7)+wnoise;
%     theta = 2*pi*rand();
%     X(i,5) = X(i,5)+cos(theta)*vnoise;
%     X(i,6) = X(i,6)+sin(theta)*vnoise;

end

if Ifcohesion
    Coh = zeros(size(R,1),1);
else
    Coh = [];
end

%% Writing the intintial condition to input_points.txt
fileID = fopen('input_points.txt','w');
A = [M 0.4.*M.*R.^2 X(:,1:3) X(:,5:7) R Coh];
fprintf(fileID,'%f %f %f %f %f %f %f 0 0 0 %f %f %f\n',A');
fclose(fileID);

fileID = fopen('input_points.csv','w');
fprintf(fileID,'m,i,x,y,z,vx,vy,vz,wx,wy,wz,r\n');
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,0,0,0,%f,%f,%f\n',A');

%%
r = 0:0.1:1; 
d = length(r);
a = a_L + (a_U-a_L).*r;
S = zeros(1,d-1);
for j = 1:d-1
    S(j) = pi*(a(j+1)^2-a(j)^2)*R_p^2;
end

rhoD = zeros(1,d-1);
for i = 1:length(M)
    for j = 1:(d-1)
        if X(i,4)<=a(j+1)
            rhoD(j) = rhoD(j) + M(i);
            break;
        end
    end
end
rhoD = rhoD./S;
plot(log(a(2:end)),log(rhoD)')
axis equal
hold on
p = polyfit(log(a(2:end)),log(rhoD(1:end))',1);