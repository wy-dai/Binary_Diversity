% satPlusDisk.m
% Assembling the disk and rubble-pile satellite.
% Output: the input granular info file for the new customized DEMBody code
% (with customized cohesion for each particle).
%
% Wen-Yue Dai, Tsinghua Univ, 2024
%--------------------------------------------------------------------------

clear
sat_filepath = 'satellite_pile\';
disk_filepath = 'disk_distr\';

sat_no = 'sat_001Mp';
disk_no = '004md_p3_m3';

M_p = 5.5e11;
M_s = 0.01;
M_d = 0.04;

a_sat = 3000;
G = 6.67384e-11;

%% Disk particle information
disk_file = [disk_filepath, disk_no, '\input_points.csv'];
disk = csvread(disk_file,1,0);

n_disk = length(disk(:,1));
disk = [disk zeros(n_disk,1)];

M_d = sum(disk(:,1));

%% satellite particle information 
sat_file = [sat_filepath, sat_no, '.csv'];
sat = csvread(sat_file,1,0);

sat = sat(:,1:(end-1));
m_sat = sat(:,1);
pos_sat = sat(:,3:5);
vel_sat = sat(:,6:8);
w_sat = sat(:,9:11);

Center = sum(m_sat.*pos_sat,1)./sum(m_sat); % center of mass (MC)
r_sat = sqrt(pos_sat(:,1).^2 + pos_sat(:,2).^2 + pos_sat(:,3).^2); % r for each particle

Vel = sqrt(G*(M_p+M_d)/a_sat);
W = Vel/a_sat;

pos_sat = pos_sat - Center;
w_sat = [0 0 W].*ones(length(m_sat),3);
vel_sat = cross(w_sat, pos_sat);

vel_sat = vel_sat + [0 Vel 0];
pos_sat = pos_sat + [a_sat 0 0];

sat(:,3:5) = pos_sat;
sat(:,6:8) = vel_sat;
sat(:,9:11) = w_sat;

% put a central particle into the bottom of the list (to mark the satellite)
[~,I] = min(r_sat);
n = length(m_sat);
sat([I,n],:) = sat([n,I],:);

%% Writing the intintial condition to input_points.txt
A = [disk; sat];
% A = disk;

fileID = fopen('input_points.txt','w');
fprintf(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f\n',A');
fclose(fileID);

fileID = fopen('input_points.csv','w');
fprintf(fileID,'m,i,x,y,z,vx,vy,vz,wx,wy,wz,r,c\n');
fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',A');
fclose(fileID);