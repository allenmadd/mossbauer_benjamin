
%% Analysis 1 JAN 22 2019 LAB DAY 3
% Measure energy shifts for transitions. 
trans_energy = [2.6072e-07,1.5731e-07,5.3194e-08,-2.5638e-08,-1.2969e-07,-2.3339e-07,-2.5481e-07, -1.529e-07, -5.1092e-08, 2.4764e-08, 1.2712e-07, 2.309e-07]; 
trans_energy_error = [1.8672e-21,  3.7132e-23, 1.8493e-22, 4.3933e-23, 6.2515e-22,  3.8137e-22, 1.9614e-21, 1.2119e-21, 1.111e-21, 2.5069e-22, 8.7447e-22, 4.9659e-22];

figure; % just for fun
errorbar((1:6),trans_energy(1:6),trans_energy_error(1:6),'ro-'); grid on; hold on;
errorbar((7:12),trans_energy(7:12),trans_energy_error(7:12),'b*-'); 
legend('Left Side, First Bitch','Right Side, Second Bitch');

% look at difference between data sets. for understanding purposes
transition_diff = trans_energy(1:6) - trans_energy(12:-1:7); % difference between sets of data
figure;
plot((1:6),transition_diff,'o-'); grid on; 
title('Difference between left and right data sets');

%% Analysis 1: part b 
% deduce size of internal magnetic field felt by nucleus. 

% energy jumps for pairs that end at the same place, but start at different
% places. We can use this energy difference between these states to deduce
% H. We already know the ground state magnetic moment and I. 
db_energy = abs([trans_energy(11)-trans_energy(9), trans_energy(10)-trans_energy(8), trans_energy(5)-trans_energy(3), trans_energy(4)-trans_energy(2)]); % energy differences between m_j starting points
db_energy_error = [sqrt((trans_energy_error(11))^2+trans_energy_error(9)^2), sqrt((trans_energy_error(10))^2+(trans_energy_error(8))^2), sqrt((trans_energy_error(5))^2+(trans_energy_error(3))^2),sqrt((trans_energy_error(2))^2+(trans_energy_error(4))^2)]; % error in energy differences

% find H for pairs (2,4), (3,5), (9,11), (10,8)
m_0 = 2.93 *10^-13; % given ground state magnetic moment
I = 1/2; % given spin of the nucleus
dm_j = 1; % difference between m_j states for starting positions. same for each pair. 
H = I*(1/(m_0*dm_j))*db_energy; 
% H_error = -I*(1/(m_0*dm_j))*1.0e-20*[0.141,0.124,0.0652,0.00575]; % manually taking 3 sig figs into account yet because matlab doesn't have a function for that, apparently.

% resolution_H = H_error./H % this is super small, our H resolution is not reasonable. Use resolution set by velocity. look back in code for velocity resolution

H_error = H*resolution_vel;

figure; % just for visualization purposes
errorbar(1:4,H,H_error,'o'); grid on; xlabel('pair'); ylabel('H Magnetic field (Gauss)');
title('Magnetic Field per pair'); 

% H is 1.0e+05 * 3.0412    3.0318    3.1209    3.1220 Gauss

% H_error 1.0e+03 * -7.5935   -7.5702   -7.7926   -7.7953 


Hweight = 1./H_error.^2; 
H_avg = (sum(Hweight.*H)./sum(Hweight)); % average H from our answers

% H_avg = 2.0778e+05

%% Analysis 1 Part a
% Deduce the magnetic moment of the first excited state of Fe57
% first find energy differences between the following mj transitions in the
% excited state where I= 3/2. 
% mj -3/2 -1/2
diff5and6 = abs(trans_energy(12)-trans_energy(11)); 
err5and6 = sqrt(trans_energy_error(12)^2+trans_energy_error(11)^2);
diff11and12 = abs(trans_energy(6) - trans_energy(5)); 
err11and12 = sqrt(trans_energy_error(6)^2+trans_energy_error(5)^2);

% mj -3/2 1/2
diff6and4 = abs(trans_energy(10) - trans_energy(8));
err6and4 = sqrt(trans_energy_error(10)^2+trans_energy_error(8)^2);
diff10and8 = abs(trans_energy(6) - trans_energy(4)); 
err10and8 = sqrt(trans_energy_error(6)^2+trans_energy_error(4)^2);

% mj -1/2 1/2 
diff5and4 = abs(trans_energy(11) - trans_energy(10));
err5and4 = sqrt(trans_energy_error(11)^2+trans_energy_error(10)^2);
diff10and11 = abs(trans_energy(5) - trans_energy(4));
err10and11 = sqrt(trans_energy_error(5)^2+trans_energy_error(4)^2);

% put energy differences into vectors for left and right sides
ediff_mm_right = [diff5and4, diff6and4, diff5and6];
diff_mm_err_right = [err5and4, err6and4, err5and6];
ediff_mm_left = [diff10and11, diff10and8, diff11and12]; 
diff_mm_err_left = [err10and11, err10and8, err11and12]; 

% FINISH THIS TONIGHT ANG, SO HELP ME GOD

dm_jvec = [1,2,1]; % difference in m_j per transition
I = 3/2; % excited state nuclear spin quantum number
% calculate mu
mu_mm_right = -ediff_mm_right*I./(dm_jvec*H_avg*m_0); % mu in terms of mu/m_0
mu_mm_left = -ediff_mm_left*I./(dm_jvec*H_avg*m_0); 
mu_mm_err_right = -diff_mm_err_right*I./(dm_jvec*H_avg*m_0); 
mu_mm_err_left = -diff_mm_err_left*I./(dm_jvec*H_avg*m_0);

% god knows why I kept these left and right parts separate, because I end up combining
% them in the end
mu_mm = [mu_mm_left, mu_mm_right];
mu_mm_err = [mu_mm_err_left, mu_mm_err_right];

mu_weight = 1./mu_mm_err.^2; 
mu_avg = (sum(mu_weight.*mu_mm)./sum(mu_weight)); % average mu/m_0 

% AWESOME, we get mu_avg = -1.7090

%% ANALYSIS 1 ACTUAL VALUES
% H = 329400 Gauss
% m_1 = 1.71 * m_0 

%% woohoo!! So for each next data set I'm p sure we just need to plug in the data and choose the values that are with the peaks
