%Madeleine and Angela try their bests
%Mossbauer Analysis 17 Jan 2019


%
%Analysis Step 2
%

%files found at path: "C:\Users\Student\Documents\MadelAngela_mossbauer_15Jan2019";
clear all;
close all;
%upload the file as a table
filename="step15_15Jan2019.csv";
mossbauer_data= readtable(filename);

%graph that boi
mossbauer_data= removevars(mossbauer_data, {'Energy', 'Channel'});
stackedplot(mossbauer_data,'.-');
%may be necessary for the gaussian distribution
mossbauer_data= table2array(mossbauer_data);

%calculate error in data
filename= "step7_real_15Jan2019.csv";
background_counts = readtable(filename);
%1) with background counts

%transmission (no background counts) error = sqrt(counts)
transmission_error= sqrt(mossbauer_data);
figure; 
errorbar([1:1:1024],mossbauer_data,transmission_error);
%2) without background counts 
background_counts = mean(mossbauer_data(512:550,1)); %some random portion with no interesting troughs in it
background_counts_error = sqrt(background_counts);

absorption_error= sqrt(background_counts+mossbauer_data);
mossbauer_data_no_background =  background_counts- mossbauer_data;
figure;
errorbar([1:1:1024], mossbauer_data_no_background, absorption_error);

%converting channels to energy
dv = (10.22+9.59)/(500);
for i = 1:1:500
    v(i) = -9.59 + (i-1)*dv;
    E1(i)=14.4*(v(i)/(3*10^8));
    v2(i)=10.22 - (i-1)*dv;
    E2(i) = 14.4*(v2(i)/(3*10^8));
end

resolution_vel = max(dv./v(715:737))

%represent energies with opposite sloped velocity
old_channel= [1:1:1024];
%new_channel = zeros(1,1024);
energy = [E2(13:end), zeros(1,23), E1, E2(1:13)];
figure;
plot(old_channel, energy,'o-');
%
%Analysis Step 3
%
lorentzian = @(b,x) ( b(1)./( (x-b(2)).^2 + b(3) ) )
%
%notes from Jan 20th: it looks like we did every other peak in class? why?
%i will just do the even peaks and then if necessary i can do the others
%

%% PEAK ANALYSIS

%PEAK point5 (secretly peak 1)
%this good make it look like this
energypoint5 = energy(91:1:121);
peakpoint5= mossbauer_data_no_background(91:121,1);
weightpoint5 = 1./(absorption_error(91:121).^2);
%do this for all of the peaks this is just peak 1
cpoint5 = [1.8330*10^-15,2.604*10^-7,(10^-18)]; %changed name to "constant 1" aka "c1"
fitpoint5=fitnlm(energypoint5, peakpoint5, lorentzian,cpoint5,'Weights',weightpoint5)
figure; plot(energypoint5,peakpoint5)
hold on; plot(energypoint5,lorentzian([fitpoint5.Coefficients{1,1},fitpoint5.Coefficients{2,1},fitpoint5.Coefficients{3,1}],energypoint5))
figure; plot(fitpoint5.Residuals{1:end,1})
%line width squared: 4.6256e-17, error: 4.605e-18 
%delta E (transition energy): 2.6072e-07, error: 1.8672e-21 

%PEAK 1 (secretly peek a boo 2)
%this good make it look like this
energy1 = energy(140:1:190);
peak1= mossbauer_data_no_background(140:190,1);
weight1= 1./(absorption_error(140:190).^2);
%do this for all of the peaks this is just peak 1
c1 = [1600*10^(-18),1.55*10^(-7),(10^-18)]; %changed name to "constant 1" aka "c1"
fit1=fitnlm(energy1, peak1, lorentzian,c1, 'Weights', weight1)
figure; plot(energy1,peak1)
hold on; plot(energy1,lorentzian([fit1.Coefficients{1,1},fit1.Coefficients{2,1},fit1.Coefficients{3,1}],energy1))
figure; plot(fit1.Residuals{1:end,1})
%line width squared for peak 1 with  error  {3.43978799312567e-17,3.80066459794880e-18}
%transition energy: 1.5731e-07  , error:  3.7132e-23


%PEAK 1point5 (secretly peak 3)
%this good make it look like this
energy1point5 = energy(203:1:230);
peak1point5= mossbauer_data_no_background(203:230,1);
weight1point5=1./(absorption_error(203:230).^2);
%do this for all of the peaks this is just peak 1
c1point5 = [1.1580*10^-15,5.316*10^-8,(10^-18)]; 
fit1point5=fitnlm(energy1point5, peak1point5, lorentzian,c1point5, 'Weights', weight1point5);
figure; plot(energy1point5,peak1point5)
hold on; plot(energy1point5,lorentzian([fit1point5.Coefficients{1,1},fit1point5.Coefficients{2,1},fit1point5.Coefficients{3,1}],energy1point5))
figure; plot(fit1point5.Residuals{1:end,1})
%line width squared:   3.2163e-17 , error:   5.8004e-18 
%transition energy: 5.3194e-08   error:  1.8493e-22

%PEAK 2 (secretly 4)
energy2=energy(247:1:275);
peak2=mossbauer_data_no_background(247:275,1);
weight2=1./(absorption_error(247:275).^2);

c2=[(9.553*10^-16),2.672*10^-8,(10^-18)];
fit2=fitnlm(energy2, peak2, lorentzian,c2, 'Weights', weight2)
figure;
plot(energy2, peak2);
hold on;
plot(energy2,lorentzian([fit2.Coefficients{1,1},fit2.Coefficients{2,1},fit2.Coefficients{3,1}],energy2))
figure;
plot(fit2.Residuals{1:end,1});
%line width:  3.084e-17  , error:  6.4532e-18
%transition energy:  -2.5638e-08 , error:   4.3933e-23 


%PEAK 2point5 (secretly peak 5)
%this good make it look like this
energy2point5 = energy(300:1:326);
peak2point5= mossbauer_data_no_background(300:326,1);
weight2point5=1./(absorption_error(300:326).^2);

%do this for all of the peaks this is just peak 1
c2point5 = [1.452*10^-15,-1.294*10^-7,(10^-18)]; 
fit2point5=fitnlm(energy2point5, peak2point5, lorentzian,c2point5, 'Weights', weight2point5)
figure; plot(energy2point5,peak2point5)
hold on; plot(energy2point5,lorentzian([fit2point5.Coefficients{1,1},fit2point5.Coefficients{2,1},fit2point5.Coefficients{3,1}],energy2point5))
figure; plot(fit2point5.Residuals{1:end,1})
%line width squared:   3.2163e-17,  energy:  5.8004e-18  
%transition energy: -1.2969e-07  , error:  6.2515e-22  

%PEAK 3 (secretly 6)
energy3=energy(350:1:385);
peak3=mossbauer_data_no_background(350:385,1);
weight3=1./(absorption_error(350:385).^2);

c3= [1.938*10^-15,-2.34*10^-7,(10^-18)];
fit3=fitnlm(energy3,peak3,lorentzian,c3, 'Weights', weight3)
figure;
plot(energy3,peak3);
hold on;
plot(energy3, lorentzian([fit3.Coefficients{1,1},fit3.Coefficients{2,1},fit3.Coefficients{3,1}],energy3));
figure;
plot(fit3.Residuals{1:end,1})
%line width: 4.4256e-17  , error  4.1702e-18 
%transition energy: -2.3339e-07 , error:   3.8137e-22

%PEAK 3point5 (secretly peak 7)
%this good make it look like this
energy3point5 = energy(607:1:630);
peak3point5= mossbauer_data_no_background(607:630,1);
weight3point5=1./(absorption_error(607:630).^2);

%do this for all of the peaks this is just peak 1
c3point5 = [1.8580*10^-15,-2.549*10^-7,(10^-18)]; 
fit3point5=fitnlm(energy3point5, peak3point5, lorentzian,c3point5, 'Weights', weight3point5)
figure; plot(energy3point5,peak3point5)
hold on; plot(energy3point5,lorentzian([fit3point5.Coefficients{1,1},fit3point5.Coefficients{2,1},fit3point5.Coefficients{3,1}],energy3point5))
figure; plot(fit3point5.Residuals{1:end,1})
%line width:4.9968e-17   error: 3.8841e-18 
%transition energy:  -2.5481e-07, error:   1.9614e-21 

%PEAK 4 (secretly 8) (aka 2nd from the left in the right-hand batch of
%peaks)
energy4=energy(654:1:689);
peak4=mossbauer_data_no_background(654:689,1);
weight4=1./(absorption_error(654:689).^2);

c4= [1.591*10^-15,-1.522*10^-7,(10^-18)];
fit4=fitnlm(energy4,peak4,lorentzian,c4, 'Weights', weight4);
figure;
plot(energy4,peak4);
hold on;
plot(energy4, lorentzian([fit4.Coefficients{1,1},fit4.Coefficients{2,1},fit4.Coefficients{3,1}],energy4));
figure;
plot(fit4.Residuals{1:end,1})
%line width: 7.0821e-17, error: 8.8673e-18 
%transition energy: -1.529e-07, error:    1.2119e-21 

%PEAK 4point5 (secretly peak 9)
%this good make it look like this
energy4point5 = energy(715:1:737);
peak4point5= mossbauer_data_no_background(715:737,1);
weight4point5=1./(absorption_error(715:737).^2);

%do this for all of the peaks this is just peak 1
c4point5 = [1.0860*10^-15,-5.144*10^-8,(10^-18)]; 
fit4point5=fitnlm(energy4point5, peak4point5, lorentzian,c4point5, 'Weights', weight4point5)
figure; plot(energy4point5,peak4point5)
hold on; plot(energy4point5,lorentzian([fit4point5.Coefficients{1,1},fit4point5.Coefficients{2,1},fit4point5.Coefficients{3,1}],energy4point5))
figure; plot(fit4point5.Residuals{1:end,1})
%line width:  4.8543e-17 , error:   8.7178e-18  
%transition energy: -5.1092e-08, error:     1.111e-21  

%PEAK 5 (secretly 10)
energy5=energy(757:1:780);
peak5=mossbauer_data_no_background(757:780,1);
weight5=1./(absorption_error(757:780).^2);

c5= [ 1.1120*10^-15,2.64*10^-8,(10^-18)];
fit5=fitnlm(energy5,peak5,lorentzian,c5, 'Weights', weight5);
figure;
plot(energy5,peak5);
hold on;
plot(energy5, lorentzian([fit5.Coefficients{1,1},fit5.Coefficients{2,1},fit5.Coefficients{3,1}],energy5));
figure;
plot(fit5.Residuals{1:end,1})
%line width: 3.1564e-17, error: 4.9576e-18 
%transtion energy:  2.4764e-08, error:    2.5069e-22   

%PEAK 5point5 (secretly peak 11)
%this good make it look like this
energy5point5 = energy(810:1:831);
peak5point5= mossbauer_data_no_background(810:831,1);
weight5point5=1./(absorption_error(810:831).^2);

%do this for all of the peaks this is just peak 1
c5point5 = [1.6480*10^-15,1.273*10^-7,(10^-18)]; 
fit5point5=fitnlm(energy5point5, peak5point5, lorentzian,c5point5, 'Weights', weight5point5)
figure; plot(energy5point5,peak5point5)
hold on; plot(energy5point5,lorentzian([fit5point5.Coefficients{1,1},fit5point5.Coefficients{2,1},fit5point5.Coefficients{3,1}],energy5point5))
figure; plot(fit5point5.Residuals{1:end,1})
%line width:  5.4416e-17, error: 6.0448e-18  
%transition energy:  1.2712e-07, error:    8.7447e-22 

%PEAK 6 (secretly 12)
energy6=energy(860:1:894);
peak6=mossbauer_data_no_background(860:894,1);
weight6=1./(absorption_error(860:894).^2);

c6= [1.949*10^-15,2.3*10^-7,(10^-18)];
fit6=fitnlm(energy6,peak6,lorentzian,c6, 'Weights', weight6)
figure;
plot(energy6,peak6);
hold on;
plot(energy6, lorentzian([fit6.Coefficients{1,1},fit6.Coefficients{2,1},fit6.Coefficients{3,1}],energy6));
figure;
plot(fit6.Residuals{1:end,1})
%line width:5.1285e-17 , error:  4.9383e-18 
%transition energy: 2.309e-07, error: 4.9659e-22    

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
diff5and6 = trans_energy(12)-trans_energy(11); 
err5and6 = trans_energy_error(12)-trans_energy_error(11);
diff11and12 = trans_energy(6) - trans_energy(5); 
err11and12 = trans_energy_error(6)-trans_energy_error(5);

% mj -3/2 1/2
diff6and4 = trans_energy(10) - trans_energy(8);
err6and4 = trans_energy_error(10)-trans_energy_error(8);
diff10and8 = trans_energy(6) - trans_energy(4); 
err10and8 = trans_energy_error(6)-trans_energy_error(4);

% mj -1/2 1/2 
diff5and4 = trans_energy(11) - trans_energy(10);
err5and4 = trans_energy_error(11)-trans_energy_error(10);
diff10and11 = trans_energy(5) - trans_energy(4);
err10and11 = trans_energy_error(5)-trans_energy_error(4);

% put energy differences into vectors for left and right sides
ediff_mm_right = [diff5and4, diff6and4, diff5and6];
diff_mm_err_right = [err5and4, err6and4, err5and6];
ediff_mm_left = -[diff10and11, diff10and8, diff11and12]; 
diff_mm_err_left = [err10and11, err10and8, err11and12]; 

% FINISH THIS TONIGHT ANG, SO HELP ME GOD

dm_jvec = [1,2,1]; % difference in m_j per transition
I = 3/2; % excited state nuclear spin quantum number
% calculate mu
mu_mm_right = ediff_mm_right*I./(dm_jvec*H_avg*m_0); % mu in terms of mu/m_0
mu_mm_left = ediff_mm_left*I./(dm_jvec*H_avg*m_0); 
mu_mm_err_right = 



mu_weight = 1./H_error.^2; 
H_avg = (sum(Hweight.*H)./sum(Hweight)); % average H from our answers



%% ANALYSIS 1 ACTUAL VALUES
% H = 329400 Gauss
% m_1 = 1.71 * m_0 


