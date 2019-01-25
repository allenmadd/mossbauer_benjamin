%Madeleine and Angela try their bests
%Mossbauer Analysis 17 Jan 2019


%files found at path: "C:\Users\Student\Documents\MadelAngela_mossbauer_15Jan2019";
close all
clear all

%% Analysis: Step 2: Assume Gaussian Distribution
% Calculate the error in your data with and without the backgound counts.

% Upload the file as a table. 
filename="step15_22Jan2019.csv";
mossbauer_data= readtable(filename);

% Graph the data
mossbauer_data= removevars(mossbauer_data, {'Energy', 'Channel'});
stackedplot(mossbauer_data,'.-');
% Rename the data as an array. 
mossbauer_data= table2array(mossbauer_data);

% Calculate the error in the data (1) with and (2) without background 
% counts. 
filename= "step7_22Jan2019.csv";
background_counts = readtable(filename);


% 1) Error including background counts. Take raw data as collected in step
% 7. 
transmission_error= sqrt(mossbauer_data); % By definition, since these data values are counts
figure; % graph error bars on data to visualize it
errorbar([1:1:1024],mossbauer_data,transmission_error);
title('Transmission Error (Raw Data with Background Counts)'); xlabel('Channel'); ylabel('Counts');

% 2) Error taking away background counts. 
% First take some random portion of data with no troughs in it. We select 
% region of channels 512 to 550. This area represents channels where no 
% absorption happened, the recorded data is only background counts. We use
% this section to represent all background counts in the rest of the data. 
background_counts = mean(mossbauer_data(512:550,1)); 
background_counts_error = sqrt(background_counts);

absorption_error= sqrt(background_counts+mossbauer_data);
mossbauer_data_no_background =  background_counts- mossbauer_data;
figure;
errorbar([1:1:1024], mossbauer_data_no_background, absorption_error);
title('Absorption Error (no background counts)'); xlabel('Channel'); ylabel('Counts');

%% Convert channels to energy. Account for Doppler Effect.
% Given miniumum (9.59 mm/s) and maximum velocity (10.22 mm/s) values for 
% the speed of the source. Find the slope of velocity. 
dv = (10.22+9.59)/(500);
% We know velocity starts at it's minimum and increases to reach its 
% maximum at the 500th channel, where it then decreases. Assign velocity 
% values that correspond to channels. 
% Doppler Effect: Use E = E_0*v/c to find the shift in the gamma ray energy
% based on our knowledge of the speed of the source (v). E_0 = 14.4 keV.
% (the peak we recognized for our discriminator in experimental setup). 
for i = 1:1:500
    v1(i) = -9.59 + (i-1)*dv;
    E1(i)= 14.4*(v1(i)/(3*10^8)); % decreasing energy in graph
    v2(i)=10.22 - (i-1)*dv;
    E2(i) = 14.4*(v2(i)/(3*10^8)); % increasing energy in graph
end

v = [v1,v2]; % combine velocity vectors to make one velocity vector
resolution_vel = max(dv./v(715:737)); % check the resolution of velocity vector. We get 0.0456. This will be important in later calculations



% There is a dead zone between channels 500-512 where the source is
% changing directions. Need to account for this when creating energy
% vector. 
old_channel= [1:1:1024]; % define vector for old channel
energy = [E2(13:end), zeros(1,23), E1, E2(1:13)];
figure; % graph results to visualize energy
plot(old_channel, energy,'o-');

%% Analysis: Step 3: Lorentzian distribution 
% The distribution of counts is actually a Lorentzian distribution. Use the
% lorentzian function below to analyze the following peaks and create a
% nonlinear fit for each peak. The nonlinear fit through matlab is
% important because it will allow us to extrapolate variables to apply to
% further analysis in step 1, which follows. 

% This is the Lorentzian function. 
lorentzian = @(b,x) ( b(1)./( (x-b(2)).^2 + b(3) ) )

% b(1) is a multiplicative prefactor
% b(2) is the energy value at which the Lorentzian peaks
% b(3) = (linewidth/2)^2

%% PEAK ANALYSIS
%notes from Jan 20th: it looks like we did every other peak in class? why?
%i will just do the even peaks and then if necessary i can do the others

%PEAK point5 (secretly peak 1)
%this good make it look like this
energypoint5 = energy(81:1:130);
peakpoint5= mossbauer_data_no_background(81:130,1);
weightpoint5 = 1./(absorption_error(81:130).^2);
%do this for all of the peaks this is just peak 1
cpoint5 = [2942*10^-18,2.623*10^-7,(10^-18)]; %changed name to "constant 1" aka "c1"
fitpoint5=fitnlm(energypoint5, peakpoint5, lorentzian,cpoint5,'Weights',weightpoint5)
figure; plot(energypoint5,peakpoint5)
hold on; plot(energypoint5,lorentzian([fitpoint5.Coefficients{1,1},fitpoint5.Coefficients{2,1},fitpoint5.Coefficients{3,1}],energypoint5))
figure; plot(fitpoint5.Residuals{1:end,1})
%line width squared: 4.5027e-17    3.1648e-18
%delta E (transition energy): 2.6098e-07    3.1596e-22

%PEAK 1 (secretly peek a boo 2)
%this good make it look like this
energy1 = energy(140:1:181);
peak1= mossbauer_data_no_background(140:181,1);
weight1= 1./(absorption_error(140:181).^2);
%do this for all of the peaks this is just peak 1
c1 = [2671*10^(-18),1.5578*10^(-7),(10^-18)]; %changed name to "constant 1" aka "c1"
fit1=fitnlm(energy1, peak1, lorentzian,c1, 'Weights', weight1)
figure; plot(energy1,peak1)
hold on; plot(energy1,lorentzian([fit1.Coefficients{1,1},fit1.Coefficients{2,1},fit1.Coefficients{3,1}],energy1))
figure; plot(fit1.Residuals{1:end,1})
%line width squared for peak 1 with  error  3.2821e-17    3.6489e-18
%transition energy: 1.5722e-07    1.4523e-23


%PEAK 1point5 (secretly peak 3)
%this good make it look like this
energy1point5 = energy(204:1:237);
peak1point5= mossbauer_data_no_background(204:237,1);
weight1point5=1./(absorption_error(204:237).^2);
%do this for all of the peaks this is just peak 1
c1point5 = [1559*10^-18,5.125*10^-8,(10^-18)]; 
fit1point5=fitnlm(energy1point5, peak1point5, lorentzian,c1point5, 'Weights', weight1point5)
figure; plot(energy1point5,peak1point5)
hold on; plot(energy1point5,lorentzian([fit1point5.Coefficients{1,1},fit1point5.Coefficients{2,1},fit1point5.Coefficients{3,1}],energy1point5))
figure; plot(fit1point5.Residuals{1:end,1})
%line width squared:   3.4055e-17     5.939e-18
%transition energy: 5.3576e-08    1.7534e-22 

%PEAK 2 (secretly 4)
energy2=energy(247:1:275);
peak2=mossbauer_data_no_background(247:275,1);
weight2=1./(absorption_error(247:275).^2);
c2=[(1759*10^-18),-2.672*10^-8,(10^-18)];
fit2=fitnlm(energy2, peak2, lorentzian,c2, 'Weights', weight2)
figure;
plot(energy2, peak2);
hold on;
plot(energy2,lorentzian([fit2.Coefficients{1,1},fit2.Coefficients{2,1},fit2.Coefficients{3,1}],energy2))
figure;
plot(fit2.Residuals{1:end,1});
%line width:  2.4694e-17    3.3434e-18  
%transition energy: -2.5176e-08    6.8182e-24  


%PEAK 2point5 (secretly peak 5)
%this good make it look like this
energy2point5 = energy(292:1:330);
peak2point5= mossbauer_data_no_background(292:330,1);
weight2point5=1./(absorption_error(292:330).^2);


c2point5 = [2399*10^-18,-1.294*10^-7,(10^-18)]; 
fit2point5=fitnlm(energy2point5, peak2point5, lorentzian,c2point5, 'Weights', weight2point5)
figure; plot(energy2point5,peak2point5)
hold on; plot(energy2point5,lorentzian([fit2point5.Coefficients{1,1},fit2point5.Coefficients{2,1},fit2point5.Coefficients{3,1}],energy2point5))
figure; plot(fit2point5.Residuals{1:end,1})
%line width squared:   3.3312e-17    2.8551e-18
%transition energy:  -1.2949e-07    2.5224e-22 

%PEAK 3 (secretly 6)
energy3=energy(346:1:385);
peak3=mossbauer_data_no_background(346:385,1);
weight3=1./(absorption_error(346:385).^2);

c3= [2981*10^-18,-2.359*10^-7,(10^-18)];
fit3=fitnlm(energy3,peak3,lorentzian,c3, 'Weights', weight3)
figure;
plot(energy3,peak3);
hold on;
plot(energy3, lorentzian([fit3.Coefficients{1,1},fit3.Coefficients{2,1},fit3.Coefficients{3,1}],energy3));
figure;
plot(fit3.Residuals{1:end,1})
%line width: 4.3652e-17    3.3551e-18 
%transition energy: -2.3369e-07    1.0728e-21

%PEAK 3point5 (secretly peak 7)
%this good make it look like this
energy3point5 = energy(595:1:635);
peak3point5= mossbauer_data_no_background(595:635,1);
weight3point5=1./(absorption_error(595:635).^2);

%do this for all of the peaks this is just peak 1
c3point5 = [3068*10^-18,-2.549*10^-7,(10^-18)]; 
fit3point5=fitnlm(energy3point5, peak3point5, lorentzian,c3point5, 'Weights', weight3point5)
figure; plot(energy3point5,peak3point5)
hold on; plot(energy3point5,lorentzian([fit3point5.Coefficients{1,1},fit3point5.Coefficients{2,1},fit3point5.Coefficients{3,1}],energy3point5))
figure; plot(fit3point5.Residuals{1:end,1})
%line width:5.5621e-17    3.6109e-18
%transition energy: -2.5558e-07    1.4417e-21 

%PEAK 4 (secretly 8) (aka 2nd from the left in the right-hand batch of
%peaks)
energy4=energy(649:1:689);
peak4=mossbauer_data_no_background(649:689,1);
weight4=1./(absorption_error(649:689).^2);

c4= [2626*10^-18,-1.522*10^-7,(10^-18)];
fit4=fitnlm(energy4,peak4,lorentzian,c4, 'Weights', weight4)
figure;
plot(energy4,peak4);
hold on;
plot(energy4, lorentzian([fit4.Coefficients{1,1},fit4.Coefficients{2,1},fit4.Coefficients{3,1}],energy4));
figure;
plot(fit4.Residuals{1:end,1})
%line width: 4.1229e-17    3.6712e-18 
%transition energy: -1.529e-07     4.593e-22

%PEAK 4point5 (secretly peak 9)
%this good make it look like this
energy4point5 = energy(705:1:747);
peak4point5= mossbauer_data_no_background(705:747,1);
weight4point5=1./(absorption_error(705:747).^2);

%do this for all of the peaks this is just peak 1
c4point5 = [1663*10^-18,-5.144*10^-8,(10^-18)]; 
fit4point5=fitnlm(energy4point5, peak4point5, lorentzian,c4point5, 'Weights', weight4point5)
figure; plot(energy4point5,peak4point5)
hold on; plot(energy4point5,lorentzian([fit4point5.Coefficients{1,1},fit4point5.Coefficients{2,1},fit4point5.Coefficients{3,1}],energy4point5))
figure; plot(fit4point5.Residuals{1:end,1})
%line width: 5.2732e-17    6.8711e-18 
%transition energy: -5.0651e-08    2.2912e-22   

%PEAK 5 (secretly 10)
energy5=energy(747:1:788);
peak5=mossbauer_data_no_background(747:788,1);
weight5=1./(absorption_error(747:788).^2);

c5= [ 1528*10^-18,2.463*10^-8,(10^-18)];
fit5=fitnlm(energy5,peak5,lorentzian,c5, 'Weights', weight5)
figure;
plot(energy5,peak5);
hold on;
plot(energy5, lorentzian([fit5.Coefficients{1,1},fit5.Coefficients{2,1},fit5.Coefficients{3,1}],energy5));
figure;
plot(fit5.Residuals{1:end,1})
%line width: 4.9525e-17      7.51e-18
%transtion energy: 2.5373e-08    5.8269e-23 

%PEAK 5point5 (secretly peak 11)
%this good make it look like this
energy5point5 = energy(811:1:841);
peak5point5= mossbauer_data_no_background(811:841,1);
weight5point5=1./(absorption_error(811:841).^2);


c5point5 = [2554*10^-18,1.292*10^-7,(10^-18)]; 
fit5point5=fitnlm(energy5point5, peak5point5, lorentzian,c5point5, 'Weights', weight5point5)
figure; plot(energy5point5,peak5point5)
hold on; plot(energy5point5,lorentzian([fit5point5.Coefficients{1,1},fit5point5.Coefficients{2,1},fit5point5.Coefficients{3,1}],energy5point5))
figure; plot(fit5point5.Residuals{1:end,1})
%line width: 4.2769e-17     3.671e-18  
%transition energy:   1.276e-07    2.6739e-21  

%PEAK 6 (secretly 12)
energy6=energy(850:1:894);
peak6=mossbauer_data_no_background(850:894,1);
weight6=1./(absorption_error(850:894).^2);

c6= [1639*10^-18,2.281*10^-7,(10^-18)];
fit6=fitnlm(energy6,peak6,lorentzian,c6, 'Weights', weight6)
figure;
plot(energy6,peak6);
hold on;
plot(energy6, lorentzian([fit6.Coefficients{1,1},fit6.Coefficients{2,1},fit6.Coefficients{3,1}],energy6));
figure;
plot(fit6.Residuals{1:end,1})
%line width:  5.0837e-17    3.5296e-18 
%transition energy: 2.3101e-07     4.889e-22  


%% Analysis 1 JAN 22 2019 LAB DAY 3
% Measure change in energy that corresponds to transition energy.  
% The transition energy corresponds to the x coordinate of the max peak in
% our nonlinear fit we created. Matlab recorded this value as b2 in the 
% info about the fit. 
% Manually take b2 values and put into vector as transition energies: 
trans_energy = [2.6072e-07,1.5731e-07,5.3194e-08,-2.5638e-08,-1.2969e-07,-2.3339e-07,-2.5481e-07, -1.529e-07, -5.1092e-08, 2.4764e-08, 1.2712e-07, 2.309e-07]; 
trans_energy_error = [1.8672e-21,  3.7132e-23, 1.8493e-22, 4.3933e-23, 6.2515e-22,  3.8137e-22, 1.9614e-21, 1.2119e-21, 1.111e-21, 2.5069e-22, 8.7447e-22, 4.9659e-22];

% Visualize transition energies per peak
figure;
title('Transition Energies per peak');
ylabel('Transition Energy (eV)'); xlabel('Peak Number'); 
errorbar((1:6),trans_energy(1:6),trans_energy_error(1:6),'ro-'); grid on; hold on;
errorbar((7:12),trans_energy(7:12),trans_energy_error(7:12),'b*-'); 
legend('Left Side, Peaks 7-12','Right Side, Peaks 1-6');

% Consider differences between data sets. Visualize it for conceptual
% understanding. 
transition_diff = trans_energy(1:6) - trans_energy(12:-1:7); % difference between sets of data
figure;
plot((1:6),transition_diff,'o-'); grid on; 
title('Difference between left and right data sets');

% range of differences between left and right data sets (2.1-3)*10^-8. Aka
% differences in data are one order of magnitude smaller than data itself. 

%% Analysis 1: Part B
% Deduce the size of the internal magnetic field felt by the nucelus
% (H_avg). 
% We need H_avg in calculations for part a, so do part b first. 

% Find the difference in energy between pairs that end at the same energy 
% level, but start at different energy levels. 
% These pairs include peaks: (2,4) (3,5) (9,11) (10,8) 
% We can use these energy differences to find H_avg since we already know 
% the ground state magnetic moment (m_0) and I.
d_energy_h = ([trans_energy(9)- trans_energy(11), trans_energy(8) - trans_energy(10), trans_energy(5)-trans_energy(3), trans_energy(4)-trans_energy(2)]); % energy differences between m_j starting points
d_energy_h_error = [sqrt((trans_energy_error(11))^2+trans_energy_error(9)^2), sqrt((trans_energy_error(10))^2+(trans_energy_error(8))^2), sqrt((trans_energy_error(5))^2+(trans_energy_error(3))^2),sqrt((trans_energy_error(2))^2+(trans_energy_error(4))^2)]; % error in energy differences





% Use these energy differences to find magnetic field H. 
m_0 = 2.93 *10^-13; % Given ground state magnetic moment. 
I = 1/2; % given spin of the nucleus
dm_j = 1; % difference between m_j states for starting positions. same for each pair. 
H = -I*(1/(m_0*dm_j))*d_energy_h; 
% H_error = -I*(1/(m_0*dm_j))*1.0e-20*[0.141,0.124,0.0652,0.00575]; % manually taking 3 sig figs into account yet because matlab doesn't have a function for that, apparently.

% Check resolution of H to see if its reasonably sized. 
% resolution_H = H_error./H 
% We find out the resolution of H, from transition energy error values
% (given as SE for b2 value by matlab), is super small. The H resolution is
% not reasonable. Instead use the resolution set by the velocity, already found above.  

% Redefine H_error based on limits of system, which is velocity resolution
% in this case. 
H_error = H*resolution_vel;

% H is 1.0e+05 * [3.0412    3.0318    3.1209    3.1220] Gauss
% H_error 1.0e+04 * [1.3855    1.3812    1.4218    1.4223] Gauss 



% Graph H values found based on energy pairs. For visualization purposes. 
figure; 
errorbar(1:4,H,H_error,'o'); grid on; xlabel('pair'); ylabel('H Magnetic field (Gauss)');
title('Magnetic Field per pair'); 
% Reflection: all the error bars are overlapping, so we cannot say that the
% magnetic fields are different. They are supposed to be the same magnetic 
% field, so this is good. 

% Weigh results by error and find average H: 
Hweight = 1./H_error.^2; 
H_avg = (sum(Hweight.*H)./sum(Hweight)); % average H from our answers

% RESULTS! 
% H_avg = 3.0778e+05
% Compare to known H = 329400 Gauss

%% Analysis 1 Part A
% Deduce the magnetic moment of the first excited state of Fe57. 
% First find the energy differences between the following mj transitions in
% the excited state where I= 3/2. 
% The pairs energy differences can be found on the right side are: (5,6),
% (6,4), and (4,5). The pairs on the left side are (11,12), (12,10),
% (10,11). 
% Note that these energy transitions have different dm_j transitions. 

% Also note: again we are using velocity resolution to define mu_error
% because velocity resolution more limiting factor. 
% Pair 1: mj -3/2 to -1/2
diff5and6 = abs(trans_energy(12)-trans_energy(11)); 
diff11and12 = abs(trans_energy(6) - trans_energy(5)); 

% Pair 2: mj -3/2 to 1/2
diff6and4 = abs(trans_energy(10) - trans_energy(8));
diff10and8 = abs(trans_energy(6) - trans_energy(4)); 

% Pair 3: mj -1/2 to 1/2 
diff5and4 = abs(trans_energy(11) - trans_energy(10));
diff10and11 = abs(trans_energy(5) - trans_energy(4));

% Organize energy differences based on left and right data sets. 
% Since there is less rhyme and reason in dealing with these pairs, I kept 
% the right and left data sets separate to help me from getting 
% disorganized. 
ediff_mm_right = [diff5and4, diff6and4, diff5and6];
ediff_mm_left = [diff10and11, diff10and8, diff11and12]; 

% Calculate magnetic moment, mu. Calculated for left and right data sets. 
dm_jvec = [1,2,1]; % Create vector for difference in m_j per transition 
I = 3/2; % Excited state nuclear spin quantum number
mu_mm_right = -ediff_mm_right*I./(dm_jvec*H_avg*m_0); % mu in terms of mu/m_0
mu_mm_left = -ediff_mm_left*I./(dm_jvec*H_avg*m_0); 


% After seeing the calculations were done correctly for left and right data 
% sets individually, I combine the data sets in the end. 
mu_mm = [mu_mm_left, mu_mm_right];
mu_mm_err = mu_mm*resolution_vel/m_0;

H_error = H*resolution_vel

% Weigh results by error and find average mu: 
mu_weight = 1./mu_mm_err.^2; 
mu_avg = (sum(mu_weight.*mu_mm)./sum(mu_weight)); % average mu/m_0 

% RESULTS! 
% mu_avg = -1.7090

% Consider probability that the accepted mu_1 value falls within our error
% bars around mu_avg. Compare to known mu_1 = 1.71 * m_0.

mu_1 = 1.71; % Accepted value in terms of 1/m_0. 
check_mu = abs(mu_avg - mu_1)./(mu_mm_err);

