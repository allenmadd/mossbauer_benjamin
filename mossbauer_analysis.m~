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


%PEAK point5 (secretly peak 1)
%this good make it look like this
energypoint5 = energy(91:1:121);
peakpoint5= mossbauer_data_no_background(91:121,1);

%do this for all of the peaks this is just peak 1
cpoint5 = [1.8330*10^-15,2.604*10^-7,(10^-18)]; %changed name to "constant 1" aka "c1"
fitpoint5=fitnlm(energypoint5, peakpoint5, lorentzian,cpoint5);
figure; plot(energypoint5,peakpoint5)
hold on; plot(energypoint5,lorentzian([fitpoint5.Coefficients{1,1},fitpoint5.Coefficients{2,1},fitpoint5.Coefficients{3,1}],energypoint5))
figure; plot(fitpoint5.Residuals{1:end,1})
%line width: 4.6256e-17, error: 4.605e-18 

%PEAK 1 (secretly peek 2)
%this good make it look like this
energy1 = energy(140:1:190);
peak1= mossbauer_data_no_background(140:190,1);

%do this for all of the peaks this is just peak 1
c1 = [1600*10^(-18),1.55*10^(-7),(10^-18)]; %changed name to "constant 1" aka "c1"
fit1=fitnlm(energy1, peak1, lorentzian,c1);
figure; plot(energy1,peak1)
hold on; plot(energy1,lorentzian([fit1.Coefficients{1,1},fit1.Coefficients{2,1},fit1.Coefficients{3,1}],energy1))
figure; plot(fit1.Residuals{1:end,1})
%line width squared for peak 1 with  error 
%{3.43978799312567e-17,3.80066459794880e-18}

%PEAK 1point5 (secretly peak 3)
%this good make it look like this
energy1point5 = energy(203:1:230);
peak1point5= mossbauer_data_no_background(203:230,1);

%do this for all of the peaks this is just peak 1
c1point5 = [1.1580*10^-15,5.316*10^-8,(10^-18)]; 
fit1point5=fitnlm(energy1point5, peak1point5, lorentzian,c1point5);
figure; plot(energy1point5,peak1point5)
hold on; plot(energy1point5,lorentzian([fit1point5.Coefficients{1,1},fit1point5.Coefficients{2,1},fit1point5.Coefficients{3,1}],energy1point5))
figure; plot(fit1point5.Residuals{1:end,1})
%line width squared:   3.2163e-17 , error:   5.8004e-18 

%PEAK 2 (secretly 4)
energy2=energy(247:1:275);
peak2=mossbauer_data_no_background(247:275,1);

c2=[(9.553*10^-16),2.672*10^-8,(10^-18)];
fit2=fitnlm(energy2, peak2, lorentzian,c2);
figure;
plot(energy2, peak2);
hold on;
plot(energy2,lorentzian([fit2.Coefficients{1,1},fit2.Coefficients{2,1},fit2.Coefficients{3,1}],energy2))
figure;
plot(fit2.Residuals{1:end,1});
%line width:  3.084e-17  , error:  6.4532e-18

%PEAK 3 (secretly 6)
energy3=energy(350:1:385);
peak3=mossbauer_data_no_background(350:385,1);

c3= [1.938*10^-15,-2.34*10^-7,(10^-18)];
fit3=fitnlm(energy3,peak3,lorentzian,c3)
figure;
plot(energy3,peak3);
hold on;
plot(energy3, lorentzian([fit3.Coefficients{1,1},fit3.Coefficients{2,1},fit3.Coefficients{3,1}],energy3));
figure;
plot(fit3.Residuals{1:end,1})
%line width: 4.4256e-17  , error  4.1702e-18 

%PEAK 4 (secretly 8) (aka 2nd from the left in the right-hand batch of
%peaks)
energy4=energy(654:1:689);
peak4=mossbauer_data_no_background(654:689,1);

c4= [1.591*10^-15,-1.522*10^-7,(10^-18)];
fit4=fitnlm(energy4,peak4,lorentzian,c4);
figure;
plot(energy4,peak4);
hold on;
plot(energy4, lorentzian([fit4.Coefficients{1,1},fit4.Coefficients{2,1},fit4.Coefficients{3,1}],energy4));
figure;
plot(fit4.Residuals{1:end,1})
%line width: 7.0821e-17, error: 8.8673e-18 



%PEAK 5 (secretly 10)
energy5=energy(757:1:780);
peak5=mossbauer_data_no_background(757:780,1);

c5= [ 1.1120*10^-15,2.64*10^-8,(10^-18)];
fit5=fitnlm(energy5,peak5,lorentzian,c5)
figure;
plot(energy5,peak5);
hold on;
plot(energy5, lorentzian([fit5.Coefficients{1,1},fit5.Coefficients{2,1},fit5.Coefficients{3,1}],energy5));
figure;
plot(fit5.Residuals{1:end,1})
%line width: 3.1564e-17, error: 4.9576e-18 

%PEAK 6 (secretly 12)
energy6=energy(860:1:894);
peak6=mossbauer_data_no_background(860:894,1);

c6= [   1.949*10^-15,2.3*10^-7,(10^-18)];
fit6=fitnlm(energy6,peak6,lorentzian,c6)
figure;
plot(energy6,peak6);
hold on;
plot(energy6, lorentzian([fit6.Coefficients{1,1},fit6.Coefficients{2,1},fit6.Coefficients{3,1}],energy6));
figure;
plot(fit6.Residuals{1:end,1})
%line width:5.1285e-17 , error:  4.9383e-18 

