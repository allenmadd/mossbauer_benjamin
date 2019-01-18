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
energy1 = energy(140:1:190);
peak1= mossbauer_data_no_background(140:190,1);

energy2=[247:1:275];
peak2=mossbauer_data_no_background(247:275,1);

lorentzian = @(b,x) ( b(1)./( (x-b(2)).^2 + b(3) ) )

b1 = [1600*10^(-18),1.55*10^(-7),(10^-18)]
fit1=fitnlm(energy1, peak1, lorentzian,b1)
figure; plot(energy1,peak1)
hold on; plot(energy1,lorentzian([fit1.Coefficients{1,1},fit1.Coefficients{2,1},fit1.Coefficients{3,1}],energy1))
figure; plot(fit1.Residuals{1:end,1})

%works but not needed:
%convert to array
%mossbauer_data= table2array(mossbauer_data);
%findpeaks(mossbauer_data);