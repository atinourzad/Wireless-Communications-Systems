%% Part a:

clc
clear
close all

%Given data:
SNR = [20 10 6];
angels = [140 100 55];

%Finding antenna weights: 
w = zeros(1,4);

for i=1:4
  w(1,i) = sum(sqrt(10.^(SNR/10)) .* exp(-1i * pi * (i-1) * sin(angels*pi/180)));
end    

%Plotting antenna pattern:
antennaIndex = [1 2 3 4];
teta = 0:5:360;
A = zeros(1, length(teta));

for i=1:length(teta)
    A(i) = sum(w .* exp(1i * pi * (antennaIndex-1) * sin(teta(i)*pi/180)));
end

polarplot(teta*pi/180, abs(A), 'Color', [19/255, 206/255, 188/255]);

%Finding total SNR:
num = 0;
signalsEnergy = 10.^(SNR/10);

for i=1:3
    num = num + signalsEnergy(1,i) * abs(sum(w .* exp(1i * pi * (antennaIndex-1) * sin(angels(i)*pi/180))))^2;
end

denum = sum(abs(w).^2);

SNR = num/denum;
fprintf('The total SNR is: %f (dB) \n ',10 * log10(SNR)); 


%% Part b:

clc
clear
close all

%Given data:
incomingAngles = [45; 0; 330; 120];

%Creating equations for finding weights:
A = zeros(4);
for i=1:4
  A(:,i) = exp(1i * pi * (i-1) * sin(incomingAngles*pi/180)); 
end

B = zeros(4,1);
B(1,1) = 1;

w = linsolve(A,B);
w = w.';
        
%Plotting antenna pattern:
antennaIndex = [1 2 3 4];
teta = 0:5:360;
A = zeros(1, length(teta));

for i=1:length(teta)
    A(i) = sum(w .* exp(1i * pi * (antennaIndex-1) * sin(teta(i)*pi/180)));
end

polarplot(teta*pi/180, abs(A), 'Color', [19/255, 206/255, 188/255]);

%Finding total SNR:
num = 0;

SNRsdB = [20 14 11 10];
signalsEnergy = 10.^(SNRsdB/10);

for i=1:4
    num = num + signalsEnergy(1,i) * abs(sum(w .* exp(1i * pi * (antennaIndex-1) * sin(incomingAngles(i)*pi/180))))^2;
end

denum1 = sum(abs(w).^2);

SNR = num/denum1;
fprintf('The total SNR is: %f (dB) \n ',10 * log10(SNR));

%% Part c:

clc
clear
close all

%Given data:
incomingAngles = [0 150];
signalEnergy = 10 ^ (20/10);

%Calculating phase difference between antennas:
phis = zeros(2,4);
antennaIndex = [1 2 3 4];

for i=1:2
   phis(i,:) = pi * (antennaIndex-1) * sin(incomingAngles(i)*pi/180);
end
phase = exp(1i .* phis);

%Defiening the goal function to minimize:
SINRmin = @(w) (abs(w(1)+w(2)*phase(2,2)+w(3)*phase(2,3)+w(4)*phase(2,4))^2+ 0.01 *(abs(w(1))^2+abs(w(2))^2+abs(w(3))^2+abs(w(4))^2))/(abs(w(1)+w(2)+w(3)+w(4))^2);
wAlgo1 = fminsearch(SINRmin, [1 1 1 1]);
wAlgo2 = ga(SINRmin, 4);

%Plotting the antenna pattern:
antennaIndex = [1 2 3 4];
teta = 0:5:360;
AAlgo1 = zeros(1, length(teta));
AAlgo2 = zeros(1, length(teta));

for i=1:length(teta)
    AAlgo1(i) = sum(wAlgo1 .* exp(1i * pi * (antennaIndex-1) * sin(teta(i)*pi/180)));
    AAlgo2(i) = sum(wAlgo2 .* exp(1i * pi * (antennaIndex-1) * sin(teta(i)*pi/180)));
end

figure
polarplot(teta*pi/180, abs(AAlgo1), 'Color', [19/255, 206/255, 188/255]);
title('Antenna Pattern Using Nelder-Mead Algorithm');

figure
polarplot(teta*pi/180, abs(AAlgo2), 'Color', [28/255, 152/255, 140/255]);
title('Antenna Pattern Using Genetic Algorithm');


%Calculating total SNR:
num1 = signalEnergy * abs(sum(wAlgo1 .* exp(1i * pi * (antennaIndex-1) * sin(incomingAngles(1)*pi/180))))^2;
denum1 = signalEnergy * abs(sum(wAlgo1 .* exp(1i * pi * (antennaIndex-1) * sin(incomingAngles(2)*pi/180))))^2 + sum(abs(wAlgo1).^2);

num2 = signalEnergy * abs(sum(wAlgo2 .* exp(1i * pi * (antennaIndex-1) * sin(incomingAngles(1)*pi/180))))^2;
denum2 = signalEnergy * abs(sum(wAlgo2 .* exp(1i * pi * (antennaIndex-1) * sin(incomingAngles(2)*pi/180))))^2 + sum(abs(wAlgo2).^2);

SNR1 = num1/denum1;
fprintf('The total SINR using Nelder-Mead algorithm is: %f (dB) \n ',10 * log10(SNR1));

SNR2 = num2/denum2;
fprintf('The total SINR using genetic algorithm is: %f (dB) \n ',10 * log10(SNR2));
