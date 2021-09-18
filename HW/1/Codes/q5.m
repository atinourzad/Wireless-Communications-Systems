%% Part a:

clc
clear
close all

d0 = 10;%m
numUser = 10^5;
n = 4;
p0 = 10^-6;%W

%User's coordinates:
x = -1000 + 2000 .* rand(1,numUser);
y = -1000 + 2000 .* rand(1,numUser);

%Finding user's distance:
validDist = [];%For users in the ring 
totalDist = [];%For all users
for userIndex=1:numUser
   r = sqrt(x(userIndex)^2 + y(userIndex)^2);
   totalDist = [totalDist,r];
   if (r<=1000) && (r>=10)
      validDist = [validDist,r];
   end
end
validDist = sort(validDist);
save('validDist', 'validDist');

%Calculating recieved power:
p0_dBm = 10 * log10(p0) + 30;
pR_dBmUsers = p0_dBm - 10 * n * log10(validDist/d0);%For users in the ring
pR_dBmTotal = p0_dBm - 10 * n * log10(totalDist/d0);%For all users

plot1 = cdfplot(pR_dBmUsers);
hold on
plot2 = cdfplot(pR_dBmTotal);
hold off

set(plot1, 'Color', [19/255, 206/255, 188/255]);
set(plot2, 'Color', [28/255, 152/255, 140/255]);
legend('In Ring Users', 'All Users');
title('CDF of Recieved Power');
%% Part b:

clc
clear
close all

d0 = 10;%m
D = 1000;%m
n = 4;
p0 = 10^-6;%W
N0 = -174; %dBm/Hz
BW = 10^6; %Hz

%Changing N0 unit from dBm/Hz to dBm and calculating noise power in dbm
pN_dBm = N0 + 10 * log10(BW);

%Calculating recieved power in dBm:
d = d0:0.5:D;
load('validDist.mat');

p0_dBm = 10 * log10(p0) + 30;
pR_dBmUsers = p0_dBm - 10 * n * log10(validDist/d0);%For users in the ring
pR_dBmTotal = p0_dBm - 10 * n * log10(d/d0);%For all values of d between 10 and 1000

%Calculating SNR in dB:
SNRUsers = pR_dBmUsers - pN_dBm;%For users in the ring
SNRTotal = pR_dBmTotal - pN_dBm; %For all values of d between 10 and 1000

%Plots:
figure
plot(log10(validDist), SNRUsers, 'Color', [19/255, 206/255, 188/255]);
title("SNR's Relation with Distance for Users in the Ring");
xlabel('Distance(in Logarithmic Scale)');
ylabel('SNR');

figure
plot(log10(d), SNRTotal, 'Color', [19/255, 206/255, 188/255]);
title("SNR's Relation with Distance for All Values of 10 < d < 1000");
xlabel('Distance(in Logarithmic Scale)');
ylabel('SNR');

figure
plot(log10(validDist), SNRUsers, 'Color', [19/255, 206/255, 188/255]);
hold on 
plot(log10(d), SNRTotal, 'Color', [28/255, 152/255, 140/255]);
legend('All Users in the Ring', 'All Values of 10 < d < 1000');
title("SNR's Relation with Distance");
xlabel('Distance(in Logarithmic Scale)');
ylabel('SNR');

%% Part c:

clc
clear
close all

n = 4;
p0 = 10^-6;%W
XMean = 0;
XSD = 5;%dB
d0 = 10;%m
N0 = -174; %dBm/Hz
BW = 10^6; %Hz

%%%%*********************************************************************************************************%%%%
%Finding CDF of Pr:

load('validDist.mat');

p0_dBm = 10 * log10(p0) + 30;
X_dB = normrnd(XMean, XSD, [1, length(validDist)]);
pR_dBmUsers = p0_dBm - 10 * n * log10(validDist/d0) + X_dB; 

figure
plot1 = cdfplot(pR_dBmUsers);
set(plot1, 'Color', [19/255, 206/255, 188/255]);
title('CDF of recieved power when users are in 10 < d < 1000 in presence of shadowing');


%%%%*********************************************************************************************************%%%%
%Plotting SNR:

%Changing N0 unit from dBm/Hz to dBm and calculating noise power in dbm
pN_dBm = N0 + 10 * log10(BW);

%Calculating SNR in dB:
SNR = pR_dBmUsers - pN_dBm;

figure
plot(log10(validDist), SNR, 'Color', [19/255, 206/255, 188/255]);
title("SNR's Relation with Distance in Presence of Shadowing");
xlabel('Distance in Logarithmic Scale');
ylabel('SNR');


%% Part d:

clc
clear 
close all

n = 4;
p0 = 10^-6;%W
d0 = 10;%m
D = 1000;%m
N0 = -174; %dBm/Hz
BW = 10^6; %Hz
SNRmin = 20;
XSD = 5;%dB

%Changing N0 unit from dBm/Hz to dBm and calculating noise power in dbm
pN_dBm = N0 + 10 * log10(BW);

%Calculating Pout:
load('validDist.mat');
p0_dBm = 10 * log10(p0) + 30;
pOut = zeros(1, length(validDist));

for i=1:length(validDist)
    pOut(1,i) = 1 - qfunc((SNRmin - p0_dBm + 10 * n * log10(validDist(i)/d0) + pN_dBm)/XSD);
end

plot(log10(validDist), pOut, 'Color', [19/255, 206/255, 188/255]);
title("Outage Probability VS. Distance");
xlabel('Distance(in Logarithmic Scale)');
ylabel('Outage Probability');

%% Part e:


clc
clear 
close all

numUser = 10^5;
n = 4;
p0 = 10^-6;%W
XMean = 0;
XSD = 5;%dB
d0 = 10;%m
D = 1000;%m
N0 = -174; %dBm/Hz
BW = 10^6; %Hz
SNRmin = 20; %dB
R = D;

%%%%*********************************************************************************************************%%%%
%Using Goldsmith formula:

%Changing N0 unit from dBm/Hz to dBm and calculating noise power in dBm:
pN_dBm = N0 + 10 * log10(BW);

%Calculating Pmin in dBm:
pMin = SNRmin + pN_dBm;

%Calculating Pr:
p0_dBm = 10 * log10(p0) + 30;
pRBar = p0_dBm - 10 * n * log10(R/d0);

%Calculating a an b:
a = (pMin - pRBar)/XSD;
b = 10 * n * log10(exp(1))/XSD;

%Calculating cell area:
C = qfunc(a) + exp((2 - 2 * a * b)/(b^2)) * qfunc((2 - a * b)/b);
CPercent = C * 100;

%%%%*********************************************************************************************************%%%%
%Finding distances with SNR more than 20dB:

%Calculating X:
load('validDist.mat');
X_dB = normrnd(XMean, XSD, [1, length(validDist)]);

%Calculating recieved power in dBm:
p0_dBm = 10 * log10(p0) + 30;
pR_dBmUsers = p0_dBm - 10 * n * log10(validDist/d0) + X_dB;

%Calculating SNR in dB:
SNR = pR_dBmUsers - pN_dBm;

%Finding average distance:
dists = 0;
desiredUsers = 0;
for i=1:length(SNR)
   if SNR(i) > 20
      desiredUsers = desiredUsers + 1;
   end
end

allUsers = length(SNR);

cellArea = desiredUsers * 100 / allUsers;













