%% Part a:

clc
clear
close all

numSamples = 10^5;
numClusters = 10;
maxDelay = 10;
minDelay = 1;

%Creating gains:
gMean = 0;
g = zeros(numClusters, numSamples);

for clusterIndex = 1:numClusters
    %Creating delays:
    tou = minDelay + (maxDelay - minDelay) * rand(1, numSamples);
    
    %Creating average power gain:
    sigma2 = (10^(-3)) * (tou.^(-4));
    
    gSD = sqrt(sigma2/2);
    g(clusterIndex,:) = normrnd(gMean, gSD) + 1i * normrnd(gMean, gSD);
end

%Overall channel gain for each realization:
gTotal = sum(g,1);

%Computing power gain(|h|^2):
h = gTotal;
hMag = abs(h) .^ 2;

%Computing average power gain(E{|h|^2}):
hMagAvg = sum(hMag,2) / numSamples;
sigmaAvg = 10 * sum(sigma2,2) / numSamples;
fprintf('The simulated power gain is %f, while the actual its value is %f.  \n ',hMagAvg, sigmaAvg); 

plot1 = cdfplot(hMag);
set(plot1, 'Color', [19/255, 206/255, 188/255]);
title('CDF of Average Power Gain');

save('gains', 'h', 'g');
save('power gain', 'hMag');
%% Part b:

clc
clear
close all

numClusters = 10;
maxFreq = 10^6;
maxDelay = 10;
minDelay = 1;

load('gains.mat');

%Calculating H(f):
oneSampleG = g(:,1);
f = 0:1000:maxFreq-1;
H = 0;

tou = minDelay + (maxDelay - minDelay) * rand(10,1);
for clusterIndex = 1:numClusters
    H = H + oneSampleG(clusterIndex,1) * exp(-1i * 2 * pi * f * tou(clusterIndex,1));
end

plot(f, real(H), 'Color', [19/255, 206/255, 188/255]);
xlabel('f');
ylabel('H(f)');
title('Channel Responce');

%% Part c:

clc
clear
close all

numSamples = 10^5;
numClusters = 10;
maxDelay = 10;
minDelay = 1;
v = 20;
fc = 3 * 10^9;
c = 3 * 10^8;

%Using gains from part a:
load('gains.mat');

%Creating different incoming rays' angle:
teta = 2 * pi * rand(numClusters,1);

%Creating Doppler frequency for all clusters:
fd = fc * v * cos(teta) / c;

%Time interval:
time = 0:10^-4:1;
timeLen = length(time);

%Overall channel gain for each realization with Doppler effect:
oneSampleG = g(:,1);
gWithDop = 0;

for clusterIndex = 1:numClusters
   gWithDop =  gWithDop + g(clusterIndex,1) .* exp(1i * 2 * pi * fd(clusterIndex, 1) * time);
end

%Computing power gain(|h|^2):
hWithDop = gWithDop;
hMagWithDop = abs(hWithDop) .^ 2;

%Plotting power gain for one realization with Doppler effect:
figure
plot(time, abs(gWithDop), 'Color', [19/255, 206/255, 188/255]);
xlabel('Time');
ylabel('Channel Gain');

%Computing average power gain(E{|h|^2}):
hMagAvg = sum(hMagWithDop,2) / timeLen;
fprintf('The simulated power gain with doppler gain is %f \n',hMagAvg); 

load('power gain');

figure
plot1 = cdfplot(hMagWithDop);
hold on
plot2 = cdfplot(hMag);
hold off

set(plot1, 'Color', [19/255, 206/255, 188/255]);
set(plot2, 'Color', [28/255, 152/255, 140/255]);
legend('With Doppler', 'Without Doppler');
title('CDF of Average Power Gain');
