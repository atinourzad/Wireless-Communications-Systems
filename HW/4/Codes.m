%% a
clc
clear 

hSquareAvg = [1.5 0.75 0.5];
numChannel = length(hSquareAvg);
signalPower = 6;

N0 = 0.2 * 10^(-6);
Bc = 3*10^3;
noisePower = N0 * Bc;

signal = 0.001;
hSquare = signal:0.01:10*max(hSquareAvg);
gamma = hSquare * (signalPower/noisePower);

pGamma = zeros(numChannel, length(hSquare));
for i=1:numChannel
    pGamma(i,:) = (1/hSquareAvg(i)) * exp(-hSquare/hSquareAvg(i));
end

%Integration:
gamma0v = 1:0.1:10;
sumP = zeros(length(gamma0v),1);

for j=1:length(gamma0v)
   gamma0 = gamma0v(j);
   sumP(j) = 0;
   
   for i = 1:numChannel
      a = gamma .* (gamma>gamma0);
      [maxValue,maxIndex] = max(a>0);
      gammaMax = a(find(a));
      pGammaMax = pGamma(i, maxIndex:length(gamma));
      PjP = (1/gamma0) - (1./gammaMax);
      sumP(j)= sumP(j) + sum(PjP.*pGammaMax) * signal;
   end
end

[minPowerValue, minPowerIndex] = min(abs(sumP-1));
gamma0Channel = gamma0v(minPowerIndex);

C = 0;
for i = 1:length(numChannel)
    a = gamma.*(gamma>gamma0Channel);
    [maxValue,maxIndex] = max(a>0);
    gammaMax = a(find(a));
    pGammaMax = pGamma(i,maxIndex:length(gamma));
    C = C + Bc * signal * sum(log2(gammaMax/gamma0Channel).*pGammaMax);
end

fprintf('channel capacity is %f. \n ',C); 

%% B
clc
clear

hSquareAvg = [1.5 0.75 0.5];
numChannel = length(hSquareAvg);
signalPower = 6;

N0 = 0.2 * 10^(-6);
Bc = 3*10^3;
noisePower = N0 * Bc;

signal = 0.001;

P = signalPower/3;

C = zeros(numChannel,1);
for k = 1:numChannel
    gammaBar = hSquareAvg(k);
    hSquare = signal:signal:10*gammaBar;
    gamma = hSquare*(P/noisePower);
    pGamma = (1/gammaBar)*exp(-hSquare/gammaBar);
    gamma0v = 0.01:0.01:1;
    
    for j = 1:length(gamma0v)
        gamma0 = gamma0v(j);
        a = gamma.*(gamma>gamma0);
        [maxValue,maxIndex] = max(a>0);
        gammaMax = a(find(a));
        pGammaMax = pGamma(maxIndex:length(gamma));
        PjP = (1/gamma0)-(1./gammaMax);
        sumP(j) = sum(PjP.*pGammaMax)*signal;
    end
    [minPowerValue, minPowerIndex] = min(abs(sumP-1));
    gamma0Channel = gamma0v(minPowerIndex);
    a = gamma.*(gamma>gamma0Channel);
    [maxValue,maxIndex] = max(a>0);
    gammaMax = a(find(a));
    pGammaMax = pGamma(maxIndex:length(gamma));
    C(k) = Bc * signal * sum(log2(gammaMax/gamma0Channel).*pGammaMax);
end
Ctot = sum(C);
fprintf('channel capacity is %f. \n ',Ctot); 