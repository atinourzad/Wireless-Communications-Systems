%%Part a

clc
clear
close all

numOfBits = 10^4;
symbolsRef = [0 0 0 0;0 1 1 2;0 2 3 1;0 3 2 3;1 0 2 1;1 1 3 3;1 2 1 0;1 3 0 2;2 0 1 3;2 1 0 1;2 2 2 2;2 3 3 0;3 0 3 2;3 1 2 0;3 2 0 3;3 3 1 1];

%Creating random 0 and 1 bits:
bits = randi([0,1], [numOfBits,4]);

%Creating coded symbols:
codedSymbols = zeros(numOfBits, 4);

for i=1:numOfBits
    bitBlock = bi2de(bits(i,:), 'left-msb');
    codedSymbols(i,:) = symbolsRef(bitBlock+1,:);
end

%Creating modulated symbols with QPSK mudolation: 0->1, 1->j, 2->-1, 3->-j 
modulatedSymbols = exp(1i*pi*codedSymbols/2);

%Creating fading gains:
hMean = 0;
hSD = sqrt(1/2);
h = normrnd(hMean, hSD, [numOfBits,4]) + 1i * normrnd(hMean, hSD, [numOfBits,4]);

%Modulated symbols go through fading gain:
fadedSymbols = modulatedSymbols .* h;
%sd = std(fadedSymbols);

%Modulating refrence symbols:
modulatedSymbolsRef = exp(1i*pi*symbolsRef/2);

%%%-----------------------------------------------Finding Probability of Error----------------------------------------%%%

%Finding Different SNRs:
SNRdB = 1:22;
SNR = 10 .^ (SNRdB/10);

pSymbolErrorFaded = zeros(length(SNR), 1);
pSymbolErrorNoisy = zeros(length(SNR), 1);

for i=1:length(SNR)
    
    %Creating AWGN noise:
    noiseMean = 0;
    noiseSD = sqrt(1/(2*SNR(i)));
    noise = normrnd(noiseMean, noiseSD, [numOfBits,4]) + 1i * normrnd(noiseMean, noiseSD, [numOfBits,4]);
    
    %Faded signals go through AWGN channel:
    fadedY = fadedSymbols + noise;
    
    %Symbols go through AWGN channel:
    noisyY = modulatedSymbols + noise;
    
    %Using ML
    detectedFadedSymbols = zeros(numOfBits,4);
    detectedNoisySymbols = zeros(numOfBits,4);
    
    for j=1:numOfBits
      minDistFaded = inf;
      minDistNoisy = inf;
      
      for k=1:16
          %Modulated refreces go through fading:
          fadedRef = modulatedSymbolsRef(k,:) .* h(j,:);
          
          %Finding argmin:
          distFaded = norm(fadedY(j,:) - fadedRef);
          distNoisy = norm(noisyY(j,:) - modulatedSymbolsRef(k,:));
          
          if distFaded < minDistFaded
              minDistFaded = distFaded;
              fadedSymbolIndex = k;
          end
          
          if distNoisy < minDistNoisy
              minDistNoisy = distNoisy;
              noisySymbolIndex = k;
          end
      end
      detectedFadedSymbols(j,:) = symbolsRef(fadedSymbolIndex,:);
      detectedNoisySymbols(j,:) = symbolsRef(noisySymbolIndex,:);
    end 
    pSymbolErrorFaded(i,1) = nnz(codedSymbols - detectedFadedSymbols)/numOfBits;
    pSymbolErrorNoisy(i,1) = nnz(codedSymbols - detectedNoisySymbols)/numOfBits;
   
end

figure 
semilogy(SNRdB, (pSymbolErrorFaded), 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, (pSymbolErrorNoisy), 'Color', [28/255, 152/255, 140/255]);
hold off
xlabel('SNR (dB)');
ylabel('P_e');
legend('With fading', 'Without fading');

%% Part b:
clc
%Using performance curve to find diversity order:

for i=length(pSymbolErrorFaded):-1:1
   if pSymbolErrorFaded(i) ~= 0 && pSymbolErrorFaded(i-1) ~= 0 && pSymbolErrorFaded(i) ~= pSymbolErrorFaded(i-1)
      index = i;
      break
   end
end

LPlot = floor(abs(10*log10(pSymbolErrorFaded(index)) - 10*log10(pSymbolErrorFaded(index-1))));

%Using system analysis to find diversity order:

%Creating diagonal X matrix:
X = zeros(64,4);

for i=1:16
   X(4*i-3:4*i,:) = diag(symbolsRef(i,:)); 
end

%Finding minimum rank:
Lmin = inf;
for i=1:4:61
   for j=i+5:4:61 
      L = rank(X(i:i+3,:) - X(j:j+3,:));
      
      if L < Lmin
         Lmin = L; 
      end
   end
end

fprintf('Diversity order from performance curve is %d, while from system analysis is %d. \n ',LPlot, Lmin); 







