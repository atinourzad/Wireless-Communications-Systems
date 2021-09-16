%% Part a

clc
clear
close all

numOfBits = 10^4;
bitsRef = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1; 1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1; 1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1];
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

%Modulating refrence symbols:
modulatedSymbolsRef = exp(1i*pi*symbolsRef/2);

%Creating fading gains:
hMean = 0;
hSD = sqrt(1/2);
h = normrnd(hMean, hSD, [4*numOfBits,4]) + 1i * normrnd(hMean, hSD, [4*numOfBits,4]);

%Modulated symbols go through fading gain:
fadedSymbols = zeros(numOfBits,4);
for i=1:numOfBits
    fadedSymbols(i,:) =  modulatedSymbols(i,:) * h(4*i-3:4*i,:).';
end

%%%-----------------------------------------------Finding Probability of Error----------------------------------------%%%

%Finding Different SNRs:
SNRdB = 0:20;
SNR = 10 .^ (SNRdB/10);

pBitErrorFadedA = zeros(length(SNR), 1);
pSymbolErrorFaded = zeros(length(SNR), 1);

for i=1:length(SNR)
    
    %Creating AWGN noise:
    noiseMean = 0;
    noiseSD = sqrt(1/(1*SNR(i)));
    noise = normrnd(noiseMean, noiseSD, [numOfBits,4]) + 1i * normrnd(noiseMean, noiseSD, [numOfBits,4]);
    
    %Faded signals go through AWGN channel:
    fadedY = fadedSymbols + noise;
    
    %Using ML
    detectedFadedBits = zeros(numOfBits,4);
    detectedFadedSymbols = zeros(numOfBits,4);
    
    for j=1:numOfBits
      minDistFaded = inf;
      
      for k=1:16
          %Modulated refreces go through fading:
          fadedRef = modulatedSymbolsRef(k,:) * h(4*j-3:4*j,:).';
          
          %Finding argmin:
          distFaded = norm(fadedY(j,:) - fadedRef);
          
          if distFaded < minDistFaded
              minDistFaded = distFaded;
              fadedSymbolIndex = k;
          end
      end
      detectedFadedBits(j,:) = bitsRef(fadedSymbolIndex,:);
      detectedFadedSymbols(j,:) = symbolsRef(fadedSymbolIndex,:);
    end 
    pBitErrorFadedA(i,1) = nnz(bits - detectedFadedBits)/numOfBits;
    pSymbolErrorFaded(i,1) = nnz(codedSymbols - detectedFadedSymbols)/numOfBits;
end

figure 
semilogy(SNRdB, (pBitErrorFadedA), 'Color', [19/255, 206/255, 188/255]);
hold on
semilogy(SNRdB, (pSymbolErrorFaded), 'Color', [28/255, 152/255, 140/255]);
hold off
xlabel('SNR (dB)');
ylabel('P_e');
xlim([0 20]);
legend('BER', 'SER');
title("Probability of Error vs. Different SNR values for Given Coding");

%Diversity Order:
for i=length(pSymbolErrorFaded):-1:1
   if pSymbolErrorFaded(i) ~= 0 && pSymbolErrorFaded(i-1) ~= 0 && pSymbolErrorFaded(i) ~= pSymbolErrorFaded(i-1)
      index = i;
      break
   end
end
%%
LPlot = (abs(10*log10(pSymbolErrorFaded(index)) - 10*log10(pSymbolErrorFaded(index-1))));
fprintf('Diversity order from performance curve is %f \n ',LPlot); 

save('BER_PartA.mat', 'pBitErrorFadedA');
save('Bits.mat', 'bits');
save('Channel Matrix.mat', 'h');

%% Part b:

clc
clear 
close all

numOfBits = 10^4;

%Loading Created Bits:
load('Bits.mat');

%Creating Refrence Bits:
totalBits = decimalToBinaryVector(0:255);
symbolsRef = zeros(256,4);

for i=1:256
   for j=1:2:8
       symbolsRef(i,ceil(j/2)) = bi2de(totalBits(i,j:j+1), 'left-msb');
   end
end

%Creating coded symbols:
codedSymbols = zeros(numOfBits, 4);

for i=1:numOfBits
    bitBlock = bi2de(bits(i,:), 'left-msb');
    codedSymbols(i,:) = symbolsRef(bitBlock+1,:);
end

%Creating modulated symbols with QPSK mudolation: 0->1, 1->j, 2->-1, 3->-j 
modulatedSymbols = exp(1i*pi*codedSymbols/2);

%Modulating refrence symbols:
modulatedSymbolsRef1 = exp(1i*pi*symbolsRef/2);
modulatedSymbolsRef2 = [-conj(modulatedSymbolsRef1(:,2)) conj(modulatedSymbolsRef1(:,1)) -conj(modulatedSymbolsRef1(:,4)) conj(modulatedSymbolsRef1(:,3))]; 

%Alamouti Sent Symbols:
modulatedSentSymbol1 = modulatedSymbols;
modulatedSentSymbol2 = [-conj(modulatedSymbols(:,2)) conj(modulatedSymbols(:,1)) -conj(modulatedSymbols(:,4)) conj(modulatedSymbols(:,3))]; 

%Loading fading gains:
load('Channel Matrix.mat');

%Modulated symbols go through fading gain:
fadedSymbols1 = zeros(numOfBits,4);
fadedSymbols2 = zeros(numOfBits,4);

for j=1:numOfBits
    fadedSymbols1(j,:) = modulatedSentSymbol1(j,:) * h(4*j-3:4*j,:).';
    fadedSymbols2(j,:) = modulatedSentSymbol2(j,:) * h(4*j-3:4*j,:).';
end
    

%%%-----------------------------------------------Finding Probability of Error----------------------------------------%%%

%Finding Different SNRs:
SNRdB = 0:20;
SNR = 10 .^ (SNRdB/10);

pSymbolErrorFaded = zeros(length(SNR), 1);
for i=1:length(SNR)
    
    %Creating AWGN noise:
    noiseMean = 0;
    noiseSD = sqrt(1/(1*SNR(i)));
    noise1 = normrnd(noiseMean, noiseSD, [numOfBits,4]) + 1i * normrnd(noiseMean, noiseSD, [numOfBits,4]);
    noise2 = normrnd(noiseMean, noiseSD, [numOfBits,4]) + 1i * normrnd(noiseMean, noiseSD, [numOfBits,4]);
    
    
    %Faded signals go through AWGN channel:
    fadedY1 = fadedSymbols1 + noise1;
    fadedY2 = fadedSymbols2 + noise2;
    
    %Using ML
    detectedFadedSymbols = zeros(numOfBits,4);
    
    for j=1:numOfBits
      minDistFaded = inf;
      
      for k=1:256
          %Modulated refreces go through fading:
          fadedRef1 = modulatedSymbolsRef1(k,:) * h(4*j-3:4*j,:).';
          fadedRef2 = modulatedSymbolsRef2(k,:) * h(4*j-3:4*j,:).';
          
          %Finding argmin:
          distFaded = sqrt(norm(fadedY1(j,:) - fadedRef1)^2 + norm(fadedY2(j,:) - fadedRef2)^2);
          
          
          if distFaded < minDistFaded
              minDistFaded = distFaded;
              fadedSymbolIndex = k;
          end
      end
      detectedFadedSymbols(j,:) = symbolsRef(fadedSymbolIndex,:);
    end 
    pSymbolErrorFaded(i,1) = nnz(codedSymbols - detectedFadedSymbols)/numOfBits;
end

figure 
semilogy(SNRdB, pSymbolErrorFaded/2, 'Color', [19/255, 206/255, 188/255]);
xlabel('SNR (dB)');
ylabel('BER');
xlim([0 20]);
title("Probability of Error Bit vs. SNR for Alamouti");

%Diversity Order:

for i=length(pSymbolErrorFaded):-1:1
   if pSymbolErrorFaded(i) ~= 0 && pSymbolErrorFaded(i-1) ~= 0 && pSymbolErrorFaded(i) ~= pSymbolErrorFaded(i-1)
      index = i;
      break
   end
end

LPlot = (abs(10*log10(pSymbolErrorFaded(index)) - 10*log10(pSymbolErrorFaded(index-1))));
fprintf('Diversity order from performance curve is %f \n ',LPlot); 

pBitErrorFadedB = pSymbolErrorFaded/2;
save('BER_PartB.mat', 'pBitErrorFadedB');

%% Part c:


clc
clear 
close all

numOfBits = 10^4;
symbolsRef = [-3-3i; -3-i; -3+3i; -3+i; -1-3i; -1-i; 1+3i; 1+i; 3-3i; 3-i; 3+3i; 3+i; 1-3i; 1-i; 1+3i; 1+i];

%Loading Created Bits:
load('Bits.mat');

%Created bits are modulated 16-QAM mudolation:
modulatedSymbols = zeros(numOfBits, 1);

for i=1:numOfBits
    bitBlock = bi2de(bits(i,:), 'left-msb');
    modulatedSymbols(i,1) = symbolsRef(bitBlock+1,1);
end

%Loading fading gains:
load('Channel Matrix.mat');

%SVD:
U = zeros(4*numOfBits,4);
sigmaMax = zeros(numOfBits,1);

for i=1:numOfBits
   [u, s, v] = svd(h(4*i-3:4*i,:));
   U(4*i-3:4*i,:) = u;
   sigmaMax(i,1) = s(1,1);
end

%%%-----------------------------------------------Finding Probability of Error----------------------------------------%%%

%Finding Different SNRs:
SNRdB = 0:20;
SNR = 10 .^ (SNRdB/10);

pSymbolErrorFaded = zeros(length(SNR), 1);
for i=1:length(SNR)
    
    %Creating AWGN noise:
    noiseMean = 0;
    noiseSD = sqrt(1/(1*SNR(i)));
    noise = normrnd(noiseMean, noiseSD, [4,numOfBits]) + 1i * normrnd(noiseMean, noiseSD, [4,numOfBits]);
    
    %Modulated symbols go through fading gain and AWGN channel:
    fadedSymbols = zeros(numOfBits,1);
    for j=1:numOfBits
        noiseTilda = U(4*j-3:4*j,:)' * noise(:,j);
        fadedSymbols(j,:) =  modulatedSymbols(j,1) * sigmaMax(j) + noiseTilda(1);
    end
    
    %Using ML
    detectedFadedSymbols = zeros(numOfBits,1);
    
    for j=1:numOfBits
      minDistFaded = inf;
      
      for k=1:16
          %Modulated refreces go through fading:
          fadedRef = sigmaMax(j) * symbolsRef(k,1);
          
          %Finding argmin:
          distFaded = norm(fadedSymbols(j,1) - fadedRef);
          
          if distFaded < minDistFaded
              minDistFaded = distFaded;
              fadedSymbolIndex = k;
          end
          
      end
      detectedFadedSymbols(j,1) = symbolsRef(fadedSymbolIndex,1);
    end 
    pSymbolErrorFaded(i,1) = nnz(modulatedSymbols - detectedFadedSymbols)/numOfBits;
end

figure 
semilogy(SNRdB, pSymbolErrorFaded/4, 'Color', [19/255, 206/255, 188/255]);
xlabel('SNR (dB)');
ylabel('BER');
xlim([0 20]);
title("Probability of Error Bit vs. SNR Using BeamForming in TX and RX");

for i=length(pSymbolErrorFaded):-1:1
   if pSymbolErrorFaded(i) ~= 0 && pSymbolErrorFaded(i-1) ~= 0 && pSymbolErrorFaded(i) ~= pSymbolErrorFaded(i-1)
      index = i;
      break
   end
end

LPlot = (abs(10*log10(pSymbolErrorFaded(index)) - 10*log10(pSymbolErrorFaded(index-1))));
fprintf('Diversity order from performance curve is %f \n ',LPlot); 

pBitErrorFadedC = pSymbolErrorFaded/4;
save('BER_PartC.mat', 'pBitErrorFadedC');

%%
load('BER_PartA.mat');
load('BER_PartB.mat');
load('BER_PartC.mat');

figure 
semilogy(SNRdB, pBitErrorFadedA);
hold on
semilogy(SNRdB, pBitErrorFadedB);
semilogy(SNRdB, pBitErrorFadedC);
xlabel('SNR (dB)');
ylabel('BER');
xlim([0 20]);
legend('BER Part a', 'BER Part b', 'BER Part c')
