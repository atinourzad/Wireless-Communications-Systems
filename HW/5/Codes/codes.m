%% Part b

clc
clear

options = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point');
[g,RBarMin] = fmincon(@rate, [3 7 14 20], [], [],[],[],[],[],@avgPowerB,options);

fprintf(' The resulting bounds for g are: \n')
fprintf(' %f \t ',g);

RBarMax = -1 * RBarMin;
fprintf('\n \n Average rate is %f. \n ',RBarMax);

deltaG = 0.001;

g1 = g(1):deltaG:g(2);
p1 = sum((10*log10(2)-g1) .* (1/sqrt(2*pi*64)) .* exp(-((g1-10).^2)/128).* deltaG);

g2 = g(2):deltaG:g(3);
p2 = sum((10*log10(6)-g2) .* (1/sqrt(2*pi*64)) .* exp(-((g2-10).^2)/128).* deltaG);

g3 = g(3):deltaG:g(4);
p3 = sum((10*log10(28)-g3) .* (1/sqrt(2*pi*64)) .* exp(-((g3-10).^2)/128).* deltaG);

g4 = g(4):deltaG:g(4)+10;
p4 = sum((10*log10(120)-g4) .* (1/sqrt(2*pi*64)) .* exp(-((g4-10).^2)/128) .* deltaG);

calcAvgPower = p1 + p2 + p3 + p4;
calcAvgPower = 10^(calcAvgPower/10);

fprintf('\n Average power is %f, showing that the constraint holds. \n ',calcAvgPower);
fprintf('\n Average power for region 1 with rate = 1 is %f. \n ',10^(p1/10));
fprintf('Average power for region 2 with rate = 2 is %f. \n ',10^(p2/10));
fprintf('Average power for region 3 with rate = 3 is %f. \n ',10^(p3/10));
fprintf('Average power for region 4 with rate = 4 is %f. \n ',10^(p4/10));

prob1 = qfunc((g(1)-10)/8) - qfunc((g(2)-10)/8);
prob2 = qfunc((g(2)-10)/8) - qfunc((g(3)-10)/8);
prob3 = qfunc((g(3)-10)/8) - qfunc((g(4)-10)/8);
prob4 = qfunc((g(4)-10)/8);

fprintf('\n Probability of using region 1 with rate = 1 is %f. \n ',prob1);
fprintf('Probability of using region 2 with rate = 2 is %f. \n ',prob2);
fprintf('Probability of using region 3 with rate = 3 is %f. \n ',prob3);
fprintf('Probability of using region 4 with rate = 4 is %f. \n ',prob4);


%% Part c
clc
clear

options = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point');
[g,RBarMin] = fmincon(@rate, [3 7 14 20], [], [],[],[],[],[],@avgPowerC,options);

fprintf(' The resulting bounds for g are: \n')
fprintf(' %f \t ',g);

RBarMax = -1 * RBarMin;
fprintf('\n \n Average rate is %f. \n ',RBarMax);

p1 = 2/(10^(g(1)/10));
p2 = 6/(10^(g(2)/10));
p3 = 28/(10^(g(3)/10));
p4 = 120/(10^(g(4)/10));
calcAvgPower = p1 * qfunc((g(1)-10)/8) + (p2-p1) * qfunc((g(2)-10)/8) + (p3-p2) * qfunc((g(3)-10)/8) + (p4-p3) * qfunc((g(4)-10)/8);

fprintf('\n Average power is %f, showing that the constraint holds. \n ',calcAvgPower);
fprintf('\n Power for region 1 with rate = 1 is %f. \n ',p1);
fprintf('Power for region 2 with rate = 2 is %f. \n ',p2);
fprintf('Power for region 3 with rate = 3 is %f. \n ',p3);
fprintf('Power for region 4 with rate = 4 is %f. \n ',p4);

prob1 = qfunc((g(1)-10)/8) - qfunc((g(2)-10)/8);
prob2 = qfunc((g(2)-10)/8) - qfunc((g(3)-10)/8);
prob3 = qfunc((g(3)-10)/8) - qfunc((g(4)-10)/8);
prob4 = qfunc((g(4)-10)/8);

fprintf('\n Probability of using region 1 with rate = 1 is %f. \n ',prob1);
fprintf('Probability of using region 2 with rate = 2 is %f. \n ',prob2);
fprintf('Probability of using region 3 with rate = 3 is %f. \n ',prob3);
fprintf('Probability of using region 4 with rate = 4 is %f. \n ',prob4);

