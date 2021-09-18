clc
clear
close all

pOut = zeros(100,1);

for i = 1:100
   pOut(i,1) = i/100; 
end

pT = 6 * qfuncinv(pOut);

plot(pOut, pT);
title('P_o_u_t vs. needed TX power(in dbm)');
xlabel('P_o_u_t');
ylabel('P_T(dbm)');