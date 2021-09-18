clc
clear
close all

n = zeros(1000,1);
R = zeros(1000,1);

for i=1:1000
   n(i,1) = i;
   R(i,1) = 1000/(power(10,2.163/i)+1);
end

plot(n, R);
title('R vs. n');
xlabel('n');
ylabel('R');