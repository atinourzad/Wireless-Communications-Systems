%%%HW6_Wireless_Ali_Shokouhifar_810198186
%q5
% I
clear all
close all
clc

tic
N=10^4;
sigma=sqrt(1/2);
H=sigma.*(randn([4*N 4])+1i*randn([4*N 4])); %h(channel response) CN variable
Ph=mean(abs(H).^2) % expected value of abs(h)^2 is 1

xb=decimalToBinaryVector(0:15);
xs=[0 0 0 0;0 1 1 2;0 2 3 1;0 3 2 3;1 0 2 1;1 1 3 3;1 2 1 0;1 3 0 2;2 0 1 3;2 1 0 1;2 2 2 2;2 3 3 0;3 0 3 2;3 1 2 0;3 2 0 3;3 3 1 1];

x1=zeros(N,4);
x1(find(rand([N 4])>0.5))=1; %x1 : bits sequence

x2=zeros(N,4);%x2 : symbols sequence
for i=1:N
    for j=1:16
        if x1(i,:)==xb(j,:)
            x2(i,:)=xs(j,:);
        end
    end
end

m=zeros(N,4);
m=exp(1j*pi*x2/2);% m : symbols sequence after modulation QPSK
X=m.';

xsm=exp(1j*pi*xs/2);%xsm : symbols after QPSK modulation 
Px=std(xsm).^2 %Power of modulated symbols is 1

sigma=sqrt(1/2);
n=sigma.*(randn([4 N])+1i*randn([4 N])); %n(noise) CN variable
% pn=mean(abs(n).^2,'all') %Power of noise is 1
nums=10;
first=0;
for k=first:nums
SNR=10.^(k/10);
N02=sqrt(1/(SNR)); % N02 is (2N0)^0.5 that is power of CN noise
n=n*N02;

Y=zeros(4,N);
for i=1:N
    Y(:,i) = H(4*(i-1)+1:4*i,:)*X(:,i) + n(:,i); % AWGN & fading channel
end

% detection using Maximum-Likelihood (ML) decoder
Xml=zeros(N,4);
for i=1:N
    for j=1:16
        arg(j)=sum(abs(Y(:,i)-H(4*(i-1)+1:4*i,:)*xsm(j,:).').^2);
    end
    [M I]=min(arg);
    Xml(i,:)=xs(I,:);
end

% corrects=0;
% for i=1:N
%     if Xml(i,:)==x2(i,:)
%         corrects=corrects+1;
%     end
% end
% errors=N-corrects;
% Pe1(k-first+1)=errors/N;

corrects=0;
for i=1:N
    for j=1:4
    if Xml(i,j)==x2(i,j)
        corrects=corrects+1;
    end
    end
end
errors=4*N-corrects;
Pe1(k-first+1)=errors/(4*N);

end
toc

y1=10^(-5);
figure
semilogy(first:nums,Pe1((1:nums-first+1)))
xlim([first nums]);ylim([y1 1]);
xlabel('SNR in dB');ylabel('Pe');
title(['Pe in part I']);

slope=(10*log10(Pe1(5))-10*log10(Pe1(4)));

%%
% II

clc

tic
N=10^4
sigma=sqrt(1/2);
H=sigma.*(randn([4*N 4])+1i*randn([4*N 4])); %h(channel response) CN variable
Ph=mean(abs(H).^2) % expected value of abs(h)^2 is 1 

xb=decimalToBinaryVector(0:15);
xs=[0 0 0 0;0 1 1 2;0 2 3 1;0 3 2 3;1 0 2 1;1 1 3 3;1 2 1 0;1 3 0 2;2 0 1 3;2 1 0 1;2 2 2 2;2 3 3 0;3 0 3 2;3 1 2 0;3 2 0 3;3 3 1 1];

x1=zeros(N,4);
x1(find(rand([N 4])>0.5))=1; %x1 : bits sequence

x2=zeros(N,4);%x2 : symbols sequence
for i=1:N
    for j=1:16
        if x1(i,:)==xb(j,:)
            x2(i,:)=xs(j,:);
        end
    end
end

m=zeros(N,4);
m=exp(1j*pi*x2/2);% m : symbols sequence after modulation QPSK
X1=m.';

xsm=exp(1j*pi*xs/2);%xsm : symbols after QPSK modulation 
Px=std(xsm).^2 %Power of modulated symbols is 1

sigma=sqrt(1/2);
n=sigma.*(randn([4 2*N])+1i*randn([4 2*N])); %n(noise) CN variable
pn=mean(abs(n).^2,'all') %Power of noise is 1
nums=10;
first=0;
for k=first:nums
SNR=10.^(k/10);
N02=sqrt(1/(SNR)); % N02 is (2N0)^0.5 that is power of CN noise
n=n*N02;

Y=zeros(4,2*N);
X=zeros(4,2);
Xml=zeros(N,4);
x=zeros(4,2);
for i=1:N
    X(:,1)=X1(:,i);
    X(:,2)=[-conj(X1(2,i));conj(X1(1,i));-conj(X1(4,i));conj(X1(3,i))];
    Y(:,2*(i-1)+1:2*i) = H(4*(i-1)+1:4*i,:)*X + n(:,2*(i-1)+1:2*i); % AWGN & fading channel
end


% detection using Maximum-Likelihood (ML) decoder
Xml=zeros(N,4);
x=zeros(4,2);
for i=1:N
    for j=1:16
        x(:,1)=xsm(j,:);
        x(:,2)=[-conj(xsm(j,2));conj(xsm(j,1));-conj(xsm(j,4));conj(xsm(j,3))];
        arg(j)=sum(abs(Y(:,2*(i-1)+1:2*i)-H(4*(i-1)+1:4*i,:)*x).^2,'all');
    end
    [M I]=min(arg);
    Xml(i,:)=xs(I,:);
end

% corrects=0;
% for i=1:N
%     if Xml(i,:)==x2(i,:)
%         corrects=corrects+1;
%     end
% end
% errors=N-corrects;
% Pe2(k-first+1)=errors/N;

corrects=0;
for i=1:N
    for j=1:4
    if Xml(i,j)==x2(i,j)
        corrects=corrects+1;
    end
    end
end
errors=4*N-corrects;
Pe2(k-first+1)=errors/(4*N);

end
toc

y1=10^(-5);
figure
semilogy(first:nums,Pe2((1:nums-first+1)))
xlim([first nums]);ylim([y1 1]);
xlabel('SNR in dB');ylabel('Pe');
title(['Pe in part II']);

%%
% III (type 1)

clc

tic
N=10^5
sigma=sqrt(1/2);
H=sigma.*(randn([4*N 4])+1i*randn([4*N 4])); %h(channel response) CN variable
Ph=mean(abs(H).^2) % expected value of abs(h)^2 is 1

% SVD of H matrix
S_matrix=zeros(4*N,4);
U_matrix=zeros(4*N,4);
V_matrix=zeros(4*N,4);
for i=1:N
    [U,S,V] = svd(H(4*(i-1)+1:4*i,:));
    S_matrix(4*(i-1)+1:4*i,:)=S;
    U_matrix(4*(i-1)+1:4*i,:)=U;
    V_matrix(4*(i-1)+1:4*i,:)=V;
end


xb=decimalToBinaryVector(0:15);
% xs=[0 0 0 0;0 1 1 2;0 2 3 1;0 3 2 3;1 0 2 1;1 1 3 3;1 2 1 0;1 3 0 2;2 0 1 3;2 1 0 1;2 2 2 2;2 3 3 0;3 0 3 2;3 1 2 0;3 2 0 3;3 3 1 1];

x1=zeros(N,4);
x1(find(rand([N 4])>0.5))=1; %x1 : bits sequence

% x2=zeros(N,4);%x2 : symbols sequence
% for i=1:N
%     for j=1:16
%         if x1(i,:)==xb(j,:)
%             x2(i,:)=xs(j,:);
%         end
%     end
% end

QAM16=[1+1j 1+3j 3+1j 3+3j -1+1j -1+3j -3+1j -3+3j 1-1j 1-3j 3-1j 3-3j -1-1j -1-3j -3-1j -3-3j];
QAM16=QAM16/(std(QAM16));
m=zeros(N,1);
for i=1:N
    for j=1:16
        if x1(i,:)==xb(j,:)
            m(i,1)=QAM16(j);% m : symbols sequence after modulation QPSK
        end
    end
end

Px=std(QAM16).^2 %Power of modulated symbols is 1

sigma=sqrt(1/2);
n=sigma.*(randn([4 N])+1i*randn([4 N])); %n(noise) CN variable
% Pn=mean(abs(n).^2,'all') %Power of noise is 1
nums=10;
first=0;
for k=first:nums
SNR=10.^(k/10);
N02=sqrt(1/(SNR)); % N02 is (2N0)^0.5 that is power of CN noise

% r1= m + n; % only AWGN

n=n*N02;
Y_tilda=zeros(4,N);
for i=1:N
    S=S_matrix(4*(i-1)+1:4*i,:);
    U=U_matrix(4*(i-1)+1:4*i,:);
    X_tilda=[m(i);0;0;0];
    n_tilda=(U')*n(:,i);
    Y_tilda(:,i) = S*X_tilda + n_tilda; % AWGN & fading channel
end

% detection using Maximum-Likelihood (ML) decoder
Xml=zeros(N,1);
for i=1:N
    for j=1:16
        S=S_matrix(4*(i-1)+1:4*i,:);
        x_tilda=[QAM16(j);0;0;0];
        arg(j)=sum(abs(Y_tilda(:,i)-S*x_tilda).^2);
    end
    [M I]=min(arg);
    Xml(i)=QAM16(I);
end

corrects=0;
for i=1:N
    if Xml(i)==m(i)
        corrects=corrects+1;
    end
end
errors=N-corrects;
Pe3(k-first+1)=errors/N;

end
toc

y1=10^(-5);
figure
semilogy(first:nums,Pe3((1:nums-first+1)))
xlim([first nums]);ylim([y1 1]);
xlabel('SNR in dB');ylabel('Pe');
title(['Pe in part III (type 1)']);

% figure
% semilogy(first:nums,Pe1((1:nums-first+1)))
% xlim([first nums]);ylim([y1 1]);
% hold on
% semilogy(first:nums,Pe2((1:nums-first+1)))
% xlim([first nums]);ylim([y1 1]);
% xlabel('SNR in dB');ylabel('Pe');
% hold on
% semilogy(first:nums,Pe3((1:nums-first+1)))
% xlim([first nums]);ylim([y1 1]);
% xlabel('SNR in dB');ylabel('Pe');
% title(['Pe in part I & II & III (type 1)']);

%%
% III (type 2)

clc
clear
tic
N=10^4
sigma=sqrt(1/2);
H=sigma.*(randn([4*N 4])+1i*randn([4*N 4])); %h(channel response) CN variable
Ph=mean(abs(H).^2) % expected value of abs(h)^2 is 1

% SVD of H matrix
S_matrix=zeros(4*N,4);
U_matrix=zeros(4*N,4);
V_matrix=zeros(4*N,4);
for i=1:N
    [U,S,V] = svd(H(4*(i-1)+1:4*i,:));
    S_matrix(4*(i-1)+1:4*i,:)=S;
    U_matrix(4*(i-1)+1:4*i,:)=U;
    V_matrix(4*(i-1)+1:4*i,:)=V;
end

xb=decimalToBinaryVector(0:15);

x1=zeros(N,4);
x1(find(rand([N 4])>0.5))=1; %x1 : bits sequence


QAM16=[1+1j 1+3j 3+1j 3+3j -1+1j -1+3j -3+1j -3+3j 1-1j 1-3j 3-1j 3-3j -1-1j -1-3j -3-1j -3-3j];
QAM16=QAM16/(std(QAM16));
m=zeros(N,1);
for i=1:N
    for j=1:16
        if x1(i,:)==xb(j,:)
            m(i,1)=QAM16(j);% m : symbols sequence after modulation 16-QAM
        end
    end
end

Px=std(QAM16).^2 %Power of modulated symbols is 1

sigma=sqrt(1/2);
n=sigma.*(randn([4 N])+1i*randn([4 N])); %n(noise) CN variable
% Pn=mean(abs(n).^2,'all') %Power of noise is 1
nums=10;
first=0;
for k=first:nums
SNR=10.^(k/10);
N02=sqrt(1/(SNR)); % N02 is (2N0)^0.5 that is power of CN noise

% r1= m + n; % only AWGN

n=n*N02;
Y_tilda=zeros(1,N);
for i=1:N
    S=S_matrix(4*(i-1)+1:4*i,:);
    U=U_matrix(4*(i-1)+1:4*i,:);
    X_tilda=m(i);
    n_tilda=(U')*n(:,i);
    Y_tilda(1,i) = S(1,1)*X_tilda + n_tilda(1); % AWGN & fading channel
end

% detection using Maximum-Likelihood (ML) decoder
Xml=zeros(N,1);
for i=1:N
    for j=1:16
         S=S_matrix(4*(i-1)+1:4*i,:);
         x_tilda=QAM16(j);
        arg(j)=sum(abs(Y_tilda(1,i)-S(1,1)*x_tilda).^2);
    end
    [M I]=min(arg);
    Xml(i)=QAM16(I);
end

corrects=0;
for i=1:N
    if Xml(i)==m(i)
        corrects=corrects+1;
    end
end
errors=N-corrects;
Pe33(k-first+1)=errors/N;

end
toc

y1=10^(-5);
figure
semilogy(first:nums,Pe33((1:nums-first+1))/4)
xlim([first nums]);ylim([y1 1]);
xlabel('SNR in dB');ylabel('Pe');
title(['Pe in part III (type 2)']);


slope=(10*log10(Pe33(5)/4)-10*log10(Pe33(4)/4))
slope_final=floor(mean(slope));
L1=abs(slope_final) % diversity Order


