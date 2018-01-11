function KeChengSheJi
%no noise
%Phase evenly distributed;
Freq=10;
Ts=0.01;
t=0:Ts:5;
A=5;
Ph=0;
Smpl=100;
%X(t)
for k=1: Smpl     
    X(k,:)=A*cos(Freq*2*pi*t+rand()*2*pi);%Phase random cosine signal 
end
set(figure(1),'name','Picture 1 of Jiajia¡®s Random Signal Analysis course design','Numbertitle','off');
%'Numbertitle','off'
subplot(3,1,1)
plot(X);
title('A cosine signal with evenly distributed phases')
[Tao,R]=XiangGuanHanshu(Ts,X);
subplot(3,1,2)
plot(Tao,R);
title('Autocorrelation function of X (t)')
[f,Y_F]=GongLvPu(R,Ts);
subplot(3,1,3)
plot(f,abs(Y_F)) ;
title('Power spectrum of X (t)')
xlabel('Freq(Hz)');
ylabel('|X(f)|')

%%%%%%%%%Gaussian distribution of white noise
XN_t=cos(Freq*2*pi*t);
for k=1: Smpl     
    N(k,:)=randn(1,length(t)).*XN_t;    
end
set(figure(2),'name','Picture 2 of Jiajia¡®s Random Signal Analysis course design','Numbertitle','off');
subplot(3,1,1)
plot(N);
title('Gaussian distribution of white noise')
[Tao,R]=XiangGuanHanshu(Ts,N);
%R(find(R==max(R)))=100000;
subplot(3,1,2)
plot(Tao,R);
title('Autocorrelation Function of Gaussian Distributed White Noise N (t)')
[f,Y_F]=GongLvPu(R,Ts);
subplot(3,1,3)
plot(f,abs(Y_F)) ;
title('Power Spectrum of Gaussian Distributed White Noise N (t)')
xlabel('Freq(Hz)');
ylabel('|X(f)|')


%%%%%Y(t)=X(t)+N(t)
Y=X+N;
set(figure(3),'name','Picture 2 of Jiajia¡®s Random Signal Analysis course design','Numbertitle','off');
subplot(3,1,1)
plot(Y);
title('Random signal plus Gaussian white noise')
subplot(3,1,2)
[Tao,R]=XiangGuanHanshu(Ts,Y);
plot(Tao,R);
title('Stochastic signal plus Gaussian white noise after the autocorrelation function')
subplot(3,1,3)
[f,Y_F]=GongLvPu(R,Ts);
plot(f,abs(Y_F)) ;
title('Random signal plus Gaussian white noise after the power spectrum')
xlabel('Freq (Hz)');
ylabel('|X(f)|');


function [f,Y_F]=GongLvPu(R,Ts)
L=length(R);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y_F_R = fft(R,NFFT)/L;
f_R =(1/Ts)/2*linspace(0,1,NFFT/2+1);
% Plot single-sided amplitude spectrum.
f_L=-1*fliplr(f_R);
Y_F_L=fliplr(Y_F_R(1:NFFT/2+1));
f=[f_L,f_R];
Y_F=[Y_F_L Y_F_R(1:NFFT/2+1)];



function [Tao,R_Tao_PingJun]=XiangGuanHanshu(Ts,Y)
%Y[Sample*time]
%Ts:Sample time
[L_smp,t_smp]=size(Y);
RM=zeros(t_smp,t_smp);
for k=1:L_smp
    RM=RM+Y'*Y;
end
RM=RM/L_smp/L_smp;
v = diag(RM,1);
Pos=1;
M=20;
Tao_M=[];
R_Tao_M=[];
for k=-(t_smp-1):t_smp-1
    v = diag(RM,k);
    R_Tao_PingJun(Pos)=sum(v)/length(v);
    Tao(Pos)=k*Ts;
    Pos=Pos+1;
    if(length(v)>M)
        Tao_M=[Tao_M k*Ts];
        R_Tao_M=[R_Tao_M v(M)];
    end
end