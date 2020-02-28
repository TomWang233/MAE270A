clc
close all
clear all
load u_rand.mat
load u1_impulse.mat
load u2_impulse.mat

y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs

y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);
ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets
t = [0:N-1]*ts - 1;
%%
%Hankel
co=['b*';'r*';'g*';'k*'];
N=[20 40 80 100];
for j=1:4
    n=N(j);
    row=[];
    Hankel=[];
    for k=0:n-1
        for i=1:n
            row=[row,[y11(mi+i+k) y12(mi+i+k);y21(mi+i+k) y22(mi+i+k)]];
            if i==n
                Hankel=[Hankel;row];
                row=[];
            end  
        end
    end
    s=sort(svd(Hankel),'descend');
    figure(3),
    semilogy(s,co(j,:),'LineWidth',2) 
    grid on;
    hold on;
    xlim([0 40]);
    ylim([1e-3 1]);
end
legend('H_2_0','H_4_0','H_8_0','H_1_0_0');
%% %Hankel2
row=[];
Hankel2=[];
n=100;
for k=0:n-1
    for i=1:n
        row=[row,[y11(mi+i+k+1) y12(mi+i+k+1);y21(mi+i+k+1) y22(mi+i+k+1)]];
        if i==n
           Hankel2=[Hankel2;row];
           row=[];
        end  
    end
end
%%  set dimension
% ns = 6;
% ns = 7;
% ns = 10;
 ns = 20;
[u,s,v]=svd(Hankel);
u1p = u(:,1:ns);
v1p = v(:,1:ns);
s1p = s(1:ns,1:ns);
On = u1p*s1p;
Cn = v1p';

A = pinv(On)*Hankel2*pinv(Cn);
disp('max(abs(eig(A)))');
max(abs(eig(A)))
B = Cn(:,1:2);
C = On(1:2,:);
D = zeros(2,2);
sys1 = ss(A,B,C,D,'Ts',ts);

[yp,tp] = impulse(sys1,2);
yp = yp*ts;
y11p = yp(:,1,1);
y21p = yp(:,2,1);
y12p = yp(:,1,2);
y22p = yp(:,2,2);
tt=['n_s=',num2str(ns)];

figure(4),
subplot(2,2,1)
plot(tp(2:end),y11p(2:end),'r*',t(mi+1:end),y11(mi+1:end),'b*','LineWidth',2)
grid on;
ylabel('y_1_1 (volts)');
axis([0 2 -0.1 0.1]);
title(tt);

subplot(2,2,3)
plot(tp(2:end),y21p(2:end),'r*',t(mi+1:end),y21(mi+1:end),'b*','LineWidth',2)
grid on;
ylabel('y_2_1(volts)');
axis([0 2 -0.1 0.1]);
title(tt);

subplot(2,2,2)
plot(tp(2:end),y12p(2:end),'r*',t(mi+1:end),y12(mi+1:end),'b*','LineWidth',2)
grid on;
ylabel('y_1_2(volts)');
axis([0 2 -0.1 0.1]);
title(tt);

subplot(2,2,4)
plot(tp(2:end),y22p(2:end),'r*',t(mi+1:end),y22(mi+1:end),'b*','LineWidth',2)
grid on;
ylabel('y_2_2(volts)');
axis([0 2 -0.1 0.1]);
title(tt);

%%   model frequency response
w0=[0:0.01:20]*2*pi;
I=eye(size(A));
MG(1:2,1:2,1:length(w0))=0;
AG(1:2,1:2,1:length(w0))=0;
for k=1:length(w0)
    w=w0(k);
    G0=C*inv(exp(1i.*ts.*w)*I-A)*B+D;
    MG(:,:,k)=abs(G0);   
    AG(:,:,k)=angle(G0)*180/pi;
end

figure(5);
subplot(2,2,1)
plot(w0,reshape(MG(1,1,1:end),1,length(w0)));
grid on;
y11pp=y11(mi:end);
L=length(y11pp);
y11f=fft(y11pp,L);
y11f(1)=y11f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,abs(y11f(1:end/2)),'r');
legend('y_1_1 Orig-sys','y_1_1 Ident-sys');

subplot(2,2,3)
plot(w0,reshape(MG(2,1,1:end),1,length(w0)))
grid on;
y21pp=y21(mi:end);
L=length(y21pp);
y21f=fft(y21pp,L);
y21f(1)=y21f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,abs(y21f(1:end/2)),'r');
legend('y_2_1 Orig-sys','y_2_1 Ident-sys');

subplot(2,2,2)
plot(w0,reshape(MG(1,2,1:end),1,length(w0)))
grid on;
y12pp=y12(mi:end);
L=length(y12pp);
y12f=fft(y12pp,L);
y12f(1)=y12f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,abs(y12f(1:end/2)),'r');
legend('y_1_2 Orig-sys','y_1_2 Ident-sys');

subplot(2,2,4)
plot(w0,reshape(MG(2,2,1:end),1,length(w0)))
grid on;
y22pp=y22(mi:end);
L=length(y22pp);
y22f=fft(y22pp,L);
y22f(1)=y22f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,abs(y22f(1:end/2)),'r');
legend('y_2_2 Orig-sys','y_2_2 Ident-sys');
%% 
figure(6);
subplot(2,2,1)
plot(w0,reshape(AG(1,1,1:end),1,length(w0)))
grid on;
y11pp=y11(mi:end);
L=length(y11pp);
y11f=fft(y11pp,L);
y11f(1)=y11f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,angle(y11f(1:end/2))*180/pi,'r');
legend('y_1_1 Orig-sys','y_1_1 Ident-sys');

subplot(2,2,3)
plot(w0,reshape(AG(2,1,1:end),1,length(w0)))
grid on;
y21pp=y21(mi:end);
L=length(y21pp);
y21f=fft(y21pp,L);
y21f(1)=y21f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,angle(y21f(1:end/2))*180/pi,'r');
legend('y_2_1 Orig-sys','y_2_1 Ident-sys');

subplot(2,2,2)
plot(w0,reshape(AG(1,2,1:end),1,length(w0)))
grid on;
y12pp=y12(mi:end);
L=length(y12pp);
y12f=fft(y12pp,L);
y12f(1)=y12f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(f(1:end/2)*2*pi,angle(y12f(1:end/2))*180/pi,'r');
legend('y_1_2 Orig-sys','y_1_2 Ident-sys');

subplot(2,2,4)
plot(w0,reshape(AG(2,2,1:end),1,length(w0)))
grid on;
y22pp=y22(mi:end);
L=length(y22pp);
y22f=fft(y22pp,L);
y22f(1)=y22f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Phase');
xlabel('w(rad/s)');
plot(f(1:end/2)*2*pi,angle(y22f(1:end/2))*180/pi,'r');
legend('y_2_2 Orig-sys','y_2_2 Ident-sys');
%% 4
y11pp=y11(mi:end);
y21pp=y21(mi:end);
y12pp=y12(mi:end);
y22pp=y22(mi:end);
u1pp=u1(mi:end);
u2pp=u2(mi:end);
y11f = fft(y11pp)./fft(u1pp);
y21f = fft(y21pp)./fft(u1pp);
y12f = fft(y12pp)./fft(u2pp);
y22f = fft(y22pp)./fft(u2pp);
N = length(y11f);
om = [0:N-1]/(ts*N); %%%% frequency vector in hertz
MGp11=abs(y11f);
AGp11=angle(y11f)*180/pi;
MGp21=abs(y21f);
AGp21=angle(y21f)*180/pi;
MGp12=abs(y12f);
AGp12=angle(y12f)*180/pi;
MGp22=abs(y22f);
AGp22=angle(y22f)*180/pi;

figure(7);
subplot(2,2,1)
plot(w0,reshape(MG(1,1,1:end),1,length(w0)));
grid on;
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),MGp11(1:end/2),'r');
legend('y_1_1 Orig-sys','y_1_1 Ident-sys');

subplot(2,2,3)
plot(w0,reshape(MG(2,1,1:end),1,length(w0)))
grid on;
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),MGp21(1:end/2),'r');
legend('y_2_1 Orig-sys','y_2_1 Ident-sys');

subplot(2,2,2)
plot(w0,reshape(MG(1,2,1:end),1,length(w0)))
grid on;
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),MGp12(1:end/2),'r');
legend('y_1_2 Orig-sys','y_1_2 Ident-sys');

subplot(2,2,4)
plot(w0,reshape(MG(2,2,1:end),1,length(w0)))
grid on;
hold on;
ylabel('Amplitude');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),MGp22(1:end/2),'r');
legend('y_2_2 Orig-sys','y_2_2 Ident-sys');
%% 
figure(8);
subplot(2,2,1)
plot(w0,reshape(AG(1,1,1:end),1,length(w0)))
grid on;
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),AGp11(1:end/2),'r');
legend('y_1_1 Orig-sys','y_1_1 Ident-sys');

subplot(2,2,3)
plot(w0,reshape(AG(2,1,1:end),1,length(w0)))
grid on;
y21pp=y21(mi:end);
L=length(y21pp);
y21f=fft(y21pp,L);
y21f(1)=y21f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),AGp21(1:end/2),'r');
legend('y_2_1 Orig-sys','y_2_1 Ident-sys');

subplot(2,2,2)
plot(w0,reshape(AG(1,2,1:end),1,length(w0)))
grid on;
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),AGp11(1:end/2),'r');
legend('y_1_2 Orig-sys','y_1_2 Ident-sys');

subplot(2,2,4)
plot(w0,reshape(AG(2,2,1:end),1,length(w0)))
grid on;
y22pp=y22(mi:end);
L=length(y22pp);
y22f=fft(y22pp,L);
y22f(1)=y22f(1);
f=1/ts/L*(0:(L-1));
hold on;
ylabel('Phase');
xlabel('w(rad/s)')
plot(2*pi*om(1:end/2),AGp11(1:end/2),'r');
legend('y_2_2 Orig-sys','y_2_2 Ident-sys');