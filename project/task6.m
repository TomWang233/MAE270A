clc
close all
clear all
%%task6
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
om = [0:N-1]/(ts*N);
%%
%Hankel
row1 = [];
Hankel1 = [];
n = 100;
for k = 0:n-1
    for i = 1:n
        row1 = [row1,[y11(mi+i+k) y12(mi+i+k);y21(mi+i+k) y22(mi+i+k)]];
        if i == n
           Hankel1 =[Hankel1;row1];
           row1 = [];
        end  
    end
end
%% %Hankel2
row = [];
Hankel2 = [];
n = 100;
for k = 0:n-1
    for i = 1:n
        row = [row,[y11(mi+i+k+1) y12(mi+i+k+1);y21(mi+i+k+1) y22(mi+i+k+1)]];
        if i == n
           Hankel2 = [Hankel2;row];
           row = [];
        end  
    end
end
[u,s,v] = svd(Hankel1);
%%  compute system
ns = 7;
u1p = u(:,1:ns);
v1p = v(:,1:ns);
s1p = s(1:ns,1:ns);
On = u1p*s1p;
Cn = v1p';
A = pinv(On)*Hankel2*pinv(Cn);
B = Cn(:,1:2);
C = On(1:2,:);
D = zeros(2,2);
wl = 0;
wu = 125;
tol = 1e-15;
ts = 1/40;
%3
[Hinf,winf]=Inf_norm_discrete(A,B,C,D,wl,wu,tol,ts)
%4
u = 1;
for w = wl:0.01:wu
    Pw = C*inv((exp(1i*w*ts)*eye(length(A))-A))*B+D;
    S1(:,u) = svd(Pw);
    u = u+1;
end

for i = 1:1:181
    hk = [y11f(i) y12f(i);y21f(i) y22f(i)];
    S2(:,i) = svd(hk);
    i = i + 1;
end
w=wl:0.01:wu;
figure(1)
plot(w,S1(1,:),'--r');
grid on;
hold on;
plot(winf,Hinf,'*');
xlabel('Frequency');
ylabel('Magnitude');
hold on;
plot(2*pi*om(1:181),S2(1,:),'b');
legend('from model','H-inf','from data');
