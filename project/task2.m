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
row1 =[];
Hankeln=[];
n=100;
for k=0:n-1
    for i=1:n
        row1 =[row1,[y11(mi+i+k) y12(mi+i+k);y21(mi+i+k) y22(mi+i+k)]];
        if i==n
           Hankeln =[Hankeln;row1];
           row1 =[];
        end  
    end
end
%% %Hankelp
row2 =[];
Hankelnp=[];
n=100;
for k=0:n-1
    for i=1:n
        row2 =[row2,[y11(mi+i+k+1) y12(mi+i+k+1);y21(mi+i+k+1) y22(mi+i+k+1)]];
        if i==n
           Hankelnp=[Hankelnp;row2];
           row2 =[];
        end  
    end
end
[u,s,v]=svd(Hankeln);
%%  set dimension
% ns=7;
ns=8;
u1p=u(:,1:ns);
v1p=v(:,1:ns);
s1p=s(1:ns,1:ns);
On=u1p*s1p;
Cn=v1p';
A=pinv(On)*Hankelnp*pinv(Cn);
B=Cn(:,1:2);
C=On(1:2,:);
D=zeros(2,2);
sys1=ss(A,B,C,D,'Ts',ts);


%%  task2
a = 0.1;
I = eye(size(A));
S = [a*I-A -B;C D];
display(['rank(S) = ', num2str(rank(S))]); 
AA=[A B;-C -D];
BB=[I zeros(size(B));zeros(size(C)) zeros(size(D))];
[V,D]=eig(AA,BB);
g_evals = sort(diag(D),'descend')
XD=real(g_evals(5:end));
YD=imag(g_evals(5:end));
abs_g_evals = abs(g_evals)
evals_of_A = eig(A);
display(evals_of_A);
xea=real(evals_of_A);
yea=imag(evals_of_A);
%%  %2
angle=0:0.01:2*pi;
R=1;
x=R*cos(angle);
y=R*sin(angle);
figure(1),
plot(x,y,'r');
axis equal 
hold on;
grid on;
plot(XD,YD,'o','linewidth',2);
plot(xea,yea,'x','linewidth',2);

%% %3
evals_c=(1/ts)*log(evals_of_A);
display(evals_c);

%% %4
% AA=[A B(:,1);-C(1,:) -0]; %% 11
% AA=[A B(:,2);-C(1,:) -0]; %% 12
% AA=[A B(:,1);-C(2,:) -0]; %% 21
 AA=[A B(:,2);-C(2,:) -0]; %% 22
BB=[I zeros(size(B(:,1)));zeros(size(C(1,:))) 0];
[V,D]=eig(AA,BB);
evals_channel = sort(diag(D),'descend')
XD=real(evals_channel(4:end));
YD=imag(evals_channel(4:end));
abs_evals_channel = abs(evals_channel);
display(abs_evals_channel);
%% Hankel
row = [];
Hankel = [];
n = 20;
for k=0:n-1
    for i=1:n
        row = [row,y22(mi+i+k)];
        if i==n
            Hankel = [Hankel;row];
            row = [];
        end
    end
end
[u,s,v] = svd(Hankel);
s_diag = diag(s);
display(s_diag);
%%  %circle
angle=0:0.01:2*pi;
R=1;
x=R*cos(angle);
y=R*sin(angle);
figure(2),
plot(x,y,'r');
axis square 
hold on;
grid on;
plot(XD,YD,'o');
plot(xea,yea,'x');


