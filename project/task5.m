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
%% RMS of output y
y1 = (u_rand.Y(3).Data)/2;
y2 = (u_rand.Y(4).Data)/2;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;

y = [y1;y2];
Ryy = zeros(2);
for q=1:N
    Ryy_temp = y(:,q)*y(:,q)';
    Ryy = Ryy + Ryy_temp;
end
Ryy = Ryy/N;

display(num2str(sqrt(trace(Ryy))),'y_rms is approximately: ');

%% equation (9) and (10)
C_lyap = C'*C;
B_lyap = B*B';
Go = dlyap(A',A,C_lyap);
Gc = dlyap(A,A',B_lyap);

P_1 = sqrt(trace(B'*Go*B));
P_2 = sqrt(trace(C*Gc*C'));

display(num2str(P_1),'using (9): ');
display(num2str(P_2),'using (10): ');

%% equation (8)
trace_sum = 0;
for k = 1:401
    h_k = [y11(k),y12(k);y21(k),y22(k)];
    trace_temp = trace(h_k'*h_k);
    trace_sum = trace_sum + trace_temp;
end
P_3 = sqrt(trace_sum);
display(num2str(P_3),'using (8): ');
