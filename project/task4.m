clc
close all
clear all
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;

%% 1
display(['mean of input1 is ', num2str(mean(u1))])
display(['mean of input2 is ', num2str(mean(u2))])

%% 2
u=[u1;u2];
Ruu=zeros(2);
Ruus=zeros([2 2 401]);
lag=200;
j=1;
for la=-200:200
    Ruu=zeros(2);
    Ruu_temp=zeros(2);
    for q=201:N-200
        Ruu_temp=u(:,la+q)*u(:,q)';
        Ruu=Ruu+Ruu_temp;
    end
    Ruu=Ruu/N;
    Ruus(:,:,j)=Ruu;
    j=j+1;
end
la=(-200:200)*ts;
figure(1),
subplot(221);
plot(la,reshape(Ruus(1,1,:),1,length(la)));
xlabel('\tau(s)');
ylabel('Ruu\_11');

subplot(222);
plot(la,reshape(Ruus(1,2,:),1,length(la)));
xlabel('\tau(s)');
ylabel('Ruu\_12');

subplot(223);
plot(la,reshape(Ruus(2,1,:),1,length(la)));
xlabel('\tau(s)');
ylabel('Ruu\_21');

subplot(224);
plot(la,reshape(Ruus(2,2,:),1,length(la)));
xlabel('\tau(s)');
ylabel('Ruu\_22');


%% 3
u=[u1;u2];
Ruu=zeros(2);
for q=1:N
    Ruu_temp=u(:,q)*u(:,q)';
    Ruu=Ruu+Ruu_temp;
end
Ruu=Ruu/N;
display(['Ruu[0] is approximate ', mat2str(Ruu,1)]);

%% 4
u=[u1;u2];
y=[y1;y2];
Ryu=zeros(2);
lag=0.2/ts;
Ryus=zeros([2 2 (0.2+2)/ts+1]);
j=1;
for la=-lag:lag*10
    Ryu=zeros(2);
    Ryu_temp=zeros(2);
    for q=lag+1:N-2/ts
        Ryu_temp=y(:,la+q)*u(:,q)';
        Ryu=Ryu+Ryu_temp;
    end
    Ryu=Ryu/N;
    Ryus(:,:,j)=Ryu;
    j=j+1;
end
la=(-lag:lag*10)*ts;

Ryusp1=Ryus(:,1,:)/var(u1);
Ryusp2=Ryus(:,2,:)/var(u2);
y1pp=reshape(Ryusp1,2,89);
y2pp=reshape(Ryusp2,2,89);

figure(2)
subplot(221)
plot(la,y1pp(1,:),'*','linewidth',2);
xlabel('\tau(s)');
ylabel('y\_11');
grid on
axis([-0.2 2 -0.1 0.1]);

subplot(223)
plot(la,y1pp(2,:),'*','linewidth',2);
grid on
axis([-0.2 2 -0.1 0.1]);
xlabel('\tau(s)');
ylabel('y\_21');

subplot(222)
plot(la,y2pp(1,:),'*','linewidth',2);
xlabel('\tau(s)');
ylabel('y\_12');
grid on
axis([-0.2 2 -0.1 0.1]);

subplot(224)
plot(la,y2pp(2,:),'*','linewidth',2);
grid on
axis([-0.2 2 -0.1 0.1]);
xlabel('\tau(s)');
ylabel('y\_22');
