function main 
clc;
clear;  
T=1;%扫描周期
N=120/T;%采样次数  
CR=15;%基站位置
BS=zeros(2,7);
BS(:,1)=[0,0];
BS(:,2)=[0,3*CR];
BS(:,3)=[3*CR,3*CR];
BS(:,4)=[3*CR,0];
X=zeros(4,N);%目标真实位置、速度 
X(:,1)=[30,2,40,2];%目标初始位置(30,40)速度(2,2)
Z=zeros(2,N);%基站对位置的观测 
a=exprnd(1,1,4);%指数分布的过程噪声
Q=diag(a);%过程噪声
b=exprnd(3,1,2);
R=60*diag(b); %观测噪声指数分布
A=[1,T,0,0;
   0,1,0,0;
   0,0,1,T;
   0,0,0,1]; %状态转移矩阵 
G=[1,0,0,0;
   0,0,1,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for t=2:N      
    X(:,t)=A*X(:,t-1)+sqrtm(Q)*randn(4,1);%目标真实轨迹     
end
MX=X(1,:);
MY=X(3,:);
for t=1:N
    yuce=SI(MX(1,t),MY(1,t));
    Z(:,t)=[yuce(1);yuce(2)]+sqrtm(R)*randn(2,1);%对目标的观测
end

%Taylor修正
for t=1:N
    track=Taylor(Z(1,t),Z(2,t));
    TL(:,t)=[track(1);track(2)]+sqrtm(R)*randn(2,1);%对目标的观测
end

% Kalman滤波 
Xkf=zeros(4,N); %滤波后的最优估计 
Xkf(:,1)=X(:,1);%初始化 
P0=eye(4);% 误差协方差阵初始化 
for i=2:N      
    %预测方程
    Xn=A*Xkf(:,i-1);%状态预测     
    P1=A*P0*A'+Q;%预测误差协方差  
    %信息方程
    K=P1*G'/(G*P1*G'+R);%Kalman增益
    %估计方程
    Xkf(:,i)=Xn+K*(TL(:,i)-G*Xn);%状态更新     
    P0=(eye(4)-K*G)*P1;%滤波误差协方差更新 
end
% 误差分析 
for i=1:N      
    Observation(i)=RMS(X(:,i),Z(:,i));%滤波前的误差  
    TaylorFilter(i)=RMS(X(:,i),TL(:,i));%Taylor修正后的误差
    KalmanFilter(i)=RMS(X(:,i),Xkf(:,i));%卡尔曼滤波后的误差 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 画图 
figure 
hold on;box on;  
plot(BS(1,:),BS(2,:),'ko');%基站
plot(X(1,:),X(3,:),'-k');%真实轨迹 
plot(Z(1,:),Z(2,:),'-mx');%观测轨迹  
plot(TL(1,:),TL(2,:),'-b');%Taylor修正轨迹
plot(Xkf(1,:),Xkf(3,:),'-r');%Kalman滤波轨迹 
legend('基站位置','真实轨迹','观测轨迹','Taylor修正轨迹','卡尔曼滤波后轨迹') ;
xlabel('横坐标 X/m'); ylabel('纵坐标 Y/m'); 
title('轨迹图');
figure 
hold on;box on;  
plot(Observation,'-ko','MarkerFace','m');
plot(TaylorFilter,'-kv','MarkerFace','b');
plot(KalmanFilter,'-ks','MarkerFace','r');
legend('滤波前误差','Taylor修正后误差','卡尔曼滤波后误差');
xlabel('观测时间/s'); ylabel('误差值'); 
title("均方误差对比");