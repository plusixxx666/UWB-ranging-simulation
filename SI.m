function X = SI(mx,my)
%SIALGORITHM 本函数用于实现无线定位中的SI算法
%7小区网络拓扑
BSN=4;
CR=15;%基站位置
BS=zeros(2,7);
BS(:,1)=[0,0];
BS(:,2)=[0,3*CR];
BS(:,3)=[3*CR,3*CR];
BS(:,4)=[3*CR,0];
MS = [mx;my];
%均值为0的高斯噪声
c=3*10^8;
Nmeasure=0.05*10^(-8);
Noise=c*Nmeasure*randn(1);
%均方时延扩展
T=0.1;%均方根时延程度
%tmax=T*d^0.5*10^0.4*randn;
% Ri1
R1 = sqrt(MS(1)^2 + MS(2)^2);
for i = 1: BSN-1,
    R(i) = sqrt((BS(1, i+1) - MS(1))^2 + (BS(2, i+1) - MS(2))^2);
end
for i = 1: BSN-1,
    Ri1(i) = R(i) - R1 ;+ T*R(i)^0.5*10^0.4*randn + Noise;
end
    
% W
W = eye(BSN-1);

% delt
for i = 1: BSN-1,
    K(i) = BS(1, i+1)^2 + BS(2, i+1)^2;
end
for i = 1: BSN-1,
    delt(i) = K(i) - Ri1(i)^2;
end

% Pd orthognol
I = eye(BSN-1);
coef = Ri1*Ri1';
Pd_o = I - (Ri1'*Ri1/coef);
    
% S
for i = 1: BSN-1,
    S(i, 1) = BS(1, i+1);
    S(i, 2) = BS(2, i+1);
end

% 输出：
    Za = 0.5*inv(S'*Pd_o*W*Pd_o*S)*S'*Pd_o*W*Pd_o*delt';
    X = Za;
