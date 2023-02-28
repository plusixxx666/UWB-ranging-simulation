function X = Taylor(mx,my)
%实现无线定位中的Taylor算法
%4基站网络拓扑
BSN=4;
CR=15;%基站位置
BS=zeros(2,7);
BS(:,1)=[0,0];
BS(:,2)=[0,3*CR];
BS(:,3)=[3*CR,3*CR];
BS(:,4)=[3*CR,0];

% TDOA协方差矩阵Q：
Q = eye(BSN-1);
MS=[mx my];

% 初始估计位置：
iEP = MS;

% h0:
for i = 1: BSN,
    MeaDist(i) = sqrt((MS(1) - BS(1,i))^2 + (MS(2) - BS(2,i))^2);
end
for i = 1: BSN-1,
    h0(i) = MeaDist(i+1) - MeaDist(1);
end

% 算法开始：
for n = 1:10,
    % Rn:
    R1 = sqrt(iEP(1)^2 + iEP(2)^2);
    for i =1: BSN-1,
        R(i) = sqrt((iEP(1) - BS(1,i+1))^2 + (iEP(2) - BS(2,i+1))^2);        
    end
    
    % ht:
    for i = 1: BSN-1,
        h(i) = h0(i) - (R(i) - R1);
    end
    ht = h';
    
    % Gt:
    for i = 1: BSN-1,
        Gt(i, 1) = -iEP(1)/R1 - (BS(1, i+1) - iEP(1))/R(i);
        Gt(i, 2) = -iEP(2)/R1 - (BS(2, i+1) - iEP(2))/R(i);
    end
    
    % delt:
    delt = inv(Gt'*inv(Q)*Gt)*Gt'*inv(Q)*ht;
    
    EP = iEP + delt';

    iEP = EP;
end
    X = EP;
