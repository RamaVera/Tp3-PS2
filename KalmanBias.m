clear all
close all
clc

load('tp3_kalman.mat');


%Auxiliares
I = [1,0;0,1];
O = zeros(size(I));
O_6 =[O,O,O;O,O,O;O,O,O];
final = 351 ;
Xsave = [];
Ysave = [];

%Parametros del modelo
A=[O,I,O;O,O,I;O,O,O];
C=[I,O,O];
%C=[O,I,O];
%C=[I,I,O];

sigma_a_dot = 10^-2 ;
q = sigma_a_dot*I;
Q =[O,O,O;O,O,O;O,O,q]; %Cov ruido de proceso

sigma_p = 100;
sigma_v = 10;
sigma_a = 1;

R = diag([sigma_p^2 sigma_p^2]);
%R = diag([sigma_v sigma_v])
%R = diag([sigma_p sigma_p sigma_v sigma_v]);

%Discretizacion
%h=1;
Ad = expm(A*h);
Qd = [q*h^5/20,q*h^4/8,q*h^3;q*h^4/8,q*h^3/3,q*h^2/2;q*h^3/6,q*h^2/2,q*h];

x0_0 = [40 -200 0 0 0 0]';
P0_0=diag([10^6 10^6 100^2 100^2 10^2 10^2 ]);

%Inicializacion

X_k_k = x0_0;
P_k_k = P0_0;

for k = 2:final
    
    etha = mvnrnd(zeros(2,1),R)';
    Ykplus = p(k)+ etha;
    %Ykplus = v(k);
    %Ykplus = [p(k);v(k)];
    
    %Prediccion
    
    X_kplus_k = Ad * X_k_k ;
    P_kplus_k = Ad * P_k_k * Ad' + Qd ;
    
    %Actualizacion
    
    K_kplus =  P_kplus_k * C' * inv( C * P_kplus_k * C' + R);
    X_kplus_kplus =  X_kplus_k + K_kplus * (Ykplus - C * X_kplus_k );
    P_kplus_kplus = (eye(size(K_kplus*C)) - K_kplus*C) * P_kplus_k ;

    Xsave = [Xsave ;(X_kplus_kplus)' ];
    Ysave = [Ysave ;(C * X_kplus_kplus)'];
end


figure

subplot(2,1,1)
hold on
plot(1:final-1,Ysave(:,1))
plot(1:final,p(:,1))

subplot(2,1,2)
hold on
plot(1:final-1,Ysave(:,2))
plot(1:final,p(:,2))
legend({'Estimacion','Posicion Medida'})

figure(2)






