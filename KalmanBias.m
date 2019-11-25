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
P0_0 = diag([10^4 10^4 10^2 10^2 10 10 ]);

for k = 1:final
    
    %Inicializacion
    if k == 1
        X_kminus_kminus = x0_0
        P_kminus_kminus = P0_0
    else
        X_kminus_kminus = X_k_k
        P_kminus_kminus = P_k_k
    end
    
    %Valor de la medicion
    etha = mvnrnd(zeros(2,1),R)';
    Yk = p(k)+ etha
    %Ykplus = v(k);
    %Ykplus = [p(k);v(k)];

    %Prediccion
    X_k_kminus = Ad * X_kminus_kminus 
    P_k_kminus = Ad * P_kminus_kminus * Ad' + Qd 
    
    %Actualizacion
    K_k =  P_k_kminus * C' * inv( C * P_k_kminus * C' + R);
    X_k_k =  X_k_kminus + K_k * (Yk - C * X_k_kminus )
    P_k_k = (eye(size(K_k*C)) - K_k*C) * P_k_kminus 

    Xsave = [Xsave ;(X_kminus_kminus)' ]
    Ysave = [Ysave ;(C * X_kminus_kminus)']
end


figure

subplot(2,1,1)
hold on
plot(1:final,Ysave(:,1))
plot(1:final,p(:,1))
legend({'Estimacion','Posicion Medida'})


subplot(2,1,2)
hold on
plot(1:final,Ysave(:,2))
plot(1:final,p(:,2))
legend({'Estimacion','Posicion Medida'})






