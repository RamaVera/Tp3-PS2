clear all
close all
clc

load('tp3_kalman.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toEvaluate= 'position';
%toEvaluate= 'velocity';
%toEvaluate= 'position and velocity';

%Auxiliares
I = [1,0;0,1];
O = zeros(size(I));
O_6 =[O,O,O;O,O,O;O,O,O];
final = length(p) ;
Xsave = [];
Xpred = [];
Ysave = [];



%Parametros del modelo continuo
sigma_a_dot = 10^-1 ;
q = sigma_a_dot^2*I;
Q = [O,O,O;O,O,O;O,O,q]; %Cov ruido de proceso
A=[O,I,O;O,O,I;O,O,O];
sigma_p = 10;
sigma_v = sqrt(10);
sigma_a = 1;
if (strcmp(toEvaluate,'position'))
    C=[I,O,O];
    R = diag([sigma_p^2 sigma_p^2]);
end
if (strcmp(toEvaluate,'velocity'))
  	C=[O,I,O];
    R = diag([sigma_v^2 sigma_v^2]);
end
if (strcmp(toEvaluate,'position and velocity'))
  	C=[I,O,O;O,I,O];
    R = diag([sigma_p^2 sigma_p^2 sigma_v^2 sigma_v^2]);
end


%Discretizacion
%h=1;
Ad = expm(A*h);
Qd = [q*h^5/20,q*h^4/8,q*h^3/6;q*h^4/8,q*h^3/3,q*h^2/2;q*h^3/6,q*h^2/2,q*h];

%Condiciones Iniciales
x0_0 = [40 -200 0 0 0 0]';
P0_0 = diag([10^4 10^4 10^2 10^2 10 10 ]);

%Test de Observabilidad
if (rank(obsv(Ad,C))== length(Ad))
    disp('El sistema es completamente observable')
else
    disp('El sistema NO es completamente observable')
end

%Algoritmo de Kalman

for k = 1:final
    
    %Inicializacion
    if k == 1
        X_kminus_kminus = x0_0;
        P_kminus_kminus = P0_0;
    else
        X_kminus_kminus = X_k_k;
        P_kminus_kminus = P_k_k;
    end
    
    %Valor de la medicion
    etha = mvnrnd(zeros(length(R),1),R)';
    Yk = p(k)+ 0*etha;
    %Ykplus = v(k);
    %Ykplus = [p(k);v(k)];

    %Prediccion
    X_k_kminus = Ad * X_kminus_kminus ;
    P_k_kminus = Ad * P_kminus_kminus * Ad' + Qd ;
    
    %Actualizacion
    K_k =  P_k_kminus * C' * inv( C * P_k_kminus * C' + R);
    X_k_k =  X_k_kminus + K_k * (Yk - C * X_k_kminus );
    P_k_k = (eye(size(K_k*C)) - K_k*C) * P_k_kminus ;
    %P_k_k = (eye(size(K_k*C)) - K_k * C)* P_k_kminus * (eye(size(K_k*C)) - K_k*C)' +  K_k * R * K_k';

    Xsave = [Xsave (X_kminus_kminus) ];
    Ysave = [Ysave (C * X_kminus_kminus)];
    Xpred = [Xpred X_k_kminus];
end
Xsave = Xsave(:,1:final);
Ysave = Ysave(:,1:final);


figure(1)

subplot(2,1,1)
hold on
plot(1:final,Ysave(1,:))
plot(1:final,p(:,1))
legend({'Estimacion','Posicion Medida'})

subplot(2,1,2)
hold on
plot(1:final,Ysave(2,:))
plot(1:final,p(:,2))
legend({'Estimacion','Posicion Medida'})

figure (2)
hold on
plot(Ysave(1,:),Ysave(2,:))
plot(p(:,1),p(:,2))
legend({'Estimacion','Posicion Medida'})

% figure(3)
% subplot(2,1,1)
% hold on
% plot(1:final,Xpred(1,:));
% plot(1:final,Xsave(1,:));
% 
% subplot(2,1,2)
% hold on
% plot(1:final,Xpred(2,:));
% plot(1:final,Xsave(2,:));




