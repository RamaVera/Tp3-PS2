clear all
close all
clc

load('tp3_kalman.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toEvaluate = 'position'; %ALWAYS

%Choose case
%Case = 'a';
Case = 'b';

%Auxiliares
I = [1,0;0,1];
O = zeros(size(I));
O_6 =[O,O,O;O,O,O;O,O,O];
final = length(p) ;
Xsave = [];
Xpred = [];
E = [];

%Parametros del modelo continuo
sigma_a_dot = 10^-1 ;
q = sigma_a_dot^2*I;
Q = [O,O,O;O,O,O;O,O,q]; %Cov ruido de proceso
A=[O,I,O;O,O,I;O,O,O];

sigma_p = sqrt(100);
sigma_v = sqrt(10);
sigma_a = sqrt(1);

if (strcmp(toEvaluate,'position'))
    C=[I,O,O];
    if (strcmp(Case,'a'))
        r = 10^3;
    end
    if (strcmp(Case,'b'))
        r = 10^-3;
    end
    R =r * diag([sigma_p^2 sigma_p^2]);      
end

etha = mvnrnd(zeros(length(R),1),R,final)';

%Discretizacion
%h=1;
Ad = expm(A*h);
Qd = [q*h^5/20,q*h^4/8,q*h^3/6;q*h^4/8,q*h^3/3,q*h^2/2;q*h^3/6,q*h^2/2,q*h];

%Condiciones Iniciales
x0_0 = [40 -200 0 0 0 0]';
P0_0 = diag([10^4 10^4 10^2 10^2 10 10 ]);

%Test de Observabilidad
L = rank(obsv(Ad,C));
if (L == length(Ad))
    disp('El sistema es completamente observable')
else
    disp('El sistema NO es completamente observable')
end
disp(['La cantidad de estados observables es: ',num2str(L)])

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
    if (strcmp(toEvaluate,'position'))
        Yk = [p(k,:)]'+ etha(:,k);
    end
    if (strcmp(toEvaluate,'velocity'))
        Yk = [v(k,:)]'+ etha(:,k);
    end
    if (strcmp(toEvaluate,'aceleration'))
        Yk = [a(k,:)]'+ etha(:,k);
    end
    if (strcmp(toEvaluate,'position and velocity'))
        Yk = [p(k,:)' ; v(k,:)']+ etha(:,k);
    end

    %Prediccion
    X_k_kminus = Ad * X_kminus_kminus ;
    P_k_kminus = Ad * P_kminus_kminus * Ad' + Qd ;
    
    %Actualizacion
    K_k =  P_k_kminus * C' * inv( C * P_k_kminus * C' + R);
    X_k_k =  X_k_kminus + K_k * (Yk - C * X_k_kminus );
    P_k_k = (eye(size(K_k*C)) - K_k*C) * P_k_kminus ;
    %P_k_k = (eye(size(K_k*C)) - K_k * C)* P_k_kminus * (eye(size(K_k*C)) - K_k*C)' +  K_k * R * K_k';

    Xsave = [Xsave (X_kminus_kminus) ];
    Xpred = [Xpred X_k_kminus];
    E =[E (Yk - C * X_k_kminus )];
end
%% Comparacion de Estados y Estimaciones


%Trayectoria
figure (1)
hold on
plot(p(:,1)' + etha(1,:),p(:,2)' + etha(2,:), 'linestyle', 'none', ...
    'marker', '.')
plot(Xsave(1,:),Xsave(2,:))
plot(p(:,1),p(:,2))
legend({'Medida','Trayectoria medida','Trayectoria real'})
title('Trayectoria')
saveas(gcf, 'trayectoria.png')

%Autocorrelacion de la innovacion
figure(2)
subplot(2,1,1)
[c,lags] = xcov(E(1,:));
stem(lags,c)
title('Autocorrelacion de la innovación en x')
subplot(2,1,2)
[c,lags] = xcov(E(2,:));
stem(lags,c)
title('Autocorrelacion de la innovación en y')
saveas(gcf, 'innovacion.png')



