clear all
close all
clc

load('tp3_kalman.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%toEvaluate= 'position';
%toEvaluate= 'velocity';
toEvaluate= 'aceleration';
%toEvaluate= 'position and velocity';

%Auxiliares
I = [1,0;0,1];
O = zeros(size(I));
O_6 =[O,O,O;O,O,O;O,O,O];
final = length(p) ;
Xsave = [];
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
    R = diag([sigma_p^2 sigma_p^2]);
end
if (strcmp(toEvaluate,'velocity'))
  	C=[O,I,O];
    R = diag([sigma_v^2 sigma_v^2]);
end
if (strcmp(toEvaluate,'aceleration'))
  	C=[O,O,I];
    R = diag([sigma_a^2 sigma_a^2]);
end
% if (strcmp(toEvaluate,'position and velocity'))
%   	C=[I,O,O;O,I,O];
%     R = diag([sigma_p^2 sigma_p^2 sigma_v^2 sigma_v^2]);
%  
% end
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
    E =[E (Yk - C * X_k_kminus )];
end
%% Comparacion de Estados y Estimaciones

%Posicion en X
figure(1)
subplot(2,1,1)
hold on
plot(1:final,Xsave(1,:))
plot(1:final,p(:,1))
legend({'Estimacion','Posicion Medida'})
title('Posición en x')

%Posicion en Y
subplot(2,1,2)
hold on
plot(1:final,Xsave(2,:))
plot(1:final,p(:,2))
legend({'Estimacion','Posicion Medida'})
title('Posición en y')
saveas(gcf, 'posicion.png')

%Velocidad en X
figure(2)
subplot(2,1,1)
hold on
plot(1:final,Xsave(3,:))
plot(1:final,v(:,1))
legend({'Estimacion','Posicion Medida'})
title('Velocidad en x')

%Velocidad en Y
subplot(2,1,2)
hold on
plot(1:final,Xsave(4,:))
plot(1:final,v(:,2))
legend({'Estimacion','Posicion Medida'})
title('Velocidad en y')
saveas(gcf, 'velocidad.png')

%Aceleracion en X
figure(3)
subplot(2,1,1)
hold on
plot(1:final,Xsave(5,:))
plot(1:final,a(:,1))
legend({'Estimacion','Posicion Medida'})
title('Aceleracion en x')

%Aceleracion en Y
subplot(2,1,2)
hold on
plot(1:final,Xsave(6,:))
plot(1:final,a(:,2))
legend({'Estimacion','Posicion Medida'})
title('Aceleracion en y')
saveas(gcf, 'aceleracion.png')

%Trayectoria
figure (4)
hold on
plot(Xsave(1,:),Xsave(2,:))
plot(p(:,1),p(:,2))
plot(p(:,1)' + etha(1,:),p(:,2)' + etha(2,:), 'linestyle', 'none', ...
    'marker', '.')
legend({'Estimacion','Posicion real', 'Posicion medida'})
title('Trayectoria')
saveas(gcf, 'trayectoria.png')

%Autocorrelacion de la innovacion
figure(5)
title('Autocorrelacion de la innovación')
subplot(2,1,1)
[c,lags] = xcov(E(1,:));
stem(lags,c)
subplot(2,1,2)
[c,lags] = xcov(E(2,:));
stem(lags,c)
saveas(gcf, 'innovacion.png')

%Trayectoria con zoom
figure (6)
hold on
plot(Xsave(1,:),Xsave(2,:))
plot(p(:,1),p(:,2))
plot(p(:,1)' + etha(1,:),p(:,2)' + etha(2,:), 'linestyle', 'none', ...
    'marker', '.')
legend({'Estimacion','Posicion real', 'Posicion medida'})
title('Trayectoria - Con aumento')
xlim([-5000 -3000])
saveas(gcf, 'trayectoria_aumento.png')


e_x = sqrt(sum((transpose(Xsave(1,:)) - p(:,1)).^2))
e_y = sqrt(sum((transpose(Xsave(2,:)) - p(:,2)).^2))
last_diff = transpose(Xsave(2,end)) - p(end,2)
max_diff_x = max(transpose(Xsave(1,:)) - p(:,1))
max_diff_y = max(transpose(Xsave(2,:)) - p(:,2))

