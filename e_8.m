%% Exercicio 8
%Reset do ambiente de trabalho
clear;
close all;

%Tempo de simulacao
ttotal = 10;

%Condicoes iniciais
L = 0.5;
M = 0.15;
l = 0.4;
m = 0.2;
k = 3;
beta = 1;
g = 9.8;
T = 0;

%Calculo das matrizes espaço de estado
J = ((M*L^2)/3) + (m * l^2);
a1 = (g * (M *(L/2) + m * l) - k) / J;
a2 = -(beta / J);

A = [0 1 ; a1 a2];
B = [0;  1/J];
C = [1 0; 0 1];
D = [ 0; 0];

%Calculo dos vectores proprios e valores proprios
[vec_prop, val_prop] = eig(A);

%Calculo numero de vectores proprios
[s ~] = size(vec_prop);

%Simular para cada vector proprio
for n = 1:s
    x0 = vec_prop(:,n);
    sim('space_state_model');
    teta = y(:,1);
    teta_p = y (:,2);
    plot(teta, teta_p, 'DisplayName', sprintf('x0 = [%0f %0f]', x0(1), x0(2)));
    hold on;
end

title('Espaço de estados');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'Latex');
legend('Location', 'NorthEast');

%%
% *Comentarios:*
%
% Verifica-se o esperado teóricamente, escolher como condições iniciais os vectores próprios do sistema leva a uma resposta de trajectória rectilínea no plano de fase.