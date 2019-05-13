%% Exercicio 13
%Reset do ambiente de trabalho
close all;

%Condicoes iniciais
L = 0.25;
M = 0.1;
k = 0.35;
beta = 0.001;
g = 9.8;
x0 = [pi/4 0];
T = 0;
l = l_c(1);

J = ((M*L^2)/3) + (m * l^2);
a1 = (g * (M *(L/2) + m * l) - k) / J;
a2 = -(beta / J);

A = [0 1 ; a1 a2];
B = [0;  1/J];
C = [1 0; 0 1];
D = [ 0; 0];

[~,fpeak] = getPeakGain(ss(A,B,C,D));

wp = fpeak;

syms m_s positive

J = M*L^2/3+m_s*l^2;
a = (k-g*((M*L)/2 + m_s*l))/J;

wn = sqrt(a);
qsi = beta/(2*sqrt(a*J.^2));

m_c = solve(wp == wn*sqrt(1-(2*(qsi^2))), m_s);
m_c1 = double(m_c(1));

fprintf('Massa real: %f kg\n', m);
fprintf('Massa estimada: %f kg\n', m_c1);

%%
% *Comentarios:*
%
% Pode-se verificar que o valor estimado e proximo do real, logo a estrategia de medicao e apropriada.