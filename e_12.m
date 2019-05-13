%% Exercicio 12
%Reset do ambiente de trabalho
close all;

l = l_c;

for n = 1:length(BPM)
    J = ((M*L^2)/3) + (m * l(n)^2);
    a1 = (g * (M *(L/2) + m * l(n)) - k) / J;
    a2 = -(beta / J);

    A = [0 1 ; a1 a2];
    B = [0;  1/J];
    C = [1 0];
    D = [0];

    figure;
    bode(ss(A, B, C, D));
    grid on;
    title(sprintf('Diagrama de bode para %d BPM (l = %6f m)', BPM(n), l(n)));
end

%%
% *Comentarios:*
%
%  Uma vez que nos encontramos na presenca de um sistema de segunda ordem com um valor de qsi < 0.707, verifica-se que o diagrama de bode de amplitude apresenta um pico de ressonancia para os valores de l dado.
%
% Observa-se que o comportamento do diagrama de bode de fase é identico para os dois casos de estudo, algo que nao é surprendente visto que a alteracao do parametro l apenas afecta a localizacao dos polos e por sua vez apenas onde o decrescimo de fase acontece.
%
% Conclui-se que para diferentes valores de l, o sistema responde de forma semelhante mas com velocidades de resposta diferentes. De notar que BPM = 150 é mais rapido.