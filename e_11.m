%% Exercicio 11
%Reset do ambiente de trabalho
close all;

ttotal = 15;

T_v = [ 0.03 0.05 ];

threshold_v = 0.1;

tetaf = figure;
dtetaf = figure;
T_appf = figure;

for n = 1:length(T_v)
    T = T_v(n);
    l = l_c(2);

    sim('relojoaria_model');

    figure(tetaf);
    hold on;
    plot(t, teta, 'DisplayName', sprintf('Binario externo: %f', T));

    figure(dtetaf);
    hold on;
    plot(t, dteta, 'DisplayName', sprintf('Binario externo: %f', T));

    figure(T_appf);
    hold on;
    plot(t, T_app, 'DisplayName', sprintf('Binario externo: %f', T));

    [pcs, locs] = findpeaks(teta);

    waT = zeros(length(locs)-1, 1);

    for i = 1:(length(locs)-1)
        TaT = t(locs(i+1) )- t(locs(i));
        waT(i) = (2*pi)/TaT;
    end

    BPM_T = (60*mean(waT))/pi;

    fprintf('Para 150 BPM pretendidos e binario externo = %f tem-se %f BPM\n', T, BPM_T);
end

figure(tetaf);
grid on;
title('Posicao angular do metronomo');
xlabel('Tempo [s]');
ylabel('\theta [rad]');
legend('Location', 'NorthEastOutside');

figure(dtetaf);
grid on;
title('Velocidade angular do metronomo');
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','Interpreter','Latex');
legend('Location', 'NorthEastOutside');

figure(T_appf);
grid on;
title('Binario aplicado no metronomo');
xlabel('Tempo [s]');
ylabel('T/J');
legend('Location', 'NorthEastOutside');

%%
% *Comentarios:*
%
% Para um binario externo de 0.03 pode-se observar que existe um descaimento da amplitude lento na resposta do sistema, acabando por estabilizar. 
%
% Para um binario externo de 0.05 pode-se observar que existe um crescimento da amplitudo lento na resposta do sistema, acabando por estabilizar.
%
% Tambem se verifica que para ambos os casos a frequencia e afectada, aumentando. Isto vai de acordo ao esperado fisicamente por esta forca aplicada ao sistema no sentido do movimento faz com o que o sistema adquira maior velocidade angular e por sua altera a frequencia de oscilacao do mesmo.