%% Exercicio 10
%Reset do ambiente de trabalho
close all;

BPM_s_nl = zeros(size(BPM));
BPM_s_nlr = zeros(size(BPM));

%simular o sistema nao linear para cada um dos BPM
for n = 1:length(BPM)
    l = l_c(n);

    sim('NL_model');

    [pcs, locs] = findpeaks(teta);

    TaNL = zeros(length(locs)-1, 1);
    waNL = zeros(length(locs)-1, 1);

    %calcular frequencia oscilatoria amortecida
    for i = 1:(length(locs)-1)
        TaNL(i) = t(locs(i+1)) - t(locs(i));
        waNL(i) = (2*pi)/TaNL(i);
    end

    BPM_s_nl(n) = (60*mean(waNL))/pi;
end

l_c_nl = l_c;

BPM_s_nlr = BPM_s_nl;

error = BPM - BPM_s_nlr;

error_tolerance = 0.01;     %toleracia de erro entro o valor pretendido e calculado

while( sum(abs(error) > error_tolerance ))
    for n = 1:length(BPM)
    l_c_nl(n) = l_c_nl(n) + (BPM_s_nlr(n)-BPM(n))*0.01*error_tolerance;

    l = l_c_nl(n);

    sim('NL_model');

    [pcs, locs] = findpeaks(teta);

    waNL = zeros(length(locs)-1, 1);

    %calculo da frequencia oscilatoria amortecida
    for i = 1:(length(locs)-1)
        TaNL= t(locs(i+1)) - t(locs(i));
        waNL(i) = (2*pi)/TaNL;
    end

    BPM_s_nlr(n) = (60*mean(waNL))/pi;
    end

    error = BPM - BPM_s_nlr;
end

for n = 1:length(BPM)
    fprintf('BPM = %d\n\n', BPM(n));
    fprintf('\t Sistema linearizado: %f BPM, l=%f m\n', BPM_s(n), l_c(n));
    fprintf('\t Sistema nao linear: %f BPM, l=%f m\n', BPM_s_nl(n), l_c(n));
    fprintf('\t Sistema nao linear refinado: %f BPM, l=%f m\n', BPM_s_nlr(n), l_c_nl(n));
end

%%
% *Comentários:*
%
% Utilizando o valor de l calculado anteriormente verifica-se um erro elevado das frequencias de oscilacao obtidas para o sistema nao linearizado. Apos o refinamento obtem-se entao frequencias de oscilacao proximas das pretendidas.