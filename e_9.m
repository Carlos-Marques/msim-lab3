%% Exercicio 9
%Reset do ambiente de trabalho
clear;
close all;

%Tempo de simulação
ttotal = 10;

%Condições iniciais
L = 0.25;
M = 0.1;
k = 0.35;
beta = 0.001;
g = 9.8;
x0 = [pi/4 0];
T = 0;

%Funções auxiliares
aux= @(m, l) k-g*(((M*L)/2) + m*l);
J = @(m, l) ((M*L^2)/3) + (m * l^2);
qsi = @(m, l) beta./(2 * sqrt(aux(m, l).*J(m,l)));
wn = @(m, l) sqrt(aux(m, l) / J(m,l));                          %frequencia das oscilacoes naturais
wa = @(m, l) wn(m, l) * sqrt(1 - qsi(m, l).^2);             %frequencia das oscilacoes amortecidas

%Valores para serem testados
m = linspace(0, 0.5, 1000);
l = linspace(0.05, L, 1000);

BPM = [50, 150];

BPM_c = zeros(length(m), length(l));

for n = 1:length(m)
    for i = 1:length(l)
        qsi_c = qsi(m(n), l(i));

        %verifica se oscilacao e possivel
        if imag(qsi_c) == 0 && real(qsi_c) < 1
            BPM_c(n, i) = (60 * wa( m(n) , l(i) ))/pi;
        else
            BPM_c(n, i) = NaN;
        end
    end
end

l_ind = zeros(length(m), length(BPM));        %vector para guardar indices de l
error = zeros(length(m), 1);                            %vector para guardar erros calculados entre BPM pretendido e calculado

%testar os valores para m de modo a encontrar os valores de l que originam BPM mais proximo do pretendido
for n = 1:length(m)
    for i = 1:length(BPM)
        [error_c, l_ind(n, i)] = min(abs(BPM_c(n, :) - BPM(i)));
        error(n) = error(n) + error_c;
    end
end

[~, m_ind] = min(error);
m_c = m(m_ind);

fprintf("Massa calulada: %f kg\n", m_c);

l_c = l(l_ind(m_ind, :));

BPM_v = BPM_c(m_ind, l_ind(m_ind, :));      %valores de BPM calculados
BPM_s = zeros(size(BPM));                           % vector para guardar BPM simulados

%simulacao do sistema com os valores de m e l calulados e verificar a frequencia 
for n = 1:length(BPM)
    m = m_c;
    l = l_c(n);

    %Calculo das matrizes espaço de estado
    J = ((M*L^2)/3) + (m * l^2);
    a1 = (g * (M *(L/2) + m * l) - k) / J;
    a2 = -(beta / J);

    A = [0 1 ; a1 a2];
    B = [0;  1/J];
    C = [1 0; 0 1];
    D = [ 0; 0];

    sim('space_state_model');
    
    teta = y(:,1);

    [pcs, locs] = findpeaks(teta);

    waL = zeros(length(locs)-1, 1);

    %calculo da frequencia oscilatoria amortecida
    for i = 1:(length(locs)-1)
        TaL = t(locs(i+1))-t(locs(i));
        waL(i) = (2*pi)/TaL;
    end

    BPM_s(n) = (60*mean(waL))/pi;

    qsi_s = qsi(m, l);
    wn_s = wn(m, l);

    fprintf('BPM = %d:\n', BPM(n));
    fprintf('\t l = %f m\n', l_c(n));
    fprintf('\t BPM teoricamente aproximado: %f\n', BPM_v(n));
    fprintf('\t BPM obtido da simulacao: %f\n', BPM_s(n));

    figure;
    plot(t, teta, 'DisplayName', '\theta (t)');
    hold on;
    %envolventes superior e inferior
    plot(t, x0(1)*exp (-(qsi_s)*(wn_s)*t), 'r', 'DisplayName', 'Envolvente superior');
    plot(t, -x0(1)*exp (-(qsi_s)*(wn_s)*t), 'r',  'DisplayName', 'Envolvente inferior');

    title(sprintf('Posicao angular do metronomo a %d BPM', BPM(n)));
    xlabel('Tempo [s]');
    ylabel('\theta [rad]');
    legend('Location', 'NorthEast');
end

%%
% *Comentários:*
%
% Observa-se que o comportamento do sistema vai de acordo com o esperado, verificando-se assim que o dimensionamento da massa e da posiçao da mesma estão correctos. 
% Verifica-se também que o decaimento, demonstrado pela envolvente, esta directamente relacionado com o valor da frequencia de oscilaçoes do sistema, observando-se entao que para um maior BPM 
% existe decaimento mais rapido. As diferenças obtidas entre os valores esperados e obtidos podem ser aproximados aumentando o numero de passos (aumentando o numero de pontos testados para m e l).