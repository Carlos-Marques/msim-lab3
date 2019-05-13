%% Modelação e Simulação - Laboratório 3
%
% Carlos Silva - 81323 
%
% Ricardo Espadinha - 84178
%
% Grupo 33 - Segunda-feira, 10:00h LSDC1
%
% 2018/2019 - 2º Semestre

%% Simulação das Equações de Estado do Sistema (Questão 5)
%Reset do ambiente de trabalho
clear;
close all;

%Declaração de Variáveis
L = 0.5;        %Comprimento da Barra (m)
M = 0.15;       %Massa da Barra (kg)
l = 0.4;        %Distância da Massa ao eixo derotação (m)
m = 0.2;        %Massa da massa adicional (kg)
k = 3;          %Constante da mola (Nm/rad)
beta = 0.1;     %Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Aceleração da Gravidade (m/s^2)

stop = 5;       %Tempo de paragem da simulação

%Condições Iniciais
x0 = [0 pi/4];

%Binário Externo
T = 0;

%Calculo do Momento de Inércia
J = ((3*m*(l^2))+(M*(L^2)))/3;

%Simulação
sim('Q5');

%Criação e Informação de cada plot
plot(tout, x1);
title('Gráfico 1 - Posição angular do Metrónomo');
xlabel('Tempo [s]');
ylabel('\theta [rad]');
grid on;

figure;
plot(tout,x2, 'r');
title('Gráfico 2 - Velocidade Angular do Metrónomo');
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

figure;
plot(x1, x2, 'g');
title('Gráfico 3 - Espaço de Estados');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

%%
% *Comentários:* 
% 
% Nesta questão, implementou-se um diagrama em SIMULINK do modelo 
% linearizado (de forma a simular as equações de estado obtidas na Questão 
% 2 da parte Teórica).
% As variáveis de estado do sistema são x1 (posição angular) e x2
% (velocidade algular).
% As condições iniciais (x0 = [0 pi/4]) e o binário externo (T = 0) estão
% de acordo com o pedido no enunciado.

%% Simulação das Equações de Estado do Sistema através do bloco de espaço de estados (Questão 6) 
%
% >>> A = [0 1; ((g*((m*l)+((M*L)/2)))-k)/J -(beta/J)]
%
% >>> B = [0; 1/J]
%
% >>> C = [1 0; 0 1]
%
% >>> D = [0; 0]
%
% Tendo em conta que: J = ((3*m*(l^2))+(M*(L^2)))/3
%
% O modelo pedido (diagrama que utiliza o bloco de espaço de estados),
% encontra-se no ficheiro 'Q6'. De forma a realizar este diagrama, foi
% necessário definir as matrizes anteriores (A, B, C e D), calculadas
% através da análise teorica. No entanto, através da mesma, obteve-se C da
% seguinte forma: C = [1 0]
% Obteve-se C dessa forma, porque a saída do distema em estupo é a posição
% angular (x1(t) = theta(t)). No entanto, por motivos de análise da
% velocidade angular (x2(t)), procedeu-se a adicionar mais uma coluna à
% matriz C para este efeito, desta forma, x1(t) = y(:,1) e x2(t) = y(:,2).
% (Por motivos de coerencia algébrica, também foi necessario alterar a
% dimensão da matriz D para 2x1)

%% Simulação das Equações de Estado do Sistema para beta=0 (Questão 7)
%Reset do ambiente de trabalho
clear;
close all;

%Declaração de Variáveis
L = 0.5;        %Comprimento da Barra (m)
M = 0.15;       %Massa da Barra (kg)
l = 0.4;        %Distância da Massa ao eixo derotação (m)
m = 0.2;        %Massa da massa adicional (kg)
k = 3;          %Constante da mola (Nm/rad)
beta = 0;       %Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Aceleração da Gravidade (m/s^2)

stop = 3;     %Tempo de paragem da simulação

%Condições Iniciais
x0 = [0 pi/4];

T = 0;  %Binário Externo

%Calculo do Momento de Inércia
J = ((3*m*(l^2))+(M*(L^2)))/3;

%Simulação do Sistema
sim('Q6');

%Informação e criação de cada plot
plot(tout, y(:,1));
title(sprintf('Gráfico 4 - Posição angular do metrónomo para \\beta = 0 (coeficiente de atrito nulo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('\theta [rad]');
grid on;

figure;
plot(tout, y(:,2), 'r');
title(sprintf('Gráfico 5 - Velocidade angular do metrónomo para \\beta = 0 (coeficiente de atrito nulo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

figure;
plot(y(:,1), y(:,2), 'g');
title(sprintf('Gráfico 6 - Espaço de Estados para \\beta = 0 (coeficiente de atrito nulo) e x0 = [0 \\pi/4]'));
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

%Teste para diferentes condições iniciais
x0_conj = [pi/2 0; pi/4 -1; 0 pi/4; -pi/3 -1];

LC = size(x0_conj); %Vector com o numero de Linhas e Colunas da Matriz

figure('units', 'normalized', 'outerposition', [0 0 0.8 0.8]); %Aumento do tamanho da imagem

for i = 1:LC(1)
    
    %Definir as condições iniciais para a simulação
    x0(1) = x0_conj(i, 1);
    x0(2) = x0_conj(i, 2);
    
    %Simulação do Sistema
    sim('Q6');
    
    %Gerar os vários espaços de fase
    plot(y(:,1), y(:,2), 'DisplayName', sprintf('x1(0) = %f rad, x2(0) = %f rad/s', x0_conj(i, 1), x0_conj(i, 2)));
    hold on;
    
end

title(sprintf('Gráfico 7 - Sobreposição de Espaços de Estado correspondentes a diferentes Condições iniciais x(0) e respetivo Campo de Vectores para \\beta = 0'));
legend('Location','northeast');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;
    
% Sobreposição do campo de vectores do sistema
H_grid = linspace(-1.5, 1.5, 15);
V_grid = linspace(-10, 10, 15);

A = [0 1; ((g*((m*l)+((M*L)/2)))-k)/J -(beta/J)];   %Matriz A

ZH = zeros(numel(H_grid), numel(V_grid));
ZV = ZH;
for i = 1:numel(H_grid)
    for j = 1:numel(V_grid)
        Z = A * [H_grid(i);V_grid(j)];
        ZH(j, i) = Z(1);
        ZV(j, i) = Z(2);
    end
end

quiver(H_grid, V_grid, ZH, ZV);

%%
% *Comentários:* 
% 
% Tal como esperado para beta = 0 (Coeficiente de Atrito nulo), verifica-se
% nos gráficos 4, 5 e 6 que o sistema permanece indefinidamente com a mesma
% amplitude de oscilação (sem qualquer amortecimento).
% Analisando o gráfico 7, observa-se que, para beta = 0, o espaço de fase 
% do sistema consiste sempre em circunferências de raio constante.
% Consequentemente, só existe um ponto de equílibrio do sistema: 
% x1(0) = x2(0) = 0.
% Através do campo de vectores obtido, consegue-se perceber as direções da
% posição e da velocidade angular associadas a cada ponto.

%% Simulação das Equações de Estado do Sistema para beta=1 (Questão 7)
%Reset do ambiente de trabalho
clear;
close all;

%Declaração de Variáveis
L = 0.5;              %Comprimento da Barra (m)
M = 0.15;             %Massa da Barra (kg)
l = 0.4;              %Distância da Massa ao eixo derotação (m)
m = 0.2;              %Massa da massa adicional (kg)
k = 3;                %Constante da mola (Nm/rad)
beta = 1;             %Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Aceleração da Gravidade (m/s^2)

stop = 3;     %Tempo de paragem da simulação

%Condições Iniciais
x0 = [0 pi/4];
T = 0;  %Binário Externo

%Calculo do Momento de Inércia
J = ((3*m*(l^2))+(M*(L^2)))/3;

%Simulação do Sistema
sim('Q6');

%Informação de cada plot
plot(tout, y(:,1));
title(sprintf('Gráfico 8 - Posição angular do metrónomo para \\beta=1 (coeficiente de atrito máximo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('\theta [rad]');
grid on;

figure;
plot(tout, y(:,2), 'r');
title(sprintf('Gráfico 9 - Velocidade angular do metrónomo para \\beta=1 (coeficiente de atrito máximo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

figure;
plot(y(:,1), y(:,2), 'g');
title(sprintf('Gráfico 10 - Espaço de Estados para \\beta=1 (coeficiente de atrito máximo) e x0 = [0 \\pi/4]'));
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;


%Teste para diferentes condições iniciais

x0_conj = [pi/2 0; pi/4 -1; 0 pi/4; -pi/3 -1];

LC = size(x0_conj); %Vector com o numero de Linhas e Colunas da Matriz

figure('units', 'normalized', 'outerposition', [0 0 0.8 0.8]);

for i = 1:LC(1)
    
    %Definir as condições iniciais para a simulação
    x0(1) = x0_conj(i, 1);
    x0(2) = x0_conj(i, 2);
    
    %Simulação do Sistema
    sim('Q6');
    
    %Gerar os vários espaços de fase
    plot(y(:,1), y(:,2), 'DisplayName', sprintf('x1(0) = %f rad, x2(0) = %f rad/s', x0_conj(i, 1), x0_conj(i, 2)));
    hold on;
    
end

title(sprintf('Gráfico 11 - Sobreposição de Espaços de Estado correspondentes a diferentes Condições iniciais x(0) e respetivo Campo de Vectores para \\beta = 0'));
legend('Location','northeast');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;
    
% Sobreposição do campo de vectores do sistema
H_grid = linspace(-1.5, 1.5, 15);
V_grid = linspace(-2.5, 2.5, 15);

A = [0 1; ((g*((m*l)+((M*L)/2)))-k)/J -(beta/J)];   %Matriz A

ZH = zeros(numel(H_grid), numel(V_grid));
ZV = ZH;
for i = 1:numel(H_grid)
    for j = 1:numel(V_grid)
        Z = A * [H_grid(i);V_grid(j)];
        ZH(j, i) = Z(1);
        ZV(j, i) = Z(2);
    end
end

quiver(H_grid, V_grid, ZH, ZV);

%%
% *Comentários:* 
% 
% Tal como esperado para beta = 1 (Coeficiente de Atrito máximo), verifica-se
% nos gráficos 8, 9 e 10 que o sistema tende lentamete para 
% x1(t) = x2(t) = 0 no limite de t,nunca apresentando oscilações consideráveis (amortecimento máximo).
% Analisando o gráfico 11, observa-se que, para beta = 1, o sistema parte
% das suas condições iniciais, descrevendo uma trajectória com a mesma
% direcção (mais à frente, observa-se que esta é proporcional à direcção de
% um dos vecores próprios da matriz A), terminando sempre noutra
% trajectória (proporcional ao outro vector próprio da matriz A).
% Através do campo de vectores obtido, consegue-se perceber as direções da
% posição e da velocidade angular associadas a cada ponto.

%% Valores e Vectores próprios da matriz A para beta = 0, 0.1 e 1 (Questão 7)
%Reset do ambiente de trabalho
clear;
close all;

%Declaração de Variáveis
L = 0.5;              %Comprimento da Barra (m)
M = 0.15;             %Massa da Barra (kg)
l = 0.4;              %Distância da Massa ao eixo derotação (m)
m = 0.2;              %Massa da massa adicional (kg)
k = 3;                %Constante da mola (Nm/rad)

betas = [0 0.1 1];    %Conjunto de Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Aceleração da Gravidade (m/s^2)

%Calculo do Momento de Inércia
J = ((3*m*(l^2))+(M*(L^2)))/3;

for i = 1:length(betas)
    
    % Determinação dos valores próprios e vectores próprios da matriz A
    A = [0 1; ((g*((m*l)+((M*L)/2)))-k)/J -(betas(i)/J)];
    
    % A função 'eig(A)' retorna os valores próprios associados a A numa
    % matriz diagonal (eigenvalues) e os vectores próprios de A em coluna
    % (eigenvectores)
    [eigenvectors, eigenvalues] = eig(A);

    % Impressão dos valores próprios na Command Window
    fprintf('beta = %2f', betas(i));
    fprintf('\n\n => Valores Próprios:\n\n');
    fprintf('\t lambda1 =');
    disp(eigenvalues(1, 1))
    fprintf('\t lambda2 =');
    disp(eigenvalues(2, 2));

    % Impressão dos vectores próprios na Command Window
    fprintf('=> Vectores próprios:\n\n \t v1 = [%.4f %.4f]\n\n \t v2 = [%.4f %.4f]\n\n', eigenvectors(:, 1), eigenvectors(:, 2));
    fprintf('---------------------------------\n');
end


%%
% *Comentários:* 
% 
% Para beta = 0:
% Analizando os Valores Próprios, verifica-se que estes são complexos
% conjugados, indicando que o plano de estado do sistema é sempre 
% representado por uma circunferência (como se comprova no gráfico 7), e 
% com parte real nula, oscilando permanentemente com a mesma amplitude 
% (marginalmente estável).
%
% Para beta = 0.1:
% Analizando os Valores Próprios, verifica-se que estes são complexos
% conjugados mas com a parte real negativa, indicando que se trata de um
% termo oscilatório amortecido. Isto comprova-se através do gráfico 3 que
% apresenta um comportamento com trajetórias "em espiral" que se aproximam da origem. 
%
% Para beta = 1:
% Analizando os Valores Próprios, verifica-se que estes são ambos reais
% negativos, indicando que o sistema é estável e que o seu plano de estado
% é sempre uma rect em torno na origem, como se pode comprovar através dos
% gráficos 10 e 11.
% Para além disso, como já foi referido anteriormente, analizando
% detalhadamente o gráfico 11, verifica-se que as trajetórias partem das
% posições iniciais com a mesma direção do Vector Próprio v2 e terminam com
% a mesma direção que o Vector Próprio V1.