%% Modela��o e Simula��o - Laborat�rio 3
%
% Carlos Silva - 81323 
%
% Ricardo Espadinha - 84178
%
% Grupo 33 - Segunda-feira, 10:00h LSDC1
%
% 2018/2019 - 2� Semestre

%% Simula��o das Equa��es de Estado do Sistema (Quest�o 5)
%Reset do ambiente de trabalho
clear;
close all;

%Declara��o de Vari�veis
L = 0.5;        %Comprimento da Barra (m)
M = 0.15;       %Massa da Barra (kg)
l = 0.4;        %Dist�ncia da Massa ao eixo derota��o (m)
m = 0.2;        %Massa da massa adicional (kg)
k = 3;          %Constante da mola (Nm/rad)
beta = 0.1;     %Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Acelera��o da Gravidade (m/s^2)

stop = 5;       %Tempo de paragem da simula��o

%Condi��es Iniciais
x0 = [0 pi/4];

%Bin�rio Externo
T = 0;

%Calculo do Momento de In�rcia
J = ((3*m*(l^2))+(M*(L^2)))/3;

%Simula��o
sim('Q5');

%Cria��o e Informa��o de cada plot
plot(tout, x1);
title('Gr�fico 1 - Posi��o angular do Metr�nomo');
xlabel('Tempo [s]');
ylabel('\theta [rad]');
grid on;

figure;
plot(tout,x2, 'r');
title('Gr�fico 2 - Velocidade Angular do Metr�nomo');
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

figure;
plot(x1, x2, 'g');
title('Gr�fico 3 - Espa�o de Estados');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

%%
% *Coment�rios:* 
% 
% Nesta quest�o, implementou-se um diagrama em SIMULINK do modelo 
% linearizado (de forma a simular as equa��es de estado obtidas na Quest�o 
% 2 da parte Te�rica).
% As vari�veis de estado do sistema s�o x1 (posi��o angular) e x2
% (velocidade algular).
% As condi��es iniciais (x0 = [0 pi/4]) e o bin�rio externo (T = 0) est�o
% de acordo com o pedido no enunciado.

%% Simula��o das Equa��es de Estado do Sistema atrav�s do bloco de espa�o de estados (Quest�o 6) 
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
% O modelo pedido (diagrama que utiliza o bloco de espa�o de estados),
% encontra-se no ficheiro 'Q6'. De forma a realizar este diagrama, foi
% necess�rio definir as matrizes anteriores (A, B, C e D), calculadas
% atrav�s da an�lise teorica. No entanto, atrav�s da mesma, obteve-se C da
% seguinte forma: C = [1 0]
% Obteve-se C dessa forma, porque a sa�da do distema em estupo � a posi��o
% angular (x1(t) = theta(t)). No entanto, por motivos de an�lise da
% velocidade angular (x2(t)), procedeu-se a adicionar mais uma coluna �
% matriz C para este efeito, desta forma, x1(t) = y(:,1) e x2(t) = y(:,2).
% (Por motivos de coerencia alg�brica, tamb�m foi necessario alterar a
% dimens�o da matriz D para 2x1)

%% Simula��o das Equa��es de Estado do Sistema para beta=0 (Quest�o 7)
%Reset do ambiente de trabalho
clear;
close all;

%Declara��o de Vari�veis
L = 0.5;        %Comprimento da Barra (m)
M = 0.15;       %Massa da Barra (kg)
l = 0.4;        %Dist�ncia da Massa ao eixo derota��o (m)
m = 0.2;        %Massa da massa adicional (kg)
k = 3;          %Constante da mola (Nm/rad)
beta = 0;       %Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Acelera��o da Gravidade (m/s^2)

stop = 3;     %Tempo de paragem da simula��o

%Condi��es Iniciais
x0 = [0 pi/4];

T = 0;  %Bin�rio Externo

%Calculo do Momento de In�rcia
J = ((3*m*(l^2))+(M*(L^2)))/3;

%Simula��o do Sistema
sim('Q6');

%Informa��o e cria��o de cada plot
plot(tout, y(:,1));
title(sprintf('Gr�fico 4 - Posi��o angular do metr�nomo para \\beta = 0 (coeficiente de atrito nulo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('\theta [rad]');
grid on;

figure;
plot(tout, y(:,2), 'r');
title(sprintf('Gr�fico 5 - Velocidade angular do metr�nomo para \\beta = 0 (coeficiente de atrito nulo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

figure;
plot(y(:,1), y(:,2), 'g');
title(sprintf('Gr�fico 6 - Espa�o de Estados para \\beta = 0 (coeficiente de atrito nulo) e x0 = [0 \\pi/4]'));
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

%Teste para diferentes condi��es iniciais
x0_conj = [pi/2 0; pi/4 -1; 0 pi/4; -pi/3 -1];

LC = size(x0_conj); %Vector com o numero de Linhas e Colunas da Matriz

figure('units', 'normalized', 'outerposition', [0 0 0.8 0.8]); %Aumento do tamanho da imagem

for i = 1:LC(1)
    
    %Definir as condi��es iniciais para a simula��o
    x0(1) = x0_conj(i, 1);
    x0(2) = x0_conj(i, 2);
    
    %Simula��o do Sistema
    sim('Q6');
    
    %Gerar os v�rios espa�os de fase
    plot(y(:,1), y(:,2), 'DisplayName', sprintf('x1(0) = %f rad, x2(0) = %f rad/s', x0_conj(i, 1), x0_conj(i, 2)));
    hold on;
    
end

title(sprintf('Gr�fico 7 - Sobreposi��o de Espa�os de Estado correspondentes a diferentes Condi��es iniciais x(0) e respetivo Campo de Vectores para \\beta = 0'));
legend('Location','northeast');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;
    
% Sobreposi��o do campo de vectores do sistema
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
% *Coment�rios:* 
% 
% Tal como esperado para beta = 0 (Coeficiente de Atrito nulo), verifica-se
% nos gr�ficos 4, 5 e 6 que o sistema permanece indefinidamente com a mesma
% amplitude de oscila��o (sem qualquer amortecimento).
% Analisando o gr�fico 7, observa-se que, para beta = 0, o espa�o de fase 
% do sistema consiste sempre em circunfer�ncias de raio constante.
% Consequentemente, s� existe um ponto de equ�librio do sistema: 
% x1(0) = x2(0) = 0.
% Atrav�s do campo de vectores obtido, consegue-se perceber as dire��es da
% posi��o e da velocidade angular associadas a cada ponto.

%% Simula��o das Equa��es de Estado do Sistema para beta=1 (Quest�o 7)
%Reset do ambiente de trabalho
clear;
close all;

%Declara��o de Vari�veis
L = 0.5;              %Comprimento da Barra (m)
M = 0.15;             %Massa da Barra (kg)
l = 0.4;              %Dist�ncia da Massa ao eixo derota��o (m)
m = 0.2;              %Massa da massa adicional (kg)
k = 3;                %Constante da mola (Nm/rad)
beta = 1;             %Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Acelera��o da Gravidade (m/s^2)

stop = 3;     %Tempo de paragem da simula��o

%Condi��es Iniciais
x0 = [0 pi/4];
T = 0;  %Bin�rio Externo

%Calculo do Momento de In�rcia
J = ((3*m*(l^2))+(M*(L^2)))/3;

%Simula��o do Sistema
sim('Q6');

%Informa��o de cada plot
plot(tout, y(:,1));
title(sprintf('Gr�fico 8 - Posi��o angular do metr�nomo para \\beta=1 (coeficiente de atrito m�ximo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('\theta [rad]');
grid on;

figure;
plot(tout, y(:,2), 'r');
title(sprintf('Gr�fico 9 - Velocidade angular do metr�nomo para \\beta=1 (coeficiente de atrito m�ximo) e x0 = [0 \\pi/4]'));
xlabel('Tempo [s]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;

figure;
plot(y(:,1), y(:,2), 'g');
title(sprintf('Gr�fico 10 - Espa�o de Estados para \\beta=1 (coeficiente de atrito m�ximo) e x0 = [0 \\pi/4]'));
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;


%Teste para diferentes condi��es iniciais

x0_conj = [pi/2 0; pi/4 -1; 0 pi/4; -pi/3 -1];

LC = size(x0_conj); %Vector com o numero de Linhas e Colunas da Matriz

figure('units', 'normalized', 'outerposition', [0 0 0.8 0.8]);

for i = 1:LC(1)
    
    %Definir as condi��es iniciais para a simula��o
    x0(1) = x0_conj(i, 1);
    x0(2) = x0_conj(i, 2);
    
    %Simula��o do Sistema
    sim('Q6');
    
    %Gerar os v�rios espa�os de fase
    plot(y(:,1), y(:,2), 'DisplayName', sprintf('x1(0) = %f rad, x2(0) = %f rad/s', x0_conj(i, 1), x0_conj(i, 2)));
    hold on;
    
end

title(sprintf('Gr�fico 11 - Sobreposi��o de Espa�os de Estado correspondentes a diferentes Condi��es iniciais x(0) e respetivo Campo de Vectores para \\beta = 0'));
legend('Location','northeast');
xlabel('\theta [rad]');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on;
    
% Sobreposi��o do campo de vectores do sistema
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
% *Coment�rios:* 
% 
% Tal como esperado para beta = 1 (Coeficiente de Atrito m�ximo), verifica-se
% nos gr�ficos 8, 9 e 10 que o sistema tende lentamete para 
% x1(t) = x2(t) = 0 no limite de t,nunca apresentando oscila��es consider�veis (amortecimento m�ximo).
% Analisando o gr�fico 11, observa-se que, para beta = 1, o sistema parte
% das suas condi��es iniciais, descrevendo uma traject�ria com a mesma
% direc��o (mais � frente, observa-se que esta � proporcional � direc��o de
% um dos vecores pr�prios da matriz A), terminando sempre noutra
% traject�ria (proporcional ao outro vector pr�prio da matriz A).
% Atrav�s do campo de vectores obtido, consegue-se perceber as dire��es da
% posi��o e da velocidade angular associadas a cada ponto.

%% Valores e Vectores pr�prios da matriz A para beta = 0, 0.1 e 1 (Quest�o 7)
%Reset do ambiente de trabalho
clear;
close all;

%Declara��o de Vari�veis
L = 0.5;              %Comprimento da Barra (m)
M = 0.15;             %Massa da Barra (kg)
l = 0.4;              %Dist�ncia da Massa ao eixo derota��o (m)
m = 0.2;              %Massa da massa adicional (kg)
k = 3;                %Constante da mola (Nm/rad)

betas = [0 0.1 1];    %Conjunto de Coeficiente de Atrito (Nms/rad)

g = 9.8;        %Acelera��o da Gravidade (m/s^2)

%Calculo do Momento de In�rcia
J = ((3*m*(l^2))+(M*(L^2)))/3;

for i = 1:length(betas)
    
    % Determina��o dos valores pr�prios e vectores pr�prios da matriz A
    A = [0 1; ((g*((m*l)+((M*L)/2)))-k)/J -(betas(i)/J)];
    
    % A fun��o 'eig(A)' retorna os valores pr�prios associados a A numa
    % matriz diagonal (eigenvalues) e os vectores pr�prios de A em coluna
    % (eigenvectores)
    [eigenvectors, eigenvalues] = eig(A);

    % Impress�o dos valores pr�prios na Command Window
    fprintf('beta = %2f', betas(i));
    fprintf('\n\n => Valores Pr�prios:\n\n');
    fprintf('\t lambda1 =');
    disp(eigenvalues(1, 1))
    fprintf('\t lambda2 =');
    disp(eigenvalues(2, 2));

    % Impress�o dos vectores pr�prios na Command Window
    fprintf('=> Vectores pr�prios:\n\n \t v1 = [%.4f %.4f]\n\n \t v2 = [%.4f %.4f]\n\n', eigenvectors(:, 1), eigenvectors(:, 2));
    fprintf('---------------------------------\n');
end


%%
% *Coment�rios:* 
% 
% Para beta = 0:
% Analizando os Valores Pr�prios, verifica-se que estes s�o complexos
% conjugados, indicando que o plano de estado do sistema � sempre 
% representado por uma circunfer�ncia (como se comprova no gr�fico 7), e 
% com parte real nula, oscilando permanentemente com a mesma amplitude 
% (marginalmente est�vel).
%
% Para beta = 0.1:
% Analizando os Valores Pr�prios, verifica-se que estes s�o complexos
% conjugados mas com a parte real negativa, indicando que se trata de um
% termo oscilat�rio amortecido. Isto comprova-se atrav�s do gr�fico 3 que
% apresenta um comportamento com trajet�rias "em espiral" que se aproximam da origem. 
%
% Para beta = 1:
% Analizando os Valores Pr�prios, verifica-se que estes s�o ambos reais
% negativos, indicando que o sistema � est�vel e que o seu plano de estado
% � sempre uma rect em torno na origem, como se pode comprovar atrav�s dos
% gr�ficos 10 e 11.
% Para al�m disso, como j� foi referido anteriormente, analizando
% detalhadamente o gr�fico 11, verifica-se que as trajet�rias partem das
% posi��es iniciais com a mesma dire��o do Vector Pr�prio v2 e terminam com
% a mesma dire��o que o Vector Pr�prio V1.