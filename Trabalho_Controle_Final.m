%% Trabalho Computacional - Controle I - Compensador Avaço - Atraso
%{
    Integrantes: Ana Luisa Basilio de Castro, 
                 Gabriel Ferreira Freitas de Paula, 
                 Leonardo Nunes Rubbioli Cordeiro, 
                 Rafael Salzer Simas
%}

%% Inicialização do código
% Limpando memória, fechando figuras
clc
clear all
close all

% Modificando formato de exibição dos números
format long

% Inicializando variavies:
s = tf('s');                    %variavel complexa s, F.T.
wn = 1.9;                       
zeta = 0.5;
k=0;                            % Inicializa o contador de soluções

% Declaração de vetores

% Declaração das Funções de Transferência (TF)
G = (2*s^2 + 4*s + 34)/(s*(s + 1.5)*(s + 5));
H = 1;

%% SISTEMA DE CONTROLE
%Tópicos posteriores constam os varios itens do trabalho e sua execução

%% Letra a) Lugar das Raízes
% Após declarar anteriormente a função de transferência (FT) G(s) e H(s) vamos usar o comando rlocus(FT) para exibir o gráfico do LR

% Determinando polos e zeros de G(s)
polos_G = pole(G);
zeros_G = zero(G);

% Determinando polos e zeros da FTMF; Obs: a FTMF será denotada por A 
A = feedback(G,H);
%Cálculo dos plos e zeros de A; 
polos_A = pole(A);
zeros_A = zero(A);

% Exibindo na Command Window os valores de polos e zeros de G(s) e A(s)
fprintf('_________________________________________\n');
fprintf('\n>> Função de Transferência G(s):\n');
G
fprintf('   Polos de G:\n');
disp(polos_G);
fprintf('   Zeros de G:\n');
disp(zeros_G);

fprintf('\n>> Função de Transferência em Malha Fechada A(s):\n');
A
fprintf('   Polos de A:\n');
disp(polos_A);
fprintf('   Zeros de A:\n');
disp(zeros_A);
fprintf('_________________________________________\n'); %Apenas para melhorar a visualização na Command Window

%Cálculo do wn e do zeta da FTMF (A) não compensada
wn_A = abs(polos_A(2));
zeta_A = abs(real(polos_A(2))/wn_A);             %pols_A(2) para usar o polo complexo conjugado

% Criando o gráfico do Lugar das Raízes de G(s)
%armazenando em um vetor o LR de G para cálculos posteriores
rl_G = rlocus(G);           

% Inicializando figura e configurando o gráfico
rlocus(G,'k');
grid on
title('Lugar das Raízes de G(s), sistema não compensado')
xlabel('Real')
ylabel('Imaginário')
hold on                     %hold para adicionar ao mesmo gráfico os polos e zeros

% Plotando os polos e zeros de A(s)
plot(real(polos_A), imag(polos_A), 'rx', 'MarkerSize', 7.5, 'LineWidth', 1); % Polos (x vermelho)
plot(real(zeros_A), imag(zeros_A), 'ro', 'MarkerSize', 7.5, 'LineWidth', 1); % Zeros (o vermelho)

% Adicionando ζ = 0.5 e wn1 = 1.9 rad/s ao mesmo gráfico
%{
    Breve explicação dos cáculos pt1 - Construção das retas
    Para construir as retas desejadas foi usado coordenadas polares, aproveitando a relação entre os angulos que obtemos 
    facilmente através dos valores de ζ e wn. A reta é contruída com referência de ângulo na origem e no sentido anti 
    horário, o ângulo será dado por theta_reta. O raio é arbitrário de acordo com a visualização desejada do
    gráfico(escala), aqui variamos de 0 a 10.
%}
beta = acos(zeta);
theta_reta = pi - beta;
theta_reta_graus = 180*theta_reta/pi;

% Criando a reta a partir da origem
%Uso de coordenadas polares e transformação matemática para obter x e y para o plot
r = linspace(0, 10, 100);      %Raio de 0 a 10 (ajustável)
x = r * cos(theta_reta);       %Conversão coord polares
y = r * sin(theta_reta);

% Plotando reta para ζ = 0.5
plot(x, y, 'c--', 'LineWidth', 0.5); %reta tracejada ciano

% Criando círculo de raio 1.9 para wn = 1.9rad/s
%{
    Breve explicação pt2 - Contrução dos círculos
    Para os circulos foi definido raio = wn. No plot será aplicado vários valores de x_circ e y_circ dentro da 
    circunferência de raio wn. Os vetores de x,y_circ são resultado da multiplicação do raio (wn) por um conjunto 
    de valores de sen e cos aplicados a variação do angulo da circunferência, que vai de 0 a 2pi, ou seja, uma 
    volta completa. 
%}
r_circ = wn;                            % Raio do círculo
theta_circ = linspace(0, 2*pi, 200);    % Variação angular para definir o círculo
x_circ = r_circ * cos(theta_circ);
y_circ = r_circ * sin(theta_circ);

% Plotando cícrulo para wn = 1.9rad/s
plot(x_circ, y_circ, 'b--', 'LineWidth', 0.5); %Círculo azul tracejado

% Adicionando ζ_A = 0.97 e wn_A = 2.11 rad/s ao mesmo gráfico
%{
    Obs: O procedimento a seguir é semelhante ao anterior, com diferença
    apenas nos valores de wn e ζ. Sendo assim, também valem os mesmo comentários
%}
beta_A = acos(zeta_A);
beta_A_graus = 180*beta_A/pi;
theta_reta_A = pi - beta_A;
theta_reta_graus_A = 180*theta_reta_A/pi;

% Criando a reta a partir da origem
x_A = r * cos(theta_reta_A);                  %Conversão coord polares
y_A = r * sin(theta_reta_A);

% Plotando reta para ζ1 = 0.97
plot(x_A, y_A, 'm--', 'LineWidth', 0.5);      %Reta tracejada magenta

% Criando círculo de raio 2.11 para wn1 = 2.11rad/s
r_circ_A = wn_A;                              %Raio do círculo
x_circ_A = r_circ_A * cos(theta_circ);
y_circ_A = r_circ_A * sin(theta_circ);

% Plotando cícrulo para wn1 = 2.11rad/s
plot(x_circ_A, y_circ_A, 'g--', 'LineWidth', 0.5); %Círculo verde tracejado

% Adicionando uma legenda manual pois legend('texto') não está rodando com rlocus
h1 = plot(NaN, NaN, 'k', 'LineWidth', 1.5);                     %Linha preta para o Lugar das Raízes
h2 = plot(NaN, NaN, 'rx', 'MarkerSize', 7.5, 'LineWidth', 1);   %Polos (x vermelho)
h3 = plot(NaN, NaN, 'ro', 'MarkerSize', 7.5, 'LineWidth', 1);   %Zeros (o azul)
h4 = plot(NaN, NaN, 'c--', 'LineWidth', 1.5);                   %Reta tracejada ciano para ζ = 0.5
h5 = plot(NaN, NaN, 'm--', 'LineWidth', 1.5);                   %Reta tracejado magenta para ζ = 0.21
h6 = plot(NaN, NaN, 'b--', 'LineWidth', 1.5);                   %Círculo tracejado azul para wn = 1.9
h7 = plot(NaN, NaN, 'g--', 'LineWidth', 1.5);                   %Círculo tracejado verde para wn1 = 2.11

%Comando legend('') usando as legendas manuais definidas anteriormente
legend([h1, h2, h3, h4, h5, h6, h7], 'Lugar das Raízes de G(s)', 'Polos de A(s)', 'Zeros de A(s)', 'Reta ζ = 0.5', 'Reta ζ_A = 0.21', 'Círculo ω_n = 1.9 rad/s', 'Círculo ω_n_A = 2.11 rad/s');

hold off                                                        %Finalizando o hold


%% Letra B) Aplicação de Degrau unitário 

% Definindo vetor de tempo
t = (0:0.01:50);

% Aplicando degrau usando função step()
resp_degrau = step(A,t);              %Armazenando na variável resp_degrau o degrau aplicado

% Plotando o gráfico com o degrau aplicado usando o vetor de tempo t definido
figure                                %Iniciando nova figura
plot(t, resp_degrau)                  %Plotando gráfico

%configurações do gráfico
grid on
xlabel('Tempo t(s)')
ylabel('Amplitude')
title('Resposta no tempo do sistema não compensado ao degrau unitário')

%% Letra C) Cálculo do sobressinal percentual, da frequência 𝑓 e de 𝜔𝑑
% Para encontrar o sobressinal, vamos varrer o vetor 'degrau_ref' e encontrar seu valor máximo usando a função max() do MatLab.

Mp = max(resp_degrau);                      %Mp é o sobressinal, valor máximo de amplitude ao degrau

% Cáculando sobressinal percentual usando 100*(c(tp) - c(∞)/ c(∞))
Mp_percentual = 100* (Mp - resp_degrau(end)/resp_degrau(end));    %Resp_degrau(end) para encontrar c(∞) e aplicar a fórmula anterior

% Exibindo na Command Window os valores de Mp e Mp_Percentual
fprintf('          Sobressinal             \n');

fprintf('\nMp = %.8f\n', Mp);
fprintf('Sobressinal Percentual = %.5f%%\n', Mp_percentual);
fprintf('_________________________________________\n');

% Cálculo de 𝑓 e de 𝜔d
%{
  Para encontrarmos f e wd vamos calcular o período entre picos; para isso é preciso encontrar os picos da resposta ao 
  degrau e depois calcular o período, o que é feito a seguir:
%}

%Armazenando os picos
[picos, T_p] = findpeaks(resp_degrau,t);

%Calculando o período T
T = (T_p(2) - T_p(1));

%Calculando frequência f e frequência amortecida wd
f = 1/T; 
wd = 2*pi*f;

% Exibindo na Command Window os valores de f e wd
fprintf('          f e wd             \n');

fprintf('\nf = %.8f Hz\n', f);   
fprintf('wd = %.8f rad/s\n', wd); 
fprintf('_________________________________________\n');

% Comparando os valores com polos dominantes
%Definindo os polos dominantes de A(s); Mais próximos ao eixo imaginário, parte Re -> 0
polos_complexos = polos_A(imag(polos_A) ~= 0);  %Salvando em um vetor todos os polos com parte imagninaria !=0

% Encontrando o polo dominante (menor parte real)
%{
   Obs: o procedimento seguinte vale para vários polos, nesse caso como só temos um par de polos complexos conjugados, 
   já poderíamos atribuir qual seria o dominante
%}

%Salvando o índice do polo dominante
[~, idx_dominante] = min(abs(real(polos_complexos)));

%Usando o indice no vetor polos_complexos definimos o polo_dominante
polo_dominante = polos_complexos(idx_dominante);      

%{ 
    COMENTÁRIO:
    Os cálculos manuais dão resultados parecidos com as fórmulas de sistemas de segunda ordem, mas a diferença acontece porque o nosso sistema,
    como visto na Letra A, é de terceira ordem. Mesmo com ganho unitário, ele tem dois polos complexos conjugados e um polo real,
    o que o torna essencialmente de terceira ordem. Os resultados ficam próximos porque os polos complexos conjugados têm maior influência no comportamento do sistema,
    fazendo com que ele se pareça com um de segunda ordem. Mas a precisão é limitada pelas próprias características de um sistema de terceira ordem.
    A parte imaginária dos polos dominantes e a frequência amortecida tem relação direta, ou seja, são equivalentes. 
    Temos a realação para s: s = σ +- jwd, sendo assim wd = parte imaginária do polo dominante. 
    O que pode ser comprovado através do erro aprx = 0 calculado a seguir:
%}
erro_wd_polo_dom = abs(imag(polo_dominante)) - wd;

%% Letra D) Aplicação de Rampa

% Criando novo vetor de tempo para aplicar a rampa
t_ramp = (0:0.01:200);

% Para aplicar a rampa aqui no Matlab vamos usar o comando step(FT) mas aplicado a (1/s)*A
resp_rampa = step((1/s)*A,t_ramp);

% Plotando o gráfico
figure                                     %nova figura
plot(t_ramp, resp_rampa, 'k') %plotando em preto
hold on
plot(t_ramp, t_ramp, '--r')
legend('Resposta do Sistema a Rampa','Rampa Unitária','Location','NorthWest');

%Configurações do gráfico
grid on
xlabel('Tempo t(s)')
ylabel('Amplitude')
title('Resposta no tempo do sistema não compensado a Rampa')
legend('Resposta do Sistema a Rampa','Rampa Unitária','Location','NorthWest');

% Cálculo da Constante de velocidade 𝐾𝑣 para 𝐺(𝑠) e Erro de Regime Permanente

%Temos que Kv = limite s->0 (sG(s)), sendo assim:
%Definindo a variável simbólica e FTMA em u para calcular o limite
syms u;
[num_G, den_G] = tfdata(G, 'v');                %Obtendo os coeficientes do numerador e denominador de G
G_lim = poly2sym(num_G,u)/poly2sym(den_G,u);    %Definindo G(u) para o calculo do limite, pois s ja está reservado par tf

% Calculando Kv 
Kv = limit(u*G_lim,u,0);                        %Usando função limit() para cáculo do Kv

% Calculando o Erro de Regime Permanente
erro_regime_perm = double(1/Kv);

% Gráfico e erro observado
resp_rampa_200 = resp_rampa(end);

% Temos que erro_observado = entrada - saida
erro_obs = t_ramp(end) - resp_rampa_200;
erro_relativo = erro_obs - erro_regime_perm;

%{
    Ao calcular o limite de s tendendo ao infinito de s*G(FTMA), obtem-se o
    valor de 4.53, o que corresponde ao erro obsevado entre a saída para
    entrada em rampa no gráfico com a referência. O erro relativo demonstra a diferença entre o 
    erro de regime permanente calculado e o observado, que é baixo.

    Observa-se que as duas formas de calcular o erro chegam no mesmo valor
    quando executado o código. Primeiramente, o erro de regime permanente é
    calculado usando o valor de Kv encontrado através do limite. É possível
    chegar ao mesmo valor de erro utilizando os pontos do gráfico. Com o
    compensador avanço-atraso, pretende-se diminuir esse erro 7 vezes
%}

% Exibindo na Command Window o valor de Kv e erro
fprintf('              Kv e Erros               \n');
fprintf('\nKv = %.4f\n', Kv);
fprintf('Erro de Regime Permanente = %.8f\n', erro_regime_perm);
fprintf('Erro Observado = %.8f\n', erro_obs);
fprintf('Erro Relativo = observ. - reg. perm = %.8f\n', erro_relativo);
fprintf('_________________________________________\n');


%% Letra E) Projetando um compensador avanço-atraso

%Cálculo do polo desejado a partir das especificações do projeto

%Primeiramente, é calculada a parte real do polo
PoloDesejadoReal = -(wn*zeta);
%Depois, calcula-se a parte imaginária do polo
PoloDesejadoImg = wn*sqrt(1-zeta^2);
%{
    Com isso, é possível escrever os polos complexos desejados, composto pela
    parte real e pela parte imaginária (positiva e negativa)
%}

PoloDesejado1 = (PoloDesejadoReal) + (PoloDesejadoImg)*1i;
PoloDesejado2 = (PoloDesejadoReal) - (PoloDesejadoImg)*1i;
%{
    Também é necessário calcular a deficiência angular. É utilizada a
    função angle() para encontrar o ângulo da diferença do polo
    desejado em relção aos polos e zeros do sistema. A função rad2deg() é
    utilizada para transformar o ângulo de radianos para graus.
%}

DeficienciaAngular = 180 - rad2deg(angle(PoloDesejado1 - polos_G(1))) - rad2deg(angle(PoloDesejado1 - polos_G(2))) - rad2deg(angle(PoloDesejado1 - polos_G(3))) + rad2deg(angle(PoloDesejado1 - zeros_G(1))) + rad2deg(angle(PoloDesejado1 - zeros_G(2)));

%A deficiência angular é mostrada no Command Window
disp(['A deficiência angular é: ', num2str(DeficienciaAngular)]); 
disp('------------------------');

%{
Primeiramente, precisamos encontrar a posição do zero e do polo da parte do
avanço do compensador. Para isso, é uitlizado um for, e o zero é
posicionado entre -4 e -0.5 no eixo real
%}
for ZeroAvanco = 0.5:0.05:4
    %{
        Como o for apresenta valores positivos, é necessário inverter o valor
        da posição do zero do avanço
    %}
    ZeroAvanco = -ZeroAvanco;
    %{
        É necessário encontrar o valor do ângulo Alfa. Ele é o ângulo da
        diferença entre o polo desejado e o zero do avanço
    %}
    Alfa = rad2deg(angle(PoloDesejado1 - ZeroAvanco)); %Novamente a função rad2deg() é utilizada
    %{
        Para que seja possível encontrar a posição do polo do avanço, é
        necessário que o ângulo Alfa seja maior que o módulo do ângulo da
        deficiência angular. Caso essa condição não seja satisfeita, é
        necessário utilizar duas redes por avanço de fase, cada uma
        contribuindo com a metade do ângulo de avanço de fase.
    %}
    if (Alfa > -DeficienciaAngular)
        %Se a condição for satisfeita, é possível calcular a posição do polo do avanço
        PoloAvanco = -(PoloDesejadoImg - (PoloDesejadoReal*tand(Alfa+DeficienciaAngular)))/(tand(Alfa+DeficienciaAngular));
        %{
            Com o  zero e o polo do avanço, é possível verificar se a parte do
            avanço do compensador está contribuindo com o ângulo correto.
            Com isso, a contribuição angular do avanço de fase deve ser o oposto da deficiência angular
        %}
        Teste = rad2deg(angle(PoloDesejado1 - ZeroAvanco)) - rad2deg(angle(PoloDesejado1 - PoloAvanco));

        %Também é possível encontrar os valores de T1 e Gama 
        T1 = -1/(ZeroAvanco); %O zero do avanço é utilizado para o cálculo de T1
        Gama = -PoloAvanco*T1; %O valor encontrado de T1 é utilizado para encontrar o valor de Gama
        %{
            Também é possível encontrar o valor de Kc, utilizando o polo
            desejado, a função de tranferência do sistema de controle e da
            parte do avanço do compensador
        %}
        Kc = abs(((PoloDesejado1-PoloAvanco)*PoloDesejado1*(PoloDesejado1+1.5)*(PoloDesejado1+5))/((PoloDesejado1 - ZeroAvanco)*(2*(PoloDesejado1)^2 + 4*PoloDesejado1 + 34)));
        %{
            O Beta também é calculado para que o erro de regime permanente seja
            reduzido em 7 vezes. Para isso, o valor do novo Kv deve ser 7 vezes
            o valor do Kv antigo (sistema de controle sem o compensador)
        %}
        Beta = (7*Gama)/Kc;
    else
        %{
            Quando Alfa não é maior que o módulo do ângulo da deficiência
            angular, é necessário utilizar duas redes por avanço de fase, onde cada
            uma realiza a metade do ângulo de avanço de fase
        %}
        DeficienciaAngularMetade = DeficienciaAngular/2; %A Deficiência angular é dividida pela metade
        
        %Com isso, é calculado a posição do polo do avanço
        PoloAvanco = -(PoloDesejadoImg - (PoloDesejadoReal*tand(Alfa+DeficienciaAngularMetade)))/(tand(Alfa+DeficienciaAngularMetade));
        %{
            Um teste é realizado para verificar se a contribuição do ângulo de
            avanço de fase está correta. Como cada rede realiza metade da
            contribuição, é necessário realizar a multiplicação por 2.
        %}
        Teste = 2*(rad2deg(angle(PoloDesejado1 - ZeroAvanco)) - rad2deg(angle(PoloDesejado1 - PoloAvanco)));
        
        %Também é possível encontrar os valores de T1 e Gama
        T1 = -1/(ZeroAvanco); %O zero do avanço é utilizado para o cálculo de T1
        Gama = -PoloAvanco*T1; %O valor encontrado de T1 é utilizado para encontrar o valor de Gama
        %{
            Também é possível encontrar o valor de Kc, utilizando o polo
            desejado, a função de tranferência do sistema de controle e da
            parte do avanço do compensador. Como existem duas redes por avanço
            de fase, o cálculo do Kc é diferente do realizado anteriormente (A
            função de transferência relacionada a parte do avanço está elevada
            ao quadrado, devido ao uso das duas redes de avanço de fase)
        %}
        Kc = abs((((PoloDesejado1-PoloAvanco)^2)*PoloDesejado1*(PoloDesejado1+1.5)*(PoloDesejado1+5))/(((PoloDesejado1 - ZeroAvanco)^2)*(2*(PoloDesejado1)^2 + 4*PoloDesejado1 + 34)));
       
        %O cálculo de Beta também é diferente devido ao uso das duas redes
        Beta = (7*(Gama^2))/Kc;
    end
    %{
        Depois disso, é possível encontrar os valores de T2 que atendem as
        especificações do projeto. Para isso, é realizado um segundo for, onde são testados valores de T2 entre 5 e 15 
    %}
    for T2 = 5:0.1:15
        %Primeiramente, é calculado o módulo do compensador atraso
        ModuloAtraso = abs((PoloDesejado1 + (1/T2))/(PoloDesejado1 + (1)/(Beta*T2))); %Utiliza-se o valor de Beta calculado anteriormente
        
        %Depois, a contribuição angular do compensador atraso é calculada
        ContribuicaoAngularAtraso = rad2deg(angle((PoloDesejado1 + (1/T2))/(PoloDesejado1 + (1)/(Beta*T2)))); %O valor de Beta também é utilizado nessa conta
        %{
            Como especificado no projeto, o módulo do compensador atraso pode
            ter módulo de até 0.96 e contribuição angular de até 5°. Com isso,
            é feito um if para verificar essa condições
        %}
        if (ModuloAtraso >= 0.96) && (ModuloAtraso <= 1) && (ContribuicaoAngularAtraso >= -5) && (ContribuicaoAngularAtraso <= 0)
            %{
                Como explicado anteriormente a função de tranferência do
                compensador avanço-atraso pode ter duas formas diferentes. Com
                isso, é necessário realizar um outro if.
            %}
            if (Alfa > -DeficienciaAngular)
                Gc = Kc*((s+(1/T1))/(s+(Gama/T1)))*((s+(1/T2))/(s+(1/(Beta*T2)))); %Neste caso, é utilizada apenas uma rede por avanço de fase
            else
                %{
                    Quando são utilizadas duas redes, a parte da função de
                    transferência relacionada ao avanço de fase é elevada ao
                    quadrado
                %}
                Gc = Kc*(((s+(1/T1))/(s+(Gama/T1)))^2)*((s+(1/T2))/(s+(1/(Beta*T2)))); 
            end
            %Com Gc, é possível encontrar a FTMF
            Gcomp = G*Gc/(1 + G*Gc);
            
            %Definindo vetor de tempo
            t = 0:0.1:50;
            
            %Aplicando degrau usando função step()
            RespostaDegrau = step(Gcomp, t);
            
            %Criando uma variável a para cálculo do tempo de assentamento
            a = 501;
            %{
                O código percorre o vetor de trás para frente, procurando o 
                último instante em que a resposta y saiu da faixa de 
                tolerância de 2% em torno do valor final
            %}
            while (RespostaDegrau(a) > 0.98) && (RespostaDegrau(a) < 1.02) %No projeto, é utilizada uma faixa de 2%
                a =  a - 1;
            end       
            %Com isso, é possível encontrar o tempo de assentamento
            ts = a*0.1;
            %Também é possível encontrar o sobressinal utilizando a função max()
            m = max(RespostaDegrau);
            %{
                Para o projeto, é pedido um sobressinal menor que 18% e maior
                que 10%. Além disso, o tempo de assentamento deve ser menor ou
                igual a 3.2s. Para verificar essas condições, um if é criado.
            %}
            if (ts <= 3.2) && (m < 1.18) && (m > 1.1)
                %Atendida essas condições, os valores são armazenados na matriz de soluções
                k = k + 1; %A variável k é um contador de soluções
                
                %Matriz que armazena as soluções
                solution(k,:) = [ZeroAvanco PoloAvanco Teste T1 Gama Kc Beta ts m Alfa a];
                %{
                    Os valores que atendem as condições do projeto também são
                    exibidos no command window
                %}
                disp('Zero do Avanço:');
                disp(ZeroAvanco);
                disp('Polo do Avanço:');
                disp(PoloAvanco);
                disp('Teste da Deficiência angular:');
                disp(Teste);
                disp('T1:');
                disp(T1);
                disp('Gama:');
                disp(Gama);
                disp('Kc:');
                disp(Kc);
                disp('Beta:');
                disp(Beta);
                disp('Tempo de assentamento:');
                disp(ts);
                disp('Sobressinal:');
                disp(m);
                disp('------------------------');
            end
        end
    end
end
%{
    A matriz de soluções é colocada em ordem crescente de acordo com o
    sobressinal (Do menor para o maior sobressinal)
%}
sortsolution = sortrows(solution, 9);

%A matriz de soluções ordenada é exibida no command window
disp(sortsolution);

%O número de soluções também é mostrada
disp(['Número de soluções: ', num2str(k)]);
disp('------------------------');

%A primeira linhha da matriz, ou seja, aquela que apresenta menor sobressinal, é exibida no command window 
disp(sortsolution(1,1:11));
disp('------------------------');
%{
    Para traçar os gráficos, escolhemos a solução que apresenta o menor
    sobressinal que atenda as especificações do projeto. Para isso, os valores
    referentes a essa solução são coletados na matriz (Como a matriz está
    ordenada, os valores se encontram na primeira linha)
%}
ZeroAvanco = sortsolution(1,1); PoloAvanco = sortsolution(1,2); Teste = sortsolution(1,3); T1 = sortsolution(1,4);
Gama = sortsolution(1,5); Kc = sortsolution(1,6); Beta = sortsolution(1,7); ts = sortsolution(1,8); m = sortsolution(1,9);
Alfa = sortsolution(1,10); a = sortsolution(1,11);
%{
    Como explicado anteriormente, o compensador pode ter uma ou duas redes por
    avanço de fase. Com isso, é necessário realizar novamente o if.
%}
if (Alfa > -DeficienciaAngular)
    %Função de transferência do compensador com uma rede de avanço de fase
    Gc = Kc*((s+(1/T1))/(s+(Gama/T1)))*((s+(1/T2))/(s+(1/(Beta*T2))));
else
    %Função de tranferência do compensador com duas redes de avanço de fase
    Gc = Kc*(((s+(1/T1))/(s+(Gama/T1)))^2)*((s+(1/T2))/(s+(1/(Beta*T2))));
end
%Encontra-se a FTMF
Gcomp = G*Gc/(1 + G*Gc);

%Aplicando degrau na FTMF usando função step()
RespostaDegrau = step(Gcomp, t);

%Traçando o gráfico do lugar das raízes do sistema compensado
figure;
%{
    A função rlocus é utilizado para traçar o LR. Como temos o compensador, Gc
    também é usado para traçar o LR
%}
rlocus(G*Gc);
grid on;
hold on;

%Adicionando o título e nomeando os eixos
title('Lugar das raízes do sistema compensado');
xlabel('Real');
ylabel('Imaginário');
%{
    Aqui, os novos polos de malha fechada dominantes são destacados no gráfico.
    Esses polos serão os polos desejados no projeto (Aqueles que fornecem Wn =
    1.9 rad/s e Amortecimento = 0.5)
%}
plot(real([PoloDesejado1 PoloDesejado2]), imag([PoloDesejado1 PoloDesejado2]), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % 'rx' para marcar com 'X' vermelho

%Adição da legenda ao gráfico do LR
h1 = plot(NaN, NaN, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend(h1, 'Novos polos de malha fechada dominantes');
hold off;
%{
    Traçando o gráfico da resposta ao degrau unitário até 50s (O vetor t é
    novamente utilizado)
%}
figure;
%Plotando o gráfico utilizando a função plot
plot(t, RespostaDegrau, '-');
grid on;

%Adicionando o título e nomeando os eixos
title('Resposta no tempo do sistema compensado ao degrau unitário');
xlabel('Tempo t(s)');
ylabel('Amplitude de c(t)');
hold on;

%Nessa parte, o sobressinal é destacado no gráfico
%{
    O valor no eixo y é o valor máximo da resposta ao degrau. Para descobrir o
    valor no eixo x desse ponto, é realizada a seguinte lógica:
%}
[m, idxMax] = max(RespostaDegrau); 
x_max = t(idxMax);
plot(x_max, m, 'r.', 'MarkerSize', 15); %Plotando o valor do sobressinal

%Aqui, o Tempo de assentamento é destacado no gráfico
plot(ts, RespostaDegrau(a+1), 'g.', 'MarkerSize', 15);
%{
    Adicionando uma legenda manual para o sobressinal e o tempo de
    assentamento
%}
h1 = plot(NaN, NaN, 'r.', 'MarkerSize', 15);
h2 = plot(NaN, NaN, 'g.', 'MarkerSize', 15);

%Comando legend() usando as legendas manuais definidas anteriormente
legend([h1, h2], 'Sobressinal', 'Tempo de Assentamento');

%Criando um vetor de tempo t1 para traçar o gráfico da respost a rampa
t_rampa = 0:0.1:200;
%{
    Para aplicar a rampa aqui no Matlab vamos usar o comando step(FT) mas
    aplicado a (1/s)*G
%}
rampa = step (Gcomp*(1/s), t_rampa);

%Plotando o gráfico
figure;
plot(t_rampa, rampa, 'b'); %Resposta a rampa unitário do sistema compensado
grid on;
hold on;
plot(t_rampa,t_rampa,'--r'); %Rampa unitária
%Configuração do gráfico
legend('Resposta do Sistema Compensado a Rampa','Rampa Unitária','Location','NorthWest');
xlabel('Tempo t(s)');
ylabel('Amplitude');
title('Resposta no tempo do sistema compensado a Rampa');
legend('Resposta do Sistema Compensado a Rampa','Rampa Unitária','Location','NorthWest');
%{
    Calculo do erro. O erro será a diferença entre o  valor da rampa unitária
    em t = 200s, que é 200, e o valor da resposta do sistema compensado a rampa
    em t = 200s
%}
Erro = 200 - rampa(2001);

%O erro é exibido no command window
disp('O erro de regime permanente do sistema compensado é:');
disp(Erro);
%{
    É possível perceber, com a execução do código, que o erro diminuiu 7 vezes
    com a implementação do compensador
%}






