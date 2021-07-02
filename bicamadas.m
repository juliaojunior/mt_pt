%-------------------------------------------------------------------------
%%----------------------- Transfer Matrix Method---------------------------
%---------------
%---------------
%----DAta: 29junho2021--------------------------------------------------
%----Versao: 1.0--------------------------------------------------------
%----Destaques:  -> O primeiro objetivo é reproduzir os resultados de --
%-------------   -> pelo menos um destes dois artigos: "Scattering -----
%-------------   -> properties ofPT-symmetriclayered periodic ----------
%-------------   -> structures" ou "Conservation relations and ---------
%-------------   -> anisotropic transmission resonances (..)" ----------
%


close all;
clc;
clear all;


%------------------------------------------------------------------------
%% DEFINE SIMULATION PARAMETERS
%------------------------------------------------------------------------

% tamanho dos vetores
LL = 1000;

%velocidade da luz
c = 299792458;
mu0 = 4*pi*1e-7;

% indice de refracao do background
n1 = 1; %sqrt(epsilon(1))*sqrt(mu(1));
epsa = 1;  % epsilon do background


%angulo de incidencia
%(em radianos)
Ai = 0*(pi/180);





%% Definicao das propriedades das camadas


% quantidade de bicamadas
Nlay = 1;

%Tamamnho das camadas
%d = 125e-6; %125 micrometros
d = 0.1;

% mu da camada 
muA = 1.0;    % não é magnético

% epsilon das camadas
e1 = 0.0001;
e2 = 0.001;
epsg = e1 - 1j*e2;  % epsilon com ganho
epsp = conj(epsg);  % epsilon com perda

n = 2 + 1j*0.2;
nc = conj(n);


%% Definicao das variaveis de loop


% frequencia
%omega = linspace(0.0,6.0,LL)*2*pi;
omega1 = linspace(0,3,LL)*2*pi*1e9; 


% tranmissão e reflexões
R_ri = zeros(LL,1);  % reflexão direita
R_le = zeros(LL,1);  % reflexão esquerda
T = zeros(LL,1);     % transmissão
freq_lattice = zeros(LL,1);  




%% loop para transmissao

Trans = 0;
freqY = 1;
for f = omega1
    
    % vetor de onda
    k0 = f/c;
    
    % vetor de onda longitudinal
    kza = k0*sqrt(epsa)*cos(Ai);
    
    
    %M = (MA*MB)^Nlay;
    M = mt2( n, nc, kza, d );

    t = 2/(M(1,1) + (c/f)*(kza/epsa)*M(1,2)+(f/c)*(epsa/kza)*M(2,1)+M(2,2));
    %Trans = t.^2;
    Trans = 1/(M(2,2));
    Re_esq = 1j*(M(1,2))/(M(2,2));
    Re_dir = -1j*(M(2,1))/(M(2,2));
    
    T(freqY) = Trans;
    R_le(freqY) = Re_esq;
    R_ri(freqY) = Re_dir;
    freqY = freqY + 1;
    
end
%freqY = 1;




%% plotagem das imagens

plot(omega1/10,(abs(T)).^2)
ylim( [ 0 1.25 ] )
xlim ( [ 0 0.9 ] * 1e9 ) 
%plot(omega1/10,(abs(T)).^2,omega1/10,(abs(R_ri)).^2)
%ylim( [ 0 1.25 ] )
%xlim ( [ 0 0.9 ] * 1e9 ) 
%plot(omega1/10,(abs(R_le)).^2)
%plot(omega1/10,(abs(R_ri)).^2)


