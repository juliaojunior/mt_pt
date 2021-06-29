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
LL = 100;

%velocidade da luz
c = 299792458;
mu0 = 4*pi*1e-7;

% indice de refracao do background
n1 = 1; %sqrt(epsilon(1))*sqrt(mu(1));
epsa = 1;  % epsilon do background


%angulo de incidencia
%(em radianos)
Ai = 5*(pi/180);





%% Definicao das propriedades das camadas


% quantidade de bicamadas
Nlay = 1;

%Tamamnho das camadas
d = 125e-6; %125 micrometros


% mu da camada 
muA = 1.0;    % não é magnético

% epsilon das camadas
e1 = 0.0001;
e2 = 0.001;
epsg = e1 - 1j*e2;  % epsilon com ganho
epsp = e1 + 1j*e2;  % epsilon com perda


%% Definicao das variaveis de loop


% frequencia
%omega = 3e13;
%omega = 1e13*linspace(0.0,6.0,LL);
omega = 1e13*linspace(0.0,6.0,LL);


R_ri = zeros(LL,1);  % reflexão direita
R_le = zeros(LL,1);  % reflexão esquerda
T = zeros(LL,1);     % transmissão
freq_lattice = zeros(LL,1);  




%% loop para transmissao

Trans = 0;
freqY = 1;
for f = omega
    
    % vetor de onda
    k0 = f/c;
    
    % vetor de onda longitudinal
    kza = k0*sqrt(epsa)*cos(Ai);
    
    MA = mt1(f,muA,epsg,Ai,d,n1);
    MB = mt1(f,muA,epsp,Ai,d,n1);
    
    M = (MA*MB)^Nlay;

    %t = 2/(M(1,1) + (c/f)*(kza/epsa));
    t = 2/(M(1,1) + (c/f)*(kza/epsa)*M(1,2)+(f/c)*(epsa/kza)*M(2,1)+M(2,2));
    Trans = t.^2;
    
    T(freqY) = Trans;
    freqY = freqY+1;
    
end
freqY = 1;




%% plotagem das imagens

plot(omega,abs(T))


