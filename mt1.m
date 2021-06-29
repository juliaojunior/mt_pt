function [ MatTr ] = mt1( omega,muA,epsilonA,theta,deltaA,nn )
c = 299792458;      %Velocidade da luz no vácuo

%alphaA = omega/c*sqrt(muA*epsilonA - sin(theta)*sin(theta));
%qA = alphaA*c/omega/muA;


beta = nn*sin(theta);
Q = omega/c*sqrt(muA*epsilonA - beta^2);
P = Q*c/omega/muA;


Ma11 = cos(Q*deltaA);
Ma12 = sin(Q*deltaA);
Ma21 = Ma12;
Ma12 = -1i/P*Ma12;
Ma21 = -1i*P*Ma21;
Ma22 = Ma11;

MatTr = [Ma11 Ma12; Ma21 Ma22];

end