function [ MatTr ] = mt2( nn,nnc, kk, L )
%c = 299792458;      %Velocidade da luz no vácuo


Delta = nn*(nn*kk)*L/2;
frag = nnc*sin(Delta)*(cos(conj(Delta)));
frag2 = nn*sin(Delta)*(cos(conj(Delta)));


alpha = ((abs(cos(Delta))^2)/2) - (abs(sin(Delta))^2)*nnc/2/nn;
beta = (1/2/(abs(nn)^2)) * (frag + conj(frag));
gamma = (1/2)*(frag2 + conj(frag2));

a = (alpha + conj(alpha)) + 1j*(beta + gamma);
b = -1j*(alpha - conj(alpha)) + (gamma - beta);
c = 1j*(alpha - conj(alpha)) + (gamma - beta);



Ma11 = conj(a);
Ma12 = 1j*b;
Ma21 = -1j*c;
Ma22 = a;

MatTr = [Ma11 Ma12; Ma21 Ma22];

end