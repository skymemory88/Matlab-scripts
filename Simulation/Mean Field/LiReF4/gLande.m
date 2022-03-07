function [gLande] = gLande(L,S)

J = L + S;

gLande = (3/2)+((S.*(S+1)-L.*(L+1))./(2.*J.*(J+1)));

end