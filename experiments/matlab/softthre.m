function out=softthre(a,tau)
out = sign(a).* max( abs(a) - tau, 0);
end