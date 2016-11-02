function [ht,lap] = sol_exacte(x,y,z,L,M)
% L>=M
[ht] = spharm(L,M,x,y,z);
lap=-L.*(L+1).*ht;
end