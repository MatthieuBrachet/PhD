function [v_fI,v_fII,v_fIII,v_fIV,v_fV,v_fVI]=ds(u_fI,u_fII,u_fIII,u_fIV,u_fV,u_fVI,n,nn)

 % - CALCUL 1/2 SOMME SUR LES 12 ARETES
% -----------------------------------
lwk=zeros(n);
v_fI=u_fI;v_fII=u_fII;v_fIII=u_fIII;v_fIV=u_fIV;v_fV=u_fV;v_fVI=u_fVI;
% I-II
lwk=0.5*(u_fI(nn,2:nn-1)+u_fII(1,2:nn-1));
v_fI(nn,2:nn-1)=lwk; 
v_fII(1,2:nn-1)=lwk;
% II-III
lwk=0.5*(u_fII(nn,2:nn-1)+u_fIII(1,2:nn-1));
v_fII(nn,2:nn-1)=lwk; 
v_fIII(1,2:nn-1)=lwk;
% % III-IV
lwk=0.5*(u_fIII(nn,2:nn-1)+u_fIV(1,2:nn-1));
v_fIII(nn,2:nn-1)=lwk; 
v_fIV(1,2:nn-1)=lwk;
% IV-I
lwk=0.5*(u_fIV(nn,2:nn-1)+u_fI(1,2:nn-1));
v_fIV(nn,2:nn-1)=lwk; 
v_fI(1,2:nn-1)=lwk;
% I-V
lwk=0.5*(u_fI(2:nn-1,nn)+u_fV(2:nn-1,1));
v_fI(2:nn-1,nn)=lwk; 
v_fV(2:nn-1,1)=lwk;
% II-V
for i=2:nn-1,
 lwk(i-1)=0.5*(u_fII(i,nn)+u_fV(nn,i));
end
for i=2:nn-1,
 v_fII(i,nn)=lwk(i-1);
 v_fV(nn,i)=lwk(i-1);
end
% III-V
for i=2:nn-1,
 lwk(i-1)=0.5*(u_fIII(i,nn)+u_fV(nn+1-i,nn)); % nn+1-i varie de nn-1 a 2, pas -1.
end
for i=2:nn-1
 v_fIII(i,nn)=lwk(i-1); 
 v_fV(nn+1-i,nn)=lwk(i-1);
end
% % IV-V
for i=2:nn-1,
 lwk(i-1)=0.5*(u_fIV(i,nn)+u_fV(1,nn+1-i));
end
for i=2:nn-1,
 v_fIV(i,nn)=lwk(i-1);
 v_fV(1,nn+1-i)=lwk(i-1);
end
% I-VI
lwk=0.5*(u_fI(2:nn-1,1)+u_fVI(2:nn-1,nn));
v_fI(2:nn-1,1)=lwk; 
v_fVI(2:nn-1,nn)=lwk;
% II-VI
for i=2:nn-1,
 lwk(i-1)=0.5*(u_fII(i,1)+u_fVI(nn,nn+1-i));
end
v_fII(2:nn-1,1)=lwk; 
v_fVI(nn,nn-1:-1:2)=lwk;
% III-VI
lwk=0.5*(u_fIII(2:nn-1,1)+u_fVI(nn-1:-1:2,1));
v_fIII(2:nn-1,1)=lwk; 
v_fVI(nn-1:-1:2,1)=lwk;
% IV-VI
for i=2:nn-1,
 lwk(i-1)=0.5*(u_fIV(i,1)+u_fVI(1,i));
end
v_fIV(2:nn-1,1)=lwk; 
v_fVI(1,2:nn-1)=lwk;
%
%
% - CALCUL 1/3 SOMME  SUR LES 8 SOMMETS
% ------------------------------------
% lgwk=zeros(3);
% I-II-V
lgwk=(1/3)*(u_fI(nn,nn)+u_fII(1,nn)+u_fV(nn,1));
v_fI(nn,nn)=lgwk;
v_fII(1,nn)=lgwk;
v_fV(nn,1)=lgwk;
% II-III-V
lgwk=(1/3)*(u_fII(nn,nn)+u_fIII(1,nn)+u_fV(nn,nn));
v_fII(nn,nn)=lgwk;
v_fIII(1,nn)=lgwk;
v_fV(nn,nn)=lgwk;
% III-IV-V
lgwk=(1/3)*(u_fIII(nn,nn)+u_fIV(1,nn)+u_fV(1,nn));
v_fIII(nn,nn)=lgwk;
v_fIV(1,nn)=lgwk;
v_fV(1,nn)=lgwk;
% IV-I-V
lgwk=(1/3)*(u_fIV(nn,nn)+u_fI(1,nn)+u_fV(1,1));
v_fIV(nn,nn)=lgwk;
v_fI(1,nn)=lgwk;
v_fV(1,1)=lgwk;
% I-II-VI
lgwk=(1/3)*(u_fI(nn,1)+u_fII(1,1)+u_fVI(nn,nn));
v_fI(nn,1)=lgwk;
v_fII(1,1)=lgwk;
v_fVI(nn,nn)=lgwk;
% II-III-VI
lgwk=(1/3)*(u_fII(nn,1)+u_fIII(1,1)+u_fVI(nn,1));
v_fII(nn,1)=lgwk;
v_fIII(1,1)=lgwk;
v_fVI(nn,1)=lgwk;
% III-IV-VI
lgwk=(1/3)*(u_fIII(nn,1)+u_fIV(1,1)+u_fVI(1,1));
v_fIII(nn,1)=lgwk;
v_fIV(1,1)=lgwk;
v_fVI(1,1)=lgwk;
% IV-I-VI
lgwk=(1/3)*(u_fIV(nn,1)+u_fI(1,1)+u_fVI(1,nn));
v_fIV(nn,1)=lgwk;
v_fI(1,1)=lgwk;
v_fVI(1,nn)=lgwk;

