clc; clear all; close all
x=-5:0.01:5;
[X,Y]=meshgrid(x,x);
lambda=X+1i.*Y;

k1= lambda;
k2= lambda.* ( 1 + 1/5.*        k1);
k3= lambda.* ( 1 + 3/40.*       k1 + 9/40       * k2);
k4= lambda.* ( 1 + 44/45 *      k1 - 56/15      * k2 + 32/9       * k3);
k5= lambda.* ( 1 + 19372/6561 * k1 - 25360/2187 * k2 + 64448/6561 * k3  - 212/729 * k4);
k6= lambda.* ( 1 + 9017/3168  * k1 - 355/33     * k2 - 46732/5247 * k3  + 49/176  * k4 - 5103/18656 * k5);

stab=1+35/384*k1+500/1113*k3+125/192*k4-2187/6784*k5+11/84*k6;


area=(abs(stab)<=1);

figure(1)
contourf(X,Y,area)


figure(2)
contourf(X,Y,abs(stab),10)
colorbar