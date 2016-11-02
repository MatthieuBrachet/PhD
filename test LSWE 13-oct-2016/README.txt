** July-04-2016
Repertoire LSWEC_CS_July_04_2016_orig est le répertoire de Matthieu Brachet pour les equations SW linearisees sur la CS. 
Programme principal= p72_RK4

** July-06-2016
Suite dans Cubed-Sphere-3: c'est le rerertoire pour LSWE.
Cas test 1: solution zonal avec terme en temps de la forme exp(-sigmat).
Il faut mettre des termes sources (de forcage).
C'est le cas test 5 dans p72_RK4.m.
Test de convergence en maillage.

Conclusion: la divergence dans l''equation sur h donne 
une erreur non negligeable.
Le terme en divergence vient avec H*div(v).
H= hauteur de repos. Plus H est grand en metres, plus l''erreur
est amplifiée. Il faudrait mesurer cela en terme de Bas-Mach...
Vitesse du son ici c= sqrt(gH).
Plus H est grand, plus le Mach est faible.
Mach fluctuant= v(t)/c.

v(t)= proportionnel a u0 = constante devant la fonction
en theta . On a phi(theta)=fonction de Galewski, voir fun10.m Cas 5.
Au lieu de theta0=theta1=pi/4, il faudrait prendre plus petit. Car pi/4
coincide exactement avec l'angle du panel.

On a donc u0=1/10 dans ce cas (Mach fluctuant).

On a h(x,t) mauvais quand H est gradn, mais la vitesse reste bonne.

Ca ne semble pas etre le Mach qui decide de l'imprecision de h,
mais simplement le produit H * div(v). 

Vooir les dessins sur face I de la divergence. L''imprecision
forte se concentre sur les bords du panel.

Aussi: div80 n''est pas meilleur que div72. Donc le probleme ne vient pas
de l'interpolation par 0. Ce n'est pas mieux avec la formule (39)
= Dxi F . g^xi + Deta F . gêta.
On conserve l''ancienne divergence.

Par contre je passe a un angle (theta0,theta1)= (3*pi/16,-3*pi/16).

Reprise du calcul avec 3*pi/16. Pour que la coupure soit DANS le panel 1.


 



