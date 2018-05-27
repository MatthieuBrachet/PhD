** Mar-25-2018
Programmation p6.m des points de la CS de facon sequentielle.
Dábbord les 8 points du cube. Puis
les 12(N-1) points des aretes. Enfin les points internes.
Comparaison avec code direct pour N=1, N=2 OK.
** Mar-26-2018
p7.m: comparaison des matrices aan et aa=[aa0;aa1,aa2....aamax].
    OK pour N=2. 
Suite p8.m: j'enleve l'ancien code. Tests dans p8 des rangs avec 
quelques exemples
N=1,2,3,4,6,8
nhs_max=0,1,2,3,4,5
relevé du rang dans chaque cas.
Une loi simple semble se dessiner.
MAIS probleme quand :
1) on met A'*A et calcul du rang. 
Ca ne marche pas !!!!
2) calcul du rang sensible a TOL quand la taille de la matrice augmente.
Donc la loi n''est pas si claire.... Ca a commencé a 
ne plus marcher a partir de N=6. Pour N=5 ca allait, cf relevé papier.
Reprendre le calcul du rang. ET structure de la matrice. Voir
le calcul du rang sur la plus petite matrice possible A'*A ou A* A'. A REVOIR.
C'est p8.m le code actuel 
