pgm1.m : erreurs d'integration des HS pour differents n et m

pgm2.m : e.i. pour les HS et differents N, pour le probleme d'optimisation (eps_weights.m)

pgm4.m : e.i. pour un N donne, un degre n de HS donne, et une methode donnee (commenter/decommenter)

pgm5.m : e.i. pour f1,f2,f3 pour avec le probleme d'optimisation (eps_weights.m)

pgm6.m : calcul des poids par l'interpolation sur toute la CS

pgm9.m : calcul des poids par la methode avec epsilon, et recherche du meilleur p_opt (diverses fonctions ecrites en .m utilisées)

pgm11.m : methode de Freeden, p.1201, part 2.2., pour la grille CS

pgm12.m : courbes d'e.i. pour les HS (n,m donnes) avec les methodes avec et sansles epsilon, en fonction de N

pgm13.m : tests des fonctions f_X de Alouges, methode avec epsilon

pgm14.m : allure du maximum des |epsilon| en fonction de N, methode avec epsilon

pgm15.m : interpolation sur les fonctions f_X de Alouges

pgm16.m : extrapolation de Richardson, formule sans epsilon

card_sin.m : fonction sinus cardinal reecrite

compute_A.m (parametre n_max, retourne la matrice A et le vecteur err_i des e.i. pour la formule sans les epsilon) : matrice A d'interpolation avec les poids d'un panel, sans les symetries

compute_A_sym.m (parametre n_max, retourne la matrice A et le vecteur err_i des e.i. pour la formule sans les epsilon) : idem mais avec les symetries.

compute_A_sym_alouges.m (parametres rayon, nt, nl, retourne la matrice A et le vecteur err_i des e.i. pour la formule sans les epsilon) : matrice d'interpolation sur les fonctions f_X

create_weights.m (parametres N, choix entre 1 et 4, retourne une matrice (N+1)*(N+1) des poids) : creation de poids avec des proprietes de symetries differentes

dist_to_set.m : utilise pour pgm11.m

eps_weights.m (parametre n_max, retourne la matrice des poids et la valeur optimale obtenue) : calcul les poids par optimisation

fun_alouges.m, fun_cst.m, fun_f1.m, fun_f2.m, fun_f3.m, fun_f4.m, fun_f5.m : valeurs des differentes fonctions sur la CS

int_funs_fornberg.m (parametre la matrice des poids, retourne le vecteur des e.i.) : vecteurs des e.i. pour f0,f1,f2,f3,f4 avec les poids choisis

int_sum_sph.m (parametres les poids, n_max, retourne une valeur scalaire) : fonction utilisée pour eps_weights.m

int_weights.m (parametres les poids, et les valeurs de la fonction sur les 6 faces, retourne l'integrale totale et les integrales sur les 6 faces) : integrale numerique sur la CS pour une fonction et des poids choisis

make_cs_grid.m (parametre N) : mod74.m avec N en parametre ; initialise la grille CS et ses differents parametres

plot_cs1_mesh.m (parametres N-1, N+1) : affiche la CS simple

plot_cs2_mesh.m : affiche l'harmonique spherique (n,m donnes) sur la CS

solve_weights.m (parametres A, err_i, k nombre d'HS a considerer, err = 1 si on veut les epsilon; = 0 sinon, sym = 1 si on veut les symetries; = 0 sinon, retourne la matrice des epsilon (ou des poids si err=0)) : resout le probleme lineaire AW=b pour la methode avec epsilon

solve_weights_alouges.m : idem mais avec interpolation sur les f_X

sph.m (parametres n,m degre et ordre, x,y,z les matrices des coordonnees, retourne la matrice des valeurs) : calcule la matrice des valeurs de l'HS (n,m) en (x,y,z)

sph_rot.m : idem avec rotation possible, teta0 et lambda0 en parametres

view_int_sph.m (parametres n_max, les poids, detail = 1 pour juste I; = 2 pour I1, I2, ..., I6, sph4=1 pour afficher par pas de 4; =0 pour afficher toutes les HS) : affiche toutes les e.i. des HS jusqu'a n_max, avec les poids choisis, selon le niveau de detail choisi
