% =====================================================
% principal_stationnaire;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante stationnaire, avec conditions de
% Dirichlet homogene
% | u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

clear all

% Donnees du probleme.
nom_maillage = 'geomRect005.msh';
affichage = true; % false; %

% Lecture du maillage et affichage.
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(nom_maillage);

% Declarations des matrices EF et vecteur second membre.
KK = sparse(Nbpt,Nbpt);
MM = sparse(Nbpt,Nbpt);
FF = zeros(Nbpt,1);

% Boucle d'assemblage.
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  % Numéro Global
  i = Numtri(l,1);
  j = Numtri(l,2);
  k = Numtri(l,3);

  S1=Coorneu(i,:);
  S2=Coorneu(j,:);
  S3=Coorneu(k,:);
  
  % calcul des matrices elementaires du triangle l 
  Kel=matK_elem(S1, S2, S3, Reftri(l));
  Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  MM(i,i) = MM(i,i) + Mel(1,1);
  MM(i,j) = MM(i,j) + Mel(1,2);
  MM(i,k) = MM(i,k) + Mel(1,3);
  MM(j,j) = MM(j,j) + Mel(2,2);
  MM(k,k) = MM(k,k) + Mel(3,3);
  MM(j,k) = MM(j,k) + Mel(2,3);
  MM(k,j) = MM(k,j) + Mel(3,2);
  MM(k,i) = MM(k,i) + Mel(3,1);
  MM(j,i) = MM(j,i) + Mel(2,1);
  
  KK(i,i) = KK(i,i) + Kel(1,1);
  KK(i,j) = KK(i,j) + Kel(1,2);
  KK(i,k) = KK(i,k) + Kel(1,3);
  KK(j,j) = KK(j,j) + Kel(2,2);
  KK(k,k) = KK(k,k) + Kel(3,3);
  KK(j,k) = KK(j,k) + Kel(2,3);
  KK(k,j) = KK(k,j) + Kel(3,2);
  KK(k,i) = KK(k,i) + Kel(3,1);
  KK(j,i) = KK(j,i) + Kel(2,1);
end

% Calcul du second membre.
FF = zeros(Nbpt,1);
for h=1:Nbpt
	% utiliser la routine f.m
    FF(h) = f(Coorneu(h,1),Coorneu(h,2));
end

LL = MM*FF;

% Matrice EF compl?te.
AA = MM + KK;

[AA, LL] = elimine(AA, LL, Refneu);

% R?solution du probl?me par inversion et sortie des r?sultats
%(validation et r?solution num?rique.
UU = AA\LL;

% Visualisation de sigma et de la solution.
if affichage
    afficheSigma(Numtri, Reftri, Coorneu);
    affiche(UU, Numtri, Coorneu, sprintf('Solution stationnaire, pas de maillage h = 0.05 - %s', nom_maillage));
end

%% Calcul des erreurs L2 et H1

UU_exact = sin(2*pi*Coorneu(:,1)/9).*sin(pi*Coorneu(:,2));

%Calcul de l'erreur L2 et H1
ERRL2 = (dot(MM*(UU_exact-UU),(UU_exact-UU)))^0.5/dot(UU_exact,UU_exact);
ERRH1 = (dot(KK*(UU_exact-UU),(UU_exact-UU)))^0.5/dot(UU_exact,UU_exact);

%% Visualisation de la convergence

UU_exact = sin(2*pi*Coorneu(:,1)/9).*sin(pi*Coorneu(:,2));
   
h1 = 0.2;
HH = [h1];
for t=1:3
    HH = [HH, HH(t)/2];
end

%Liste de log(1/h) pour h(n+1) = h(n)/2 et h(0) = 0.2
HH = arrayfun(@(x) logar(x),HH);
    
% Calcul de l erreur L2
ERRL2 = (dot(MM*(UU_exact-UU),(UU_exact-UU)))^0.5/dot(UU_exact,UU_exact);
    
%On réalise le calcul pour différents pas h décroissants tels que
%h(n+1) = h(n)/2, les valeurs sont stockées dans la liste suivante. Une
%automatisation pour 5 valeurs ne semblait pas pertinente, on pourra
%vérifier que les valeurs ci-dessous sont bien celles trouvées lors de
%l'éxecution de l'algorithme.
ERRL2 = [0.0011,8.5408e-05,7.4385e-06,6.5977e-07];
ERRL2 = arrayfun(@(x) log(x),ERRL2);

% Calcul de l erreur H1
ERRH1 = (dot(KK*(UU_exact-UU),(UU_exact-UU)))^0.5/dot(UU_exact,UU_exact);
    
%On réalise le calcul pour différents pas h décroissants tels que
%h(n+1) = h(n)/2, les valeurs sont stockées dans la liste suivante
ERRH1 = [0.0048,6.7499e-04,8.7595e-05,1.0984e-05];
ERRH1 = arrayfun(@(x) log(x),ERRH1);
    
figure
    
p = polyfit(HH,ERRL2,1);
x1 = linspace(log(1/0.025),log(1/0.2));
y1 = polyval(p,x1);

p = polyfit(HH,ERRH1,1);
x2 = linspace(log(1/0.025),log(1/0.2));
y2 = polyval(p,x2);
plot(HH,ERRL2,'b-o',x1,y1,'b--',HH,ERRH1,'r-o',x2,y2,'r--');
hold
    
title('Graphe de convergence L2 et H1');
xlabel('log(1/h)');
ylabel('log(e)');
legend('erreur norme L2','interpolation lin L2 coeff = -3.56','erreur norme H1','interpolation lin H1 coeff = -2.93');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022



