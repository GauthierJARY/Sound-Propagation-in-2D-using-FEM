% =====================================================
% principal_CFL;
%
% une routine permettant d'évaluer la dépendance de la condition CFL
% en fonction du pas de maillage
%
% =====================================================

%% Definition of mesh files and associated mesh steps.

meshFilePath = ["geomRect02.msh", "geomRect015.msh", "geomRect01.msh", "geomRect0075.msh" "geomRect005.msh", "geomRect0025.msh"];

steps = [0.2, 0.15, 0.1, 0.075, 0.05, 0.025];

%% Computing CFL for every meshes.

nbMesh = length(meshFilePath);
if nbMesh ~= length(meshFilePath)
    error('Number of mesh files and number of steps differ.');
end

% Allocation.
cfl = zeros(nbMesh, 1);
cflcond = zeros(nbMesh, 1);
% Boucle sur les maillages.
for i=1:nbMesh
    
    [Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(meshFilePath(i));
    
    % Assemblage de M et K et calcul de la CFL.
    [M,K] = assembleMK(Coorneu, Refneu, Numtri, Reftri);
    [Mcond,Kcond]=assembleMCondK(Coorneu,Refneu, Numtri, Reftri);
    %Calcul du rayon spectral = Plus grande Valeur propre du problème
    %généralisé
    VP = eigs(K, M, 1);
    VPcond = eigs(Kcond, Mcond, 1); 
    
    cfl(i) = 2/sqrt(VP);
    cflcond(i) = 2/sqrt(VPcond);
end


%% Affichage de la cfl en fonction de h.

p = polyfit(steps,cfl,1);
pcond = polyfit(steps,cflcond,1);
X = linspace(steps(length(steps)),0.2);
Y = polyval(p,X);
Ycond = polyval(pcond,X);
plot(steps,cfl','b-o',X,Y,'b--',steps,cflcond','g-o',X,Ycond,'g--');
title('Condition de stabilité dite CFL fonction du pas de maillage h');
xlabel('h');
ylabel('CFL');
legend('CFL calculé','interpolation linéaire de coeff 0.213','CFL calculé avec Condensation','interpolation linéaire de coeff 0.266');

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022

