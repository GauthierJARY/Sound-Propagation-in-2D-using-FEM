function afficheSigma(Numtri, Reftri, Coorneu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% affiche:
% permet d'afficher la fonction Sigma (en tant que fonction par morceaux) 
% sur le maillage (Numtri, Coorneu)
%
% SYNOPSIS : afficheSigma(Numtri, Reftri, Coorneu)
%          
% INPUT * Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%       * Reftri : liste des références des triangles.
%       * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%
% OUTPUT une fenetre graphique
%
% NOTE 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy = max(Coorneu(:,2)) - min(Coorneu(:,2));
dx = max(Coorneu(:,1)) - min(Coorneu(:,1));
ratio = dx / dy;

figure; 
fig = gcf;
fig.Position(3) = ratio * fig.Position(4);
fig.Position(3:4) = 0.5 * fig.Position(3:4);

% Fonction sigma et coordonnées des noeuds en discontinu.
[Nbtri, Nloc] = size(Numtri);
Sigma = zeros(Nbtri * Nloc, 1);
CoorneuDisc = zeros(Nbtri * Nloc, 2);
NumtriDisc = zeros(Nbtri, 3);
for i=1:Nbtri
    
    if Reftri(i) ~= 1 && Reftri(i) ~= 2
        error('Un des triangles a une réference differente de 1 et 2.');
    end
    
    for j=1:Nloc
        disc_idx = (i - 1) * Nloc + j;
        xy = Coorneu(Numtri(i, j), :);
        if Reftri(i) == 1
            Sigma(disc_idx) = sigma_1(xy(1), xy(2));
        else
            Sigma(disc_idx) = sigma_2(xy(1), xy(2));
        end
        CoorneuDisc(disc_idx, :) = xy;
        NumtriDisc(i, j) = disc_idx;
    end
end

% control on the input args
if (nargin<4), titre = ''; end;
trisurf(NumtriDisc, CoorneuDisc(:,1), CoorneuDisc(:,2), Sigma);
pbaspect([ratio 1 1])
view(2);
shading interp
% shading faceted
% shading flat
colorbar;

% ajouter eventuellement un titre
title("Sigma");
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022



