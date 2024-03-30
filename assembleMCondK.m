function [MCond, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembleMCondK :
% assemble les matrices de masse condensée et de raideur globales en P1 lagrange.
%
% SYNOPSIS [MCond, K] = assembleMCondK(Coorneu, Numtri, Reftri)
%          
% INPUT  * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        * Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        * Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        * Reftri : reference des triangles (matrice entiere Nbtri x 1)
%
% OUTPUT * M matrice de masse globale (vecteur Nbpt)
%        * K matrice de raideur globale (matrice NbptxNbpt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbpt = size(Coorneu, 1);
Nbtri = size(Numtri, 1);

% Declarations des matrices EF.
KK = sparse(Nbpt,Nbpt);
MCondDiag = zeros(Nbpt, 1);

% Boucle d'assemblage sur les triangles.
for l=1:Nbtri
    
    S1 = Coorneu(Numtri(l,1), :);
    S2 = Coorneu(Numtri(l,2), :);
    S3 = Coorneu(Numtri(l,3), :);
    
    % Assemblage de la matrice de rigidité.
    % A COMPLETER
    % Numéro Global
    i = Numtri(l,1);
    j = Numtri(l,2);
    k = Numtri(l,3);
    Kel=matK_elem(S1,S2,S3,Reftri(l));
    KK(i,i) = KK(i,i) + Kel(1,1);
    KK(i,j) = KK(i,j) + Kel(1,2);
    KK(i,k) = KK(i,k) + Kel(1,3);
    KK(j,j) = KK(j,j) + Kel(2,2);
    KK(k,k) = KK(k,k) + Kel(3,3);
    KK(j,k) = KK(j,k) + Kel(2,3);
    KK(k,j) = KK(k,j) + Kel(3,2);
    KK(k,i) = KK(k,i) + Kel(3,1);
    KK(j,i) = KK(j,i) + Kel(2,1);
    
    % Assemblage de la diagonale de la matrice de masse.
    x1 = S1(1); y1 = S1(2);
    x2 = S2(1); y2 = S2(2);
    x3 = S3(1); y3 = S3(2);
    det_Bl = abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    MCondDiag(i,1)=MCondDiag(i,1) + det_Bl/6;
    MCondDiag(j,1)=MCondDiag(j,1) + det_Bl/6;
    MCondDiag(k,1)=MCondDiag(k,1) + det_Bl/6;
    
end % for l

% Application de la pseudo élimination.
K=KK;
for l=1:Nbpt
    % A COMPLETER
    if Refneu(l) == 1 || Refneu(l) == 2
        K(:,l) = 0;
        K(l,:) = 0;
        K(l,l) = KK(l,l);
    end
end

% Transformation de la diagonale en une matrice sparse.
MCond = spdiags(MCondDiag, 0, Nbpt, Nbpt);
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022

