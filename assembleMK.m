function [MME, KKE] = assembleMK(Coorneu, Refneu, Numtri, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembleMK :
% assemble les matrices de masse et de raideur globales en P1 lagrange.
%
% SYNOPSIS [M, K] = assembleMK(Coorneu, Numtri, Reftri)
%
% INPUT  * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        * Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        * Numtri : liste de triangles
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        * Reftri : reference des triangles (matrice entiere Nbtri x 1)
%
% OUTPUT * M matrice de masse globale (matrice NbptxNbpt)
%        * K matrice de raideur globale (matrice NbptxNbpt)
%
% NOTE (1) les matrices de masse et de raideur sont définies au format
% sparse.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbpt = size(Coorneu, 1);
Nbtri = size(Numtri, 1);

% Declarations des matrices EF.
KK = sparse(Nbpt,Nbpt);
MM = sparse(Nbpt,Nbpt);

% Boucle d'assemblage sur les triangles.
for l=1:Nbtri
    
    S1 = Coorneu(Numtri(l,1), :);
    S2 = Coorneu(Numtri(l,2), :);
    S3 = Coorneu(Numtri(l,3), :);
    
    % Numéro Global
    i = Numtri(l,1);
    j = Numtri(l,2);
    k = Numtri(l,3);
    
    % Calcul des matrices elementaires du triangle l.
    Kel=matK_elem(S1, S2, S3, Reftri(l));
    Mel=matM_elem(S1, S2, S3);
    
    % Assemblage des matrices globales.
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

MME = MM;
KKE = KK;
% Application de la pseudo élimination.
for l=1:Nbpt
    if Refneu(l) == 1 || Refneu(l) == 2
        KKE(:,l) = 0;
        KKE(l,:) = 0;
        KKE(l,l) = KK(l,l);
    end
    if Refneu(l) == 1 || Refneu(l) == 2
        MME(l,1) = 0;
        MME(:,l) = 0;
        MME(l,:) = 0;
        MME(l,l) = MM(l,l);
    end
end

MME = sparse(MME);
KKE = sparse(KKE);

end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022

