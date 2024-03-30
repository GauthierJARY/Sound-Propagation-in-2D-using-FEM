function Kel = matK_elem(S1, S2, S3, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matK_elem :
% calcul la matrices de raideur elementaire en P1 lagrange.
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%       * Reftri : reference du triangle.
%
% OUTPUT * Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul utilise une formule de quadrature de Gauss-Legendre.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Reftri ~= 1 && Reftri ~= 2
    error('Un des triangles a une réference differente de 1 et 2.');
end

% preliminaires, pour faciliter la lecture.
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Points et poids de quadrature.
S_hat(:,1) = [1/6; 1/6];
S_hat(:,2) = [2/3; 1/6];
S_hat(:,3) = [1/6, 2/3];
w0 = 1/6;

% Gradients des fonctions de base sur le triangle de reference.
grad_w_hat(:,1) = [-1 ; -1];
grad_w_hat(:,2) = [1 ; 0];
grad_w_hat(:,3) = [0 ; 1];

% Transformation géométrique associée au triangle courant.
B_l = [x2-x1 x3-x1; y2-y1 y3-y1];
S_l = [x1 ; y1];

% Boucle sur les fonctions de bases locales.
Kel = zeros(3,3);
det_Bl = abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
inv_t_Bl = inv(B_l');
%Parcours des lignes de Kel
for i=1:3
    %Parcours des colonnes de Kel
    for j=1:3
        %Formule de quadrature à 3 points
        somme = 0;
        for k=1:3
            if Reftri == 1
                somme = somme + sigma_1(B_l*S_hat(:,k)+S_l)*(inv_t_Bl*grad_w_hat(:,i))'*inv_t_Bl*grad_w_hat(:,j);
            else
                QQ=B_l*S_hat(:,k)+S_l;
                somme = somme + sigma_2(QQ(1,1),QQ(2,1))*(inv_t_Bl*grad_w_hat(:,i))'*inv_t_Bl*grad_w_hat(:,j);
            end
        end
	Kel(i,j) = det_Bl*w0*somme;
    end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022


