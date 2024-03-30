function [Interp] = interpTriP1(Coorneu, Numtri, Coorpnts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpTriP1 :
% Calcul les coefficients d'interpolation d'un solution P1 en des points
% d'observation
%
% SYNOPSIS [Interp] = interpTriP1(Coorneu, Numtri, Coorpnts)
%          
% INPUT  * Coorneu : coordonnees des sommets (matrice réelle Nbpt x 2)
%        * Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        * Coorpnts : coordonnées des points d'observation (matrice réelle Nbpt obs x 2)
%
% OUTPUT * Interp : coefficient d'inerpolation sous la forme d'une matrice
% rectangulaire sparse à appliquer à une solution définie sur le maillage.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps = 1e-6;
Nbtri = size(Numtri, 1);
Nbneu = size(Coorneu, 1);
Nbpt = size(Coorpnts, 1);

II = zeros(Nbpt*3,1);
JJ = zeros(Nbpt*3,1);
VV_I = zeros(Nbpt*3,1);

for i = 1:Nbpt
    
    x = Coorpnts(i, 1);
    y = Coorpnts(i, 2);
    
    % Init. loop variables.
    triangleIndex = 0;
    interpolationCoefs = zeros(3, 1);
    j = 1;
    
    % Loop on mesh triangles.
    while triangleIndex == 0 && j <= Nbtri
        
        % Extracting mesh node coordinates.
        S1 = Coorneu(Numtri(j,1), :);
        S2 = Coorneu(Numtri(j,2), :);
        S3 = Coorneu(Numtri(j,3), :);
        
        % Inflating element for robustness.
        G = (S1 + S2 + S2) / 3.0;
        S1 = S1 + eps * (S1 - G);
        S2 = S2 + eps * (S2 - G);
        S3 = S3 + eps * (S3 - G);
        
        % Computing parametric coordinate in current element.
        A = [(S2 - S1)', (S3 - S1)'];
        hatS = A \ ([x; y] - S1');
        hatx = hatS(1);
        haty = hatS(2);
        
        % Testing if pnts is inside current element.
        if hatx >= -eps && haty >= -eps && 1.0 - hatx - haty >= -eps
            triangleIndex = j;
            interpolationCoefs(1) = 1.0 - hatx - haty;
            interpolationCoefs(2) = hatx;
            interpolationCoefs(3) = haty;
        end
        
        j = j + 1;
    end
    
    % Setting interpolator values if found.
    if triangleIndex == 0
        warning('Cannot find interpolation coefficient for point #' ...
            + string(i) + '.');
    else
        index = (3*(i-1)+1):(3*(i-1)+3);
        II(index) = i;
        JJ(index) = Numtri(triangleIndex, :);
        VV_I(index) = interpolationCoefs;
    end
    
end

if all(II == 0)
    Interp = sparse(Nbpt, Nbneu);
else
    Interp = sparse(II, JJ, VV_I, Nbpt, Nbneu);
end

end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022


