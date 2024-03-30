function [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage_cond :
% Propage les conditions initiales à partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter)
%
% INPUT  * M, K : matrice de masse condensée et de rigidité
%        * interpU0, interpU1 : interpolée en espace des conditions initiales
%        * dt : pas de temps du schéma.
%        * niter : nombre d'itérations
%
% OUTPUT * Us : saolution discrète à tous les pas de temps
%               avec les conditions initiales (matrice Nbpt x niter + 2)
%        * Kinetic, Potential : énergies cinétique et potentielle (vecteur niter)
%        * Times : vecteur temps (vecteur niter)
%
% NOTE (1) On suppose que la matrice de masse est diagonale de sorte que
% l'opération M \ b est efficace.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation.
Nbpts = length(interpU0);
Us = zeros(Nbpts, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);

% Conditions initiales.
Us(:, 1) = interpU0; %Correspond au vecteur V(0), on applique u0 à tous les points.
Us(:, 2) = (speye(Nbpts)-(dt^2)/2*inv(M)*K)*interpU0+dt*interpU1; %Correspond au vecteur V1, avec les conditions initiales calculées Q2.4

% Interiations.
for i = 1:niter
        
    % Calcul des ?nergies.
    Kinetic(i) = 0.5*(1/dt)^2*(Us(:,i+1)-Us(:,i))'*(M-(dt^2/4)*K)*(Us(:,i+1)-Us(:,i));
    Potential(i) = (Us(:,i+1)+Us(:,i))'*K*(Us(:,i+1)+Us(:,i))/8;
    
    % Calcul de la solution par résolution directe du système linéaire
    % (suppose que l'opération M \ b soit efficace, i.e. M diagonale).
    Us(:,i+2)= 2*Us(:,i+1) -Us(:,i)+ M \ (  -(dt^2)*K*Us(:,i+1)  );
    
    Times(i) = i * dt;

end
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022

