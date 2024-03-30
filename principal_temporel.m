% =====================================================
% principal_temporel;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante des ondes en r?gime temporel, avec conditions de
% Dirichlet homogene
% | d^2_{tt} u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

%% Reading mesh file and assembling matrices.
clear all

meshFilePath = "geomRect005.msh";
% massType = 'Exacte';  
massType = 'Cond';

[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri]=lecture_msh(meshFilePath);

if strcmp(massType, 'Exacte')
    [M, K] = assembleMK(Coorneu, Refneu, Numtri, Reftri);
elseif strcmp(massType, 'Cond')
   [M, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri);
end

%% Computing CLF condition and effective time step.
% Calcul de la CFL
VP = eigs(K,M,1);
cfl = 2/sqrt(VP);
 
% Calcul du pas de temps.
cfl_factor = 0.95;
dt = cfl_factor * cfl;


%% Computing initial conditions.

interpU0 = zeros(Nbpt, 1);
interpU1 = zeros(Nbpt, 1);
for i = 1:Nbpt
    interpU0(i,1) = exp(-50*((Coorneu(i,1)-3)^2 + (Coorneu(i,2)-1)^2));
end


%% Propagating.

Tmax = 3;
niter = floor(Tmax / dt) + 1;
% solverType = 'Cholesky'; 
solverType = 'Cond';

if strcmp(solverType, 'Cholesky')
    [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter);
elseif strcmp(solverType, 'Cond')
    [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter);
end


%% Plots of energy and solution.

figure
hold on;
plot(Times, Potential)
plot(Times, Kinetic)
plot(Times, Kinetic + Potential)
xlim([min(Times) max(Times)])
xlabel('Time')
ylabel('Energy')
title('Energie du système au cours du temps, pas h = 0.05')
legend('P', 'K', 'E', 'Location', 'SouthEast')

affiche(Us(:, end), Numtri, Coorneu, "Solution approchée à milieu homogène, pas h = 0.05");
afficheSigma(Numtri, Reftri, Coorneu);


%% Interpolation at point.

CoordsInterpPnts = [2, 1; 4, 1];
NbInterpPnts = size(CoordsInterpPnts, 1);

% Calcul de la matrice des co?fficients d'interpolation.
interpolationOp = interpTriP1(Coorneu, Numtri, CoordsInterpPnts);

u_interpolee = interpolationOp*Us(:,end);

plot([0 Times'], u_interpolee(1,1:end-1), 'r', [0 Times'], u_interpolee(1,1:end-1), 'b');
title("Evolution de la solution aux point d'intérêt (2,1) et (4,1)");
legend("Point (2,1)", "Point (4,1)")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022

