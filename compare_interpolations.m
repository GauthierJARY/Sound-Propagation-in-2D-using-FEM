%% Feuille de comparaison des solutions obtenues avec condensation de masse et calcul exact interpolées 

clear all

meshFilePath = "geomRect005.msh";

[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri]=lecture_msh(meshFilePath);

[M, K] = assembleMK(Coorneu, Refneu, Numtri, Reftri);
[MC, KC] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri);

VP = eigs(K,M,1);
VP2 = eigs(KC,MC,1);
cfl = 2/sqrt(VP);
cfl2 = 2/sqrt(VP2);

cfl_factor = 0.95;
dt = cfl_factor * cfl;
dt2 = cfl_factor * cfl2;

interpU0 = zeros(Nbpt, 1);
interpU1 = zeros(Nbpt, 1);
for i = 1:Nbpt
    interpU0(i,1) = exp(-50*((Coorneu(i,1)-3)^2 + (Coorneu(i,2)-1)^2));
end

Tmax = 3;
niter = floor(Tmax / dt) + 1;
niter2 = floor(Tmax / dt2) + 1;

[Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter);
[UsC, KineticC, PotentialC, TimesC] = propage_cond(MC, KC, interpU0, interpU1, dt2, niter2);

%%
CoordsInterpPnts = [2, 1; 4, 1];
NbInterpPnts = size(CoordsInterpPnts, 1);

% Calcul de la matrice des co?fficients d'interpolation.
interpolationOp = interpTriP1(Coorneu, Numtri, CoordsInterpPnts);

u_interpolee = interpolationOp*Us(:,:);
u_interpoleeC = interpolationOp*UsC(:,:);

plot([0 Times'], u_interpolee(1,1:end-1), 'r', [0 TimesC'], u_interpoleeC(1,1:end-1), 'b');
title("Evolution de la solution de matrice de exacte ou cendensée, au point d'intérêt (2,1)");
legend("Solution méthode exacte", "Solution méthode condensée")