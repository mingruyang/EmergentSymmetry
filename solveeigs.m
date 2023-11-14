datadir = '/path to the MPS tensor for certain bond dimension/D';
writedir = '/write path for certain bond dimension/D';
d = 2; % dimension of 1-site physical Hilbert space
N = 2; % number of sites G is supported from, currently only N = 1,2,3,4 have been implemented
p = 0; % momentum
poddteven = true; % quantum number for the parity (p) and time reversal (t) symmetry
D = 10; % bond dimension of MPS

% read the MPS tensors which has been optimized by VUMPS
AC = load([datadir num2str(D) '/AC.mat']).AC;
AL = load([datadir num2str(D) '/AL.mat']).AL;
AR = load([datadir num2str(D) '/AR.mat']).AR;
C = load([datadir num2str(D) '/C.mat']).C;
R = C*C';
L = C'*C;
AR = {AR{1};AR{2}}; % re-ordering the cell matrices
ALt = reshape(cell2mat(AL),[D,d,D]);
ARt = reshape(cell2mat(AR),[D,d,D]);

% solve the eigenvalue problem for D = 10
AfunG = @(x) applytoGp_ncon(x, ALt, ARt, L, R, N, d, D, p, poddteven);
[VG, lambdaG] = eigs(AfunG,d^(2*N),d^(2*N));
eigm = diag(lambdaG);
format long;
disp(eigm);
save([writedir num2str(D) '/eigpoddteven_' num2str(N) 'site_0_variational_Gwhole.mat'],'eigm');
save([writedir num2str(D) '/Vpoddteven_' num2str(N) 'site_0_variational_Gwhole.mat'],'VG');

neigs = 4; % number of non-zero eigenvalues, read from the D = 10 result
% solve the eigenvalue problem other bond dimensions
for D = [20:10:100,120,160,240,400]
    AC = load([datadir num2str(D) '/AC.mat']).AC;
    AL = load([datadir num2str(D) '/AL.mat']).AL;
    AR = load([datadir num2str(D) '/AR.mat']).AR;
    C = load([datadir num2str(D) '/C.mat']).C;
    AR = {AR{1};AR{2}};
    ALt = reshape(cell2mat(AL),[D,d,D]);
    ARt = reshape(cell2mat(AR),[D,d,D]);
    R = C*C';
    L = C'*C;

    AfunG = @(x) applytoGp_ncon(x, ALt, ARt, L, R, N, d, D, p, poddteven);
    [VG, lambdaG] = eigs(AfunG,d^(2*N),neigs);
    eigm = diag(lambdaG);
    format long;
    disp(eigm);
    save([writedir num2str(D) '/eigpoddteven_' num2str(N) 'site_0_variational_Gwhole.mat'],'eigm');
    save([writedir num2str(D) '/Vpoddteven_' num2str(N) 'site_0_variational_Gwhole.mat'],'VG');
end

% show the coefficients of the Pauli strings in G
k = 3; % choose the 3rd largest eigenvalue
G = reshape(VG(:,k),[2^N,2^N]);
op = pauliproductbasis(N);
[g,place] = Gsignificant(G, N, 1e-16);
coef = Gtranslator(g, place, op, N);
abscoef = abs(coef);
stem(1:2^(2*N),abscoef);

for i = 1:length(abscoef)
    if(abscoef(i) > 1E-4)
        text(i,abscoef(i),dec2base(i-1,4,N));
        disp(['+(' num2str(coef(i),'%.15f') ')[' dec2base(i-1,4,N) ']']);
    end
end
xlabel('$n$','interpreter','latex','FontSize',24);
ylabel('$|c(n)|$','interpreter','latex','FontSize',24);
xlim([1 2^(2*N)]);
