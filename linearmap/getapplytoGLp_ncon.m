function y = getapplytoGLp_ncon(AL, G, R, N, D, p, g)
    GAL = getGL_ncon(AL, G, N, D);
    ALc = conj(AL);
    GAL = GAL - g*eye(D);
    afun = @(x) ELLregcomp(x, AL, ALc, R, D, p);
    EGL = reshape(gmres(afun,reshape(GAL,[D^2,1]),[],1E-14,1000),[D,D]);
    if N == 1
        y = ncon({EGL,R,AL,ALc},{[3 4],[5 6],[4 -2 5],[3 -1 6]},[3 4 5 6],[-1 -2]);
        y = exp(-1i*p)*y;
    elseif N == 2
        y = ncon({EGL,R,AL,AL,ALc,ALc},{[5 6],[9 10],[6 -3 7],[7 -4 9],[5 -1 8],[8 -2 10]},[5 6 7 8 9 10],[-1 -2 -3 -4]);
        y = exp(-2*1i*p)*y;
    elseif N == 3
        y = ncon({EGL,R,AL,AL,AL,ALc,ALc,ALc},{[7 8],[13 14],[8 -4 9],[9 -5 12],[12 -6 13],[7 -1 10],[10 -2 11],[11 -3 14]},[7 8 9 10 11 12 13 14],[-1 -2 -3 -4 -5 -6]);
        y = exp(-3*1i*p)*y;
    elseif N == 4
        y = ncon({EGL,R,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[17 18],[9 10],[18 -5 15],[15 -6 13],[13 -7 11],[11 -8 9],[17 -1 16],[16 -2 14],[14 -3 12],[12 -4 10]});
        y = exp(-4*1i*p)*y;
    end
end

function y = ELLregcomp(x, AL, ALc, R, D, p)
    FL = reshape(x,[D,D]);
    yt = FL;
    yt = yt - exp(-1i*p)*ncon({FL,AL,ALc},{[3 4],[4 5 -2],[3 5 -1]},[3 4 5],[-1 -2]);
    yt = yt + (exp(-1i*p)*trace(FL*R))*eye(D);
    y = reshape(yt,[D^2,1]);
end
