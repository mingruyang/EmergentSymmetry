function sq = structurefactor2site_ncon(O1, O2, AL, AR, AC, L, R, d, D, p)
    ALc = conj(AL);
    ARc = conj(AR);
    ACc = conj(AC);
    O1e = ncon({O1,AL,AC,ALc,ACc},{[3 5 2 7],[1 2 8],[8 7 6],[1 3 4],[4 5 6]});
    O2e = ncon({O2,AL,AC,ALc,ACc},{[3 5 2 7],[1 2 8],[8 7 6],[1 3 4],[4 5 6]});
    O1n = O1 - O1e*reshape(kron(eye(d),eye(d)),[d,d,d,d]);
    O2n = O2 - O2e*reshape(kron(eye(d),eye(d)),[d,d,d,d]);
    sq1 = ncon({O1n,O2n,AL,AC,ALc,ACc},{[5 6 2 4],[10 8 5 6],[1 2 3],[3 4 7],[1 10 9],[9 8 7]});
    sq1 = sq1 + exp(-1i*p)*ncon({O1n,O2n,AL,AC,AR,ALc,ACc,ARc},{[13 7 2 4],[11 9 7 6],[1 2 3],[3 4 5],[5 6 8],[1 13 12],[12 11 10],[10 9 8]});
    sq1 = sq1 + exp(1i*p)*ncon({O1n,O2n,AL,AC,AR,ALc,ACc,ARc},{[7 13 4 2],[9 11 6 7],[8 6 5],[5 4 3],[3 2 1],[8 9 10],[10 11 12],[12 13 1]});

    O1L = ncon({O1n,AL,AL,ALc,ALc},{[7 9 4 6],[3 4 5],[5 6 -2],[3 7 8],[8 9 -1]},[3 4 5 6 7 8 9],[-1 -2]);
    O2R = ncon({O2n,AL,AC,ALc,ACc},{[9 7 6 4],[-1 6 5],[5 4 3],[-2 9 8],[8 7 3]},[3 4 5 6 7 8 9],[-1 -2]);
    O2L = ncon({O2n,AC,AR,ACc,ARc},{[7 9 4 6],[3 4 5],[5 6 -2],[3 7 8],[8 9 -1]},[3 4 5 6 7 8 9],[-1 -2]);
    O1R = ncon({O1n,AR,AR,ARc,ARc},{[9 7 6 4],[-1 6 5],[5 4 3],[-2 9 8],[8 7 3]},[3 4 5 6 7 8 9],[-1 -2]);
    afun = @(x) ELLregcomp(x, AL, ALc, R, D, p);
    sq2 = exp(-2*1i*p)*trace(reshape(gmres(afun,reshape(O1L,[D^2,1]),[],1E-14,1000),[D,D])*O2R);
    afun = @(x) ERRregcomp(x, AR, ARc, L, D, p);
    sq3 = exp(2*1i*p)*trace(O2L*reshape(gmres(afun,reshape(O1R,[D^2,1]),[],1E-14,1000),[D,D]));
    sq = sq1+sq2+sq3;
end

function y = ELLregcomp(x, AL, ALc, R, D, p)
    FL = reshape(x,[D,D]);
    yt = FL;
    yt = yt - exp(-1i*p)*ncon({FL,AL,ALc},{[3 4],[4 5 -2],[3 5 -1]},[3 4 5],[-1 -2]);
    yt = yt + (exp(-1i*p)*trace(FL*R))*eye(D);
    y = reshape(yt,[D^2,1]);
end

function y = ERRregcomp(x, AR, ARc, L, D, p)
    FR = reshape(x,[D,D]);
    yt = FR;
    yt = yt - exp(1i*p)*ncon({FR,AR,ARc},{[3 4],[-1 5 3],[-2 5 4]},[3 4 5],[-1 -2]);
    yt = yt + (exp(1i*p)*trace(L*FR))*eye(D);
    y = reshape(yt,[D^2,1]);
end