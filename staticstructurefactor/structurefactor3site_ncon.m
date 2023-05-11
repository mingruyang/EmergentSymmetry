function sq = structurefactor3site_ncon(O1, O2, AL, AR, AC, L, R, d, D, p)
    ALc = conj(AL);
    ARc = conj(AR);
    ACc = conj(AC);
    O1e = ncon({O1,AL,AC,AR,ALc,ACc,ARc},{[3 7 11 2 6 10],[1 2 5],[5 6 9],[9 10 12],[1 3 4],[4 7 8],[8 11 12]});
    O2e = ncon({O2,AL,AC,AR,ALc,ACc,ARc},{[3 7 11 2 6 10],[1 2 5],[5 6 9],[9 10 12],[1 3 4],[4 7 8],[8 11 12]});
    O1n = O1 - O1e*reshape(kron(kron(eye(d),eye(d)),eye(d)),[d,d,d,d,d,d]);
    O2n = O2 - O2e*reshape(kron(kron(eye(d),eye(d)),eye(d)),[d,d,d,d,d,d]);
    sq1 = ncon({O1n,O2n,AL,AC,AR,ALc,ACc,ARc},{[9 8 7 2 4 6],[10 12 14 9 8 7],[1 2 3],[3 4 5],[5 6 15],[1 10 11],[11 12 13],[13 14 15]});
    sq1 = sq1 + exp(-1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,ALc,ALc,ALc,ALc,R},{[19 10 9 2 4 6],[17 15 13 10 9 8],[1 2 3],[3 4 5],[5 6 7],[7 8 11],[1 19 18],[18 17 16],[16 15 14],[14 13 12],[11 12]});
    sq1 = sq1 + exp(-2*1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AL,ALc,ALc,ALc,ALc,ALc,R},{[22 20 11 2 4 6],[18 16 14 11 8 10],[1 2 3],[3 4 5],[5 6 7],[7 8 9],[9 10 12],[1 22 21],[21 20 19],[19 18 17],[17 16 15],[15 14 13],[12 13]});
    sq1 = sq1 + exp(1i*p)*ncon({O1n,O2n,AR,AR,AR,AR,ARc,ARc,ARc,ARc,L},{[9 10 19 6 4 2],[13 15 17 8 9 10],[3 2 1],[5 4 3],[7 6 5],[11 8 7],[18 19 1],[16 17 18],[14 15 16],[12 13 14],[12 11]});
    sq1 = sq1 + exp(2*1i*p)*ncon({O1n,O2n,AR,AR,AR,AR,AR,ARc,ARc,ARc,ARc,ARc,L},{[11 20 22 6 4 2],[14 16 18 10 8 11],[3 2 1],[5 4 3],[7 6 5],[9 8 7],[12 10 9],[21 22 1],[19 20 21],[17 18 19],[15 16 17],[13 14 15],[13 12]});

    O1L = ncon({O1n,AL,AL,AL,ALc,ALc,ALc},{[7 5 3 9 11 13],[8 9 10],[10 11 12],[12 13 -2],[8 7 6],[6 5 4],[4 3 -1]},[3 4 5 6 7 8 9 10 11 12 13],[-1 -2]);
    O2R = ncon({O2n,AL,AL,AC,ALc,ALc,ACc},{[3 5 7 13 11 9],[-1 13 12],[12 11 10],[10 9 8],[-2 3 4],[4 5 6],[6 7 8]},[3 4 5 6 7 8 9 10 11 12 13],[-1 -2]);
    O2L = ncon({O2n,AC,AR,AR,ACc,ARc,ARc},{[7 5 3 9 11 13],[8 9 10],[10 11 12],[12 13 -2],[8 7 6],[6 5 4],[4 3 -1]},[3 4 5 6 7 8 9 10 11 12 13],[-1 -2]);
    O1R = ncon({O1n,AR,AR,AR,ARc,ARc,ARc},{[3 5 7 13 11 9],[-1 13 12],[12 11 10],[10 9 8],[-2 3 4],[4 5 6],[6 7 8]},[3 4 5 6 7 8 9 10 11 12 13],[-1 -2]);
    afun = @(x) ELLregcomp(x, AL, ALc, R, D, p);
    sq2 = exp(-3*1i*p)*trace(reshape(gmres(afun,reshape(O1L,[D^2,1]),[],1E-14,1000),[D,D])*O2R);
    afun = @(x) ERRregcomp(x, AR, ARc, L, D, p);
    sq3 = exp(3*1i*p)*trace(O2L*reshape(gmres(afun,reshape(O1R,[D^2,1]),[],1E-14,1000),[D,D]));
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