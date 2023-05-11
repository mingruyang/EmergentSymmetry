function sq = structurefactor4site_ncon(O1, O2, AL, AR, AC, L, R, d, D, p)
    ALc = conj(AL);
    ARc = conj(AR);
    ACc = conj(AC);
    O1e = ncon({O1,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[4 8 12 16 2 6 10 14],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[1 4 5],[5 8 9],[9 12 13],[13 16 15]});
    O2e = ncon({O2,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[4 8 12 16 2 6 10 14],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[1 4 5],[5 8 9],[9 12 13],[13 16 15]});
    O1n = O1 - O1e*reshape(kron(kron(kron(eye(d),eye(d)),eye(d)),eye(d)),[d,d,d,d,d,d,d,d]);
    O2n = O2 - O2e*reshape(kron(kron(kron(eye(d),eye(d)),eye(d)),eye(d)),[d,d,d,d,d,d,d,d]);
    sq1 = ncon({O1n,O2n,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[4 9 14 18 2 7 12 17],[5 10 15 19 4 9 14 18],[1 2 3],[3 7 8],[8 12 13],[13 17 20],[1 5 6],[6 10 11],[11 15 16],[16 19 20]});
    sq1 = sq1 + exp(-1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ACc},{[4 8 13 18 2 6 11 16],[9 14 19 22 8 13 18 21],[1 2 3],[3 6 7],[7 11 12],[12 16 17],[17 21 23],[1 4 5],[5 9 10],[10 14 15],[15 19 20],[20 22 23]});
    sq1 = sq1 + exp(-2*1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ALc,ACc},{[4 8 12 17 2 6 10 15],[13 18 22 25 12 17 20 24],[1 2 3],[3 6 7],[7 10 11],[11 15 16],[16 20 21],[21 24 26],[1 4 5],[5 8 9],[9 13 14],[14 18 19],[19 22 23],[23 25 26]});
    sq1 = sq1 + exp(-3*1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ALc,ALc,ACc},{[4 8 12 16 2 6 10 14],[17 21 25 28 16 19 23 27],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[15 19 20],[20 23 24],[24 27 29],[1 4 5],[5 8 9],[9 12 13],[13 17 18],[18 21 22],[22 25 26],[26 28 29]});
    sq1 = sq1 + exp(1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ACc},{[8 13 18 22 6 11 16 21],[4 9 14 19 2 8 13 18],[1 2 3],[3 6 7],[7 11 12],[12 16 17],[17 21 23],[1 4 5],[5 9 10],[10 14 15],[15 19 20],[20 22 23]});
    sq1 = sq1 + exp(2*1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ALc,ACc},{[12 17 22 25 10 15 20 24],[4 8 13 18 2 6 12 17],[1 2 3],[3 6 7],[7 10 11],[11 15 16],[16 20 21],[21 24 26],[1 4 5],[5 8 9],[9 13 14],[14 18 19],[19 22 23],[23 25 26]});
    sq1 = sq1 + exp(3*1i*p)*ncon({O1n,O2n,AL,AL,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ALc,ALc,ACc},{[16 21 25 28 14 19 23 27],[4 8 12 17 2 6 10 16],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[15 19 20],[20 23 24],[24 27 29],[1 4 5],[5 8 9],[9 12 13],[13 17 18],[18 21 22],[22 25 26],[26 28 29]});
    
    O1L = ncon({O1n,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[6 10 14 17 4 8 12 16],[3 4 5],[5 8 9],[9 12 13],[13 16 -2],[3 6 7],[7 10 11],[11 14 15],[15 17 -1]},[3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],[-1 -2]);
    O2R = ncon({O2n,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[17 14 10 6 16 12 8 4],[-1 16 13],[13 12 9],[9 8 5],[5 4 3],[-2 17 15],[15 14 11],[11 10 7],[7 6 3]},[3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],[-1 -2]);
    O2L = ncon({O2n,AC,AR,AR,AR,ACc,ARc,ARc,ARc},{[6 10 14 17 4 8 12 16],[3 4 5],[5 8 9],[9 12 13],[13 16 -2],[3 6 7],[7 10 11],[11 14 15],[15 17 -1]},[3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],[-1 -2]);
    O1R = ncon({O1n,AR,AR,AR,AR,ARc,ARc,ARc,ARc},{[17 14 10 6 16 12 8 4],[-1 16 13],[13 12 9],[9 8 5],[5 4 3],[-2 17 15],[15 14 11],[11 10 7],[7 6 3]},[3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],[-1 -2]);
    afun = @(x) ELLregcomp(x, AL, ALc, R, D, p);
    sq2 = exp(-4*1i*p)*trace(reshape(gmres(afun,reshape(O1L,[D^2,1]),[],1E-14,1000),[D,D])*O2R);
    afun = @(x) ERRregcomp(x, AR, ARc, L, D, p);
    sq3 = exp(4*1i*p)*trace(O2L*reshape(gmres(afun,reshape(O1R,[D^2,1]),[],1E-14,1000),[D,D]));
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