function c = correlator3site_ncon(O1, O2, AL, AC, AR, R, d, r)
    ALc = conj(AL);
    ACc = conj(AC);
    ARc = conj(AR);
    O1e = ncon({O1,AL,AC,AR,ALc,ACc,ARc},{[3 7 11 2 6 10],[1 2 5],[5 6 9],[9 10 12],[1 3 4],[4 7 8],[8 11 12]});
    O2e = ncon({O2,AL,AC,AR,ALc,ACc,ARc},{[3 7 11 2 6 10],[1 2 5],[5 6 9],[9 10 12],[1 3 4],[4 7 8],[8 11 12]});
    O1n = O1 - O1e*reshape(kron(kron(eye(d),eye(d)),eye(d)),[d,d,d,d,d,d]);
    O2n = O2 - O2e*reshape(kron(kron(eye(d),eye(d)),eye(d)),[d,d,d,d,d,d]);

    c = zeros(r,1);
    c(1) = ncon({O1n,O2n,AL,AC,AR,ALc,ACc,ARc},{[9 8 7 2 4 6],[10 12 14 9 8 7],[1 2 3],[3 4 5],[5 6 15],[1 10 11],[11 12 13],[13 14 15]});
    c(2) = ncon({O1n,O2n,AL,AL,AL,AL,ALc,ALc,ALc,ALc,R},{[19 10 9 2 4 6],[17 15 13 10 9 8],[1 2 3],[3 4 5],[5 6 7],[7 8 11],[1 19 18],[18 17 16],[16 15 14],[14 13 12],[11 12]});
    c(3) = ncon({O1n,O2n,AL,AL,AL,AL,AL,ALc,ALc,ALc,ALc,ALc,R},{[22 20 11 2 4 6],[18 16 14 11 8 10],[1 2 3],[3 4 5],[5 6 7],[7 8 9],[9 10 12],[1 22 21],[21 20 19],[19 18 17],[17 16 15],[15 14 13],[12 13]});
    O1L = ncon({O1n,AL,AL,AL,ALc,ALc,ALc},{[7 5 3 9 11 13],[8 9 10],[10 11 12],[12 13 -2],[8 7 6],[6 5 4],[4 3 -1]},[3 4 5 6 7 8 9 10 11 12 13],[-1 -2]);
    O2R = ncon({O2n,AL,AL,AC,ALc,ALc,ACc},{[3 5 7 13 11 9],[-1 13 12],[12 11 10],[10 9 8],[-2 3 4],[4 5 6],[6 7 8]},[3 4 5 6 7 8 9 10 11 12 13],[-1 -2]);
    for i = 4:r
        c(i) = trace(O1L*O2R);
        O1L = leftApplyT(O1L, AL, ALc);
    end
end

function y = leftApplyT(x, AL, ALc)
    y = ncon({x,AL,ALc},{[3 4],[4 5 -2],[3 5 -1]},[3 4 5],[-1 -2]);
end