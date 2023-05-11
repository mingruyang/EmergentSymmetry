function c = correlator2site_ncon(O1, O2, AL, AC, AR, d, r)
    ALc = conj(AL);
    ACc = conj(AC);
    ARc = conj(AR);
    O1e = ncon({O1,AL,AC,ALc,ACc},{[3 5 2 7],[1 2 8],[8 7 6],[1 3 4],[4 5 6]});
    O2e = ncon({O2,AL,AC,ALc,ACc},{[3 5 2 7],[1 2 8],[8 7 6],[1 3 4],[4 5 6]});
    O1n = O1 - O1e*reshape(kron(eye(d),eye(d)),[d,d,d,d]);
    O2n = O2 - O2e*reshape(kron(eye(d),eye(d)),[d,d,d,d]);

    c = zeros(r,1);
    c(1) = ncon({O1n,O2n,AL,AC,ALc,ACc},{[5 6 2 4],[10 8 5 6],[1 2 3],[3 4 7],[1 10 9],[9 8 7]});
    c(2) = ncon({O1n,O2n,AL,AC,AR,ALc,ACc,ARc},{[13 7 2 4],[11 9 7 6],[1 2 3],[3 4 5],[5 6 8],[1 13 12],[12 11 10],[10 9 8]});
    O1L = ncon({O1n,AL,AL,ALc,ALc},{[7 9 4 6],[3 4 5],[5 6 -2],[3 7 8],[8 9 -1]},[3 4 5 6 7 8 9],[-1 -2]);
    O2R = ncon({O2n,AL,AC,ALc,ACc},{[9 7 6 4],[-1 6 5],[5 4 3],[-2 9 8],[8 7 3]},[3 4 5 6 7 8 9],[-1 -2]);
    for i = 3:r
        c(i) = trace(O1L*O2R);
        O1L = leftApplyT(O1L, AL, ALc);
    end
end

function y = leftApplyT(x, AL, ALc)
    y = ncon({x,AL,ALc},{[3 4],[4 5 -2],[3 5 -1]},[3 4 5],[-1 -2]);
end