function c = correlator4site_ncon(O1, O2, AL, AC, d, r)
    ALc = conj(AL);
    ACc = conj(AC);
    O1e = ncon({O1,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[4 8 12 16 2 6 10 14],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[1 4 5],[5 8 9],[9 12 13],[13 16 15]});
    O2e = ncon({O2,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[4 8 12 16 2 6 10 14],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[1 4 5],[5 8 9],[9 12 13],[13 16 15]});
    O1n = O1 - O1e*reshape(kron(kron(kron(eye(d),eye(d)),eye(d)),eye(d)),[d,d,d,d,d,d,d,d]);
    O2n = O2 - O2e*reshape(kron(kron(kron(eye(d),eye(d)),eye(d)),eye(d)),[d,d,d,d,d,d,d,d]);

    c = zeros(r,1);
    c(1) = ncon({O1n,O2n,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[4 9 14 18 2 7 12 17],[5 10 15 19 4 9 14 18],[1 2 3],[3 7 8],[8 12 13],[13 17 20],[1 5 6],[6 10 11],[11 15 16],[16 19 20]});
    c(2) = ncon({O1n,O2n,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ACc},{[4 8 13 18 2 6 11 16],[9 14 19 22 8 13 18 21],[1 2 3],[3 6 7],[7 11 12],[12 16 17],[17 21 23],[1 4 5],[5 9 10],[10 14 15],[15 19 20],[20 22 23]});
    c(3) = ncon({O1n,O2n,AL,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ALc,ACc},{[4 8 12 17 2 6 10 15],[13 18 22 25 12 17 20 24],[1 2 3],[3 6 7],[7 10 11],[11 15 16],[16 20 21],[21 24 26],[1 4 5],[5 8 9],[9 13 14],[14 18 19],[19 22 23],[23 25 26]});
    c(4) = ncon({O1n,O2n,AL,AL,AL,AL,AL,AL,AC,ALc,ALc,ALc,ALc,ALc,ALc,ACc},{[4 8 12 16 2 6 10 14],[17 21 25 28 16 19 23 27],[1 2 3],[3 6 7],[7 10 11],[11 14 15],[15 19 20],[20 23 24],[24 27 29],[1 4 5],[5 8 9],[9 12 13],[13 17 18],[18 21 22],[22 25 26],[26 28 29]});
    O1L = ncon({O1n,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[6 10 14 17 4 8 12 16],[3 4 5],[5 8 9],[9 12 13],[13 16 -2],[3 6 7],[7 10 11],[11 14 15],[15 17 -1]},[3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],[-1 -2]);
    O2R = ncon({O2n,AL,AL,AL,AC,ALc,ALc,ALc,ACc},{[17 14 10 6 16 12 8 4],[-1 16 13],[13 12 9],[9 8 5],[5 4 3],[-2 17 15],[15 14 11],[11 10 7],[7 6 3]},[3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],[-1 -2]);
    for i = 5:r
        c(i) = trace(O1L*O2R);
        O1L = leftApplyT(O1L, AL, ALc);
    end
end

function y = leftApplyT(x, AL, ALc)
    y = ncon({x,AL,ALc},{[3 4],[4 5 -2],[3 5 -1]},[3 4 5],[-1 -2]);
end