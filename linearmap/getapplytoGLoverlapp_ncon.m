function y = getapplytoGLoverlapp_ncon(AL, G, R, N, d, p, g, yg)
    if N == 1
        y = zeros(d,d);
    elseif N == 2
        ALc = conj(AL);
        y = ncon({G,R,AL,AL,AL,ALc,ALc,ALc},{[2 6 3 -7],[13 14],[1 2 5],[5 6 9],[9 -11 13],[1 3 4],[4 -8 10],[10 -12 14]},[1 2 3 4 5 6 9 10 13 14],[-8 -12 -7 -11]);
        y = y - g*yg;
        y = exp(-1i*p)*y;
    elseif N == 3
        ALc = conj(AL);
        y = exp(-1i*p)*ncon({G,R,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[2 6 11 3 -7 -12],[18 19],[1 2 5],[5 6 9],[9 11 15],[15 -16 18],[1 3 4],[4 -8 10],[10 -13 14],[14 -17 19]},[1 2 3 4 5 6 9 10 11 14 15 18 19],[-8 -13 -17 -7 -12 -16]);
        y = y - (exp(-1i*p)*g)*yg;
        y = y + exp(-2*1i*p)*ncon({G,R,AL,AL,AL,AL,AL,ALc,ALc,ALc,ALc,ALc},{[2 6 11 3 7 -10],[21 22],[1 2 5],[5 6 12],[12 11 13],[13 -16 17],[17 -20 21],[1 3 4],[4 7 8],[8 -9 14],[14 -15 18],[18 -19 22]},[1 2 3 4 5 6 7 8 12 11 13 14 17 18 21 22],[-9 -15 -19 -10 -16 -20]);
        y = y - (exp(-2*1i*p)*g)*yg;
    elseif N == 4
        ALc = conj(AL);
        y = exp(-1i*p)*ncon({G,R,AL,AL,AL,AL,AL,ALc,ALc,ALc,ALc,ALc},{[2 6 11 16 3 -7 -12 -17],[23 24],[1 2 4],[4 6 10],[10 11 15],[15 16 20],[20 -21 23],[1 3 5],[5 -8 9],[9 -13 14],[14 -18 19],[19 -22 24]},[1 2 3 4 5 6 9 10 11 14 15 16 19 20 23 24],[-8 -13 -18 -22 -7 -12 -17 -21]);
        y = y - (exp(-1i*p)*g)*yg;
        y = y + exp(-2*1i*p)*ncon({G,R,AL,AL,AL,AL,AL,AL,ALc,ALc,ALc,ALc,ALc,ALc},{[2 6 10 15 3 7 -11 -16],[26 27],[1 2 4],[4 6 9],[9 10 14],[14 15 19],[19 -20 23],[23 -24 26],[1 3 5],[5 7 8],[8 -12 13],[13 -17 18],[18 -21 22],[22 -25 27]},[1 2 3 4 5 6 7 8 9 10 13 14 15 18 19 22 23 26 27],[-12 -17 -21 -25 -11 -16 -20 -24]);
        y = y - (exp(-2*1i*p)*g)*yg;
        y = y + exp(-3*1i*p)*ncon({G,R,AL,AL,AL,AL,AL,AL,AL,ALc,ALc,ALc,ALc,ALc,ALc,ALc},{[2 6 10 14 3 7 11 -15],[29 30],[1 2 4],[4 6 9],[9 10 13],[13 14 18],[18 -19 22],[22 -23 26],[26 -27 29],[1 3 5],[5 7 8],[8 11 12],[12 -16 17],[17 -20 21],[21 -24 25],[25 -28 30]},[1 2 3 4 5 6 7 8 9 10 11 12 13 14 17 18 21 22 25 26 29 30],[-16 -20 -24 -28 -15 -19 -23 -27]);
        y = y - (exp(-3*1i*p)*g)*yg;
    end    
end
