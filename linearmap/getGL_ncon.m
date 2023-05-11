function GL = getGL_ncon(AL, G, N, D)
    GL = zeros(D);
    ALc = conj(AL);
    if N == 1
        GL = ncon({G,AL,ALc},{[2 3],[1 2 -5],[1 3 -4]},[1 2 3],[-4 -5]);
    elseif N == 2
        GL = ncon({G,AL,AL,ALc,ALc},{[3 7 2 6],[1 3 5],[5 7 -9],[1 2 4],[4 6 -8]},[1 2 3 4 5 6 7],[-8 -9]);
    elseif N == 3
        GL = ncon({G,AL,AL,AL,ALc,ALc,ALc},{[3 7 11 2 6 10],[1 3 5],[5 7 9],[9 11 -13],[1 2 4],[4 6 8],[8 10 -12]},[1 2 3 4 5 6 7 8 9 10 11],[-12 -13]);
    elseif N == 4
        GL = ncon({G,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[2 6 10 14 3 7 11 15],[1 2 4],[4 6 8],[8 10 12],[12 14 -17],[1 3 5],[5 7 9],[9 11 13],[13 15 -16]},[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15],[-16 -17]);
    end
end
