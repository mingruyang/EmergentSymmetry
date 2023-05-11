function y = getapplytoGRp_ncon(AR, G, L, N, D, p, g)
    GAR = getGR_ncon(AR, G, N, D);
    ARc = conj(AR);
    GAR = GAR - g*eye(D);
    afun = @(x) ERRregcomp(x, AR, ARc, L, D, p);
    EGR = reshape(gmres(afun,reshape(GAR,[D^2,1]),[],1E-14,1000),[D,D]);
    if N == 1
        y = ncon({EGR,L,AR,ARc},{[5 6],[3 4],[4 -2 5],[3 -1 6]},[3 4 5 6],[-1 -2]);
        y = exp(1i*p)*y;
    elseif N == 2
        y = ncon({EGR,L,AR,AR,ARc,ARc},{[9 10],[5 6],[6 -3 7],[7 -4 9],[5 -1 8],[8 -2 10]},[5 6 7 8 9 10],[-1 -2 -3 -4]);
        y = exp(2*1i*p)*y;
    elseif N == 3
        y = ncon({EGR,L,AR,AR,AR,ARc,ARc,ARc},{[13 14],[7 8],[8 -4 9],[9 -5 12],[12 -6 13],[7 -1 10],[10 -2 11],[11 -3 14]},[7 8 9 10 11 12 13 14],[-1 -2 -3 -4 -5 -6]);
        y = exp(3*1i*p)*y;
    elseif N == 4
        y = ncon({EGR,L,AR,AR,AR,AR,ARc,ARc,ARc,ARc},{[9 10],[17 18],[18 -5 15],[15 -6 13],[13 -7 11],[11 -8 9],[17 -1 16],[16 -2 14],[14 -3 12],[12 -4 10]});
        y = exp(4*1i*p)*y;
    end
end

function y = ERRregcomp(x, AR, ARc, L, D, p)
    FR = reshape(x,[D,D]);
    yt = FR;
    yt = yt - exp(1i*p)*ncon({FR,AR,ARc},{[3 4],[-1 5 3],[-2 5 4]},[3 4 5],[-1 -2]);
    yt = yt + (exp(1i*p)*trace(L*FR))*eye(D);
    y = reshape(yt,[D^2,1]);
end
