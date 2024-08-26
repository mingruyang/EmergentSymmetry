function y = applytoGp_ncon(x, AL, AR, L, R, N, d, D, p, podd, teven)
    sy = [0 -1i;1i 0];
    if podd == 1 && teven == 1
        if N == 1
            G = reshape(x,d*ones(1,2*N));
        elseif N == 2
            G = reshape(x,[d,d,d,d]);
            Gp = permute(G,[2 1 4 3]);
            if p == 0
                GG = reshape((G-Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^mod(N,2)*Gp)/2,[d^N,d^N]);
            end
            t = kron(sy,sy);
            Gt = t*conj(GG)*t;
            G = reshape((GG + Gt')/2,d*ones(1,2*N));
        elseif N == 3
            G = reshape(x,[d,d,d,d,d,d]);
            Gp = permute(G,[3 2 1 6 5 4]);
            if p == 0
                GG = reshape((G-Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^mod(N,2)*Gp)/2,[d^N,d^N]);
            end
            t = kron(kron(sy,sy),sy);
            Gt = t*conj(GG)*t;
            G = reshape((GG + Gt')/2,d*ones(1,2*N));
        elseif N == 4
            G = reshape(x,[d,d,d,d,d,d,d,d]);
            Gp = permute(G,[4 3 2 1 8 7 6 5]);
            if p == 0
                GG = reshape((G-Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^mod(N,2)*Gp)/2,[d^N,d^N]);
            end
            t = kron(kron(kron(sy,sy),sy),sy);
            Gt = t*conj(GG)*t;
            G = reshape((GG + Gt')/2,d*ones(1,2*N));
        end
    elseif podd == 0 && teven == 0
        if N == 1
            G = reshape(x,d*ones(1,2*N));
        elseif N == 2
            G = reshape(x,[d,d,d,d]);
            Gp = permute(G,[2 1 4 3]);
            if p == 0
                GG = reshape((G+Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
            end
            t = kron(sy,sy);
            Gt = t*conj(GG)*t;
            G = reshape((GG - Gt')/2,d*ones(1,2*N));
        elseif N == 3
            G = reshape(x,[d,d,d,d,d,d]);
            Gp = permute(G,[3 2 1 6 5 4]);
            if p == 0
                GG = reshape((G+Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
            end
            t = kron(kron(sy,sy),sy);
            Gt = t*conj(GG)*t;
            G = reshape((GG - Gt')/2,d*ones(1,2*N));
        elseif N == 4
            G = reshape(x,[d,d,d,d,d,d,d,d]);
            Gp = permute(G,[4 3 2 1 8 7 6 5]);
            if p == 0
                GG = reshape((G+Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
            end
            t = kron(kron(kron(sy,sy),sy),sy);
            Gt = t*conj(GG)*t;
            G = reshape((GG - Gt')/2,d*ones(1,2*N));
        end
    elseif podd == 0 && teven == -1
        if N == 1
            G = reshape(x,d*ones(1,2*N));
        elseif N == 2
            G = reshape(x,[d,d,d,d]);
            Gp = permute(G,[2 1 4 3]);
            if p == 0
                GG = reshape((G+Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
            end
            G = reshape(GG,d*ones(1,2*N));
        elseif N == 3
            G = reshape(x,[d,d,d,d,d,d]);
            Gp = permute(G,[3 2 1 6 5 4]);
            if p == 0
                GG = reshape((G+Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
            end
            G = reshape(GG,d*ones(1,2*N));
        elseif N == 4
            G = reshape(x,[d,d,d,d,d,d,d,d]);
            Gp = permute(G,[4 3 2 1 8 7 6 5]);
            if p == 0
                GG = reshape((G+Gp)/2,[d^N,d^N]);
            else
                GG = reshape((G+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
            end
            G = reshape(GG,d*ones(1,2*N));
        end
    elseif podd == -1 && teven == -1
        G = reshape(x,d*ones(1,2*N)); 
    end

    %G = reshape(x,d*ones(1,2*N));
    ALc = conj(AL);
    if N == 1
        g = ncon({G,R,AL,ALc},{[2 3],[5 4],[1 2 5],[1 3 4]},[1 2 3 4 5]);
    elseif N == 2
        g = ncon({G,R,AL,AL,ALc,ALc},{[3 7 2 6],[9 8],[1 3 5],[5 7 9],[1 2 4],[4 6 8]},[1 2 3 4 5 6 7 8 9]);
    elseif N == 3
        g = ncon({G,R,AL,AL,AL,ALc,ALc,ALc},{[3 7 11 2 6 10],[13 12],[1 3 5],[5 7 9],[9 11 13],[1 2 4],[4 6 8],[8 10 12]},[1 2 3 4 5 6 7 8 9 10 11 12 13]);
    elseif N == 4
        g = ncon({G,R,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[2 6 10 14 3 7 11 15],[16 17],[1 2 4],[4 6 8],[8 10 12],[12 14 16],[1 3 5],[5 7 9],[9 11 13],[13 15 17]});
    end
    yt = 2*getapplytoGLp_ncon(AL, G, R, N, D, p, g);
    yt = yt + 2*getapplytoGRp_ncon(AR, G, L, N, D, p, g);
    if N == 1
        yg = ncon({R,AL,ALc},{[3 4],[5 -2 3],[5 -1 4]},[3 4 5],[-1 -2]);
    elseif N == 2
        yg = ncon({R,AL,AL,ALc,ALc},{[8 9],[5 -3 6],[6 -4 8],[5 -1 7],[7 -2 9]},[5 6 7 8 9],[-1 -2 -3 -4]);
    elseif N == 3
        yg = ncon({R,AL,AL,AL,ALc,ALc,ALc},{[12 13],[7 -4 8],[8 -5 10],[10 -6 12],[7 -1 9],[9 -2 11],[11 -3 13]},[7 8 9 10 11 12 13],[-1 -2 -3 -4 -5 -6]);
    elseif N == 4
        yg = ncon({R,AL,AL,AL,AL,ALc,ALc,ALc,ALc},{[9 10],[17 -5 15],[15 -6 13],[13 -7 11],[11 -8 9],[17 -1 16],[16 -2 14],[14 -3 12],[12 -4 10]});
    end
    yt = yt + getapplytoGLoverlapp_ncon(AL, G, R, N, d, p, g, yg);
    yt = yt + getapplytoGRoverlapp_ncon(AR, G, L, N, d, p, g, yg);
    yt = yt + getapplytoGoverlap_ncon(AL, G, R, N, g, yg);
    yt = yt + getapplytoGLoverlapp_ncon_transpose(AL, G, R, N, d, p, g, yg);
    yt = yt + getapplytoGRoverlapp_ncon_transpose(AR, G, L, N, d, p, g, yg);
    yt = yt + getapplytoGoverlap_ncon_transpose(AL, G, R, N, g, yg);
    y = reshape(yt,[d^(2*N),1]);

    %if podd == 1 && teven == 1
    %    if N == 1
    %        y = reshape(yt,[d^(2*N),1]);
    %    elseif N == 2
    %        Gp = permute(yt,[2 1 4 3]);
    %        if p == 0
    %            GG = reshape((yt-Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^mod(N,2)*Gp)/2,[d^N,d^N]);
    %        end
    %        t = kron(sy,sy);
    %        Gt = t*conj(GG)*t;
    %        y = reshape((GG + Gt')/2,[d^(2*N),1]);
    %    elseif N == 3
    %        Gp = permute(yt,[3 2 1 6 5 4]);
    %        if p == 0
    %            GG = reshape((yt-Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^mod(N,2)*Gp)/2,[d^N,d^N]);
    %        end
    %        t = kron(kron(sy,sy),sy);
    %        Gt = t*conj(GG)*t;
    %        y = reshape((GG + Gt')/2,[d^(2*N),1]);
    %    elseif N == 4
    %        Gp = permute(yt,[4 3 2 1 8 7 6 5]);
    %        if p == 0
    %            GG = reshape((yt-Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^mod(N,2)*Gp)/2,[d^N,d^N]);
    %        end
    %        t = kron(kron(kron(sy,sy),sy),sy);
    %        Gt = t*conj(GG)*t;
    %        y = reshape((GG + Gt')/2,[d^(2*N),1]);
    %    end    
    %elseif podd == 0 && teven == 0
    %    if N == 1
    %        y = reshape(yt,[d^(2*N),1]);
    %    elseif N == 2
    %        Gp = permute(yt,[2 1 4 3]);
    %        if p == 0
    %            GG = reshape((yt+Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
    %        end
    %        t = kron(sy,sy);
    %        Gt = t*conj(GG)*t;
    %        y = reshape((GG - Gt')/2,[d^(2*N),1]);
    %    elseif N == 3
    %        Gp = permute(yt,[3 2 1 6 5 4]);
    %        if p == 0
    %            GG = reshape((yt+Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
    %        end
    %        t = kron(kron(sy,sy),sy);
    %        Gt = t*conj(GG)*t;
    %        y = reshape((GG - Gt')/2,[d^(2*N),1]);
    %    elseif N == 4
    %        Gp = permute(yt,[4 3 2 1 8 7 6 5]);
    %        if p == 0
    %            GG = reshape((yt+Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
    %        end
    %        t = kron(kron(kron(sy,sy),sy),sy);
    %        Gt = t*conj(GG)*t;
    %        y = reshape((GG - Gt')/2,[d^(2*N),1]);
    %    end 
    %elseif podd == 0 && teven == -1
    %    if N == 1
    %        y = reshape(yt,[d^(2*N),1]);
    %    elseif N == 2
    %        Gp = permute(yt,[2 1 4 3]);
    %        if p == 0
    %            GG = reshape((yt+Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
    %        end
    %        y = reshape(GG,[d^(2*N),1]);
    %    elseif N == 3
    %        Gp = permute(yt,[3 2 1 6 5 4]);
    %        if p == 0
    %            GG = reshape((yt+Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
    %        end
    %        y = reshape(GG,[d^(2*N),1]);
    %    elseif N == 4
    %        Gp = permute(yt,[4 3 2 1 8 7 6 5]);
    %        if p == 0
    %            GG = reshape((yt+Gp)/2,[d^N,d^N]);
    %        else
    %            GG = reshape((yt+(-1)^(mod(N,2)+1)*Gp)/2,[d^N,d^N]);
    %        end
    %        y = reshape(GG,[d^(2*N),1]);
    %    end 
    %elseif podd == -1 && teven == -1
    %    y = reshape(yt,[d^(2*N),1]);
    %end
end
