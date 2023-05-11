function sq = structurefactor(O1, O2, AL, AR, AC, L, R, d, D, p)
    O1e = 0.0;
    O2e = 0.0;
    for i = 1:d
        for ii = 1:d
            O1e = O1e + O1(i,ii)*trace(AC{i}*AC{ii}');
            O2e = O2e + O2(i,ii)*trace(AC{i}*AC{ii}');
        end
    end
    O1n = O1 - O1e*eye(d);
    O2n = O2 - O2e*eye(d);
    sq1 = 0.0;
    for i = 1:d
        for ii = 1:d
            for iii = 1:d
                sq1 = sq1 + O1n(i,ii)*O2n(ii,iii)*trace(AC{i}*AC{iii}');
            end
        end
    end
    O1L = zeros(D);
    O1R = zeros(D);
    O2L = zeros(D);
    O2R = zeros(D);
    for i = 1:d
        for ii = 1:d
            O1L = O1L + O1n(i,ii)*AL{ii}'*AL{i};
            O2R = O2R + O2n(i,ii)*AC{i}*AC{ii}';
            O2L = O2L + O2n(i,ii)*AC{ii}'*AC{i};
            O1R = O1R + O1n(i,ii)*AR{i}*AR{ii}';
        end
    end
    afun = @(x) ELLregcomp(x, AL, R, d, D, p);
    sq2 = exp(-1i*p)*trace(reshape(gmres(afun,reshape(O1L,[D^2,1]),[],1E-14,1000),[D,D])*O2R);
    afun = @(x) ERRregcomp(x, AR, L, d, D, p);
    sq3 = exp(1i*p)*trace(O2L*reshape(gmres(afun,reshape(O1R,[D^2,1]),[],1E-14,1000),[D,D]));
    sq = sq1+sq2+sq3;
end

function y = ELLregcomp(x, AL, R, d, D, p)
    FL = reshape(x,[D,D]);
    yt = FL;
    for i = 1:d
        yt =  yt - exp(-1i*p)*AL{i}'*FL*AL{i};
    end
    yt = yt + (exp(-1i*p)*trace(FL*R))*eye(D);
    y = reshape(yt,[D^2,1]);
end

function y = ERRregcomp(x, AR, L, d, D, p)
    FR = reshape(x,[D,D]);
    yt = FR;
    for i = 1:d
        yt =  yt - exp(1i*p)*AR{i}*FR*AR{i}';
    end
    yt = yt + (exp(1i*p)*trace(L*FR))*eye(D);
    y = reshape(yt,[D^2,1]);
end
