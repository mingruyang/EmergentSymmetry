function c = correlator1site(O1, O2, AL, AC, d, D, r)
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
    c = zeros(r,1);
    for i = 1:d
        for ii = 1:d
            for iii = 1:d
                c(1) = c(1) + O1n(i,ii)*O2n(ii,iii)*trace(AC{i}*AC{iii}');
            end
        end
    end
    O1L = zeros(D);
    O2R = zeros(D);
    for i = 1:d
        for ii = 1:d
            O1L = O1L + O1n(i,ii)*AL{ii}'*AL{i};
            O2R = O2R + O2n(i,ii)*AC{i}*AC{ii}';
        end
    end
    for i = 2:r
        c(i) = trace(O1L*O2R);
        O1L = leftApplyT(O1L, AL, d, D);
    end
end

function y = leftApplyT(x, AL, d, D)
    y = zeros(D,D);
    for i = 1:d
        y =  y + AL{i}'*x*AL{i};
    end
end
