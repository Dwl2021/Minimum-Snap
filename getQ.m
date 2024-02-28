function Q = getQ(n_seg, n_order, T)
    Q = [];
    fac = @(x) x*(x-1)*(x-2)*(x-3);
    for k = 1:n_seg
        %#####################################################
        % STEP 1.1: calculate Q_k of the k-th segment 
        Q_k=zeros(n_order+1,n_order+1);
        for i=0:n_order
            for l=0:n_order
                if (i<4) || (l<4)
                    continue;
                else
                    Q_k(i+1,l+1)=fac(i)*fac(l)/(i+l-7)*T(k)^(i+l-7);
                end
            end
        end
        Q = blkdiag(Q, Q_k);
    end
end
