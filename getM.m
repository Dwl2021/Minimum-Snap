function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        M_k = zeros(4,8);
        M_k(1,1) = 1;
        M_k(2,2) = 1;
        M_k(3,3) = 2;
        M_k(4,4) = 6;
        t = ts(k);
        M_k = [M_k;getCoeff(t)];      % 8*8
        M = blkdiag(M, M_k);    % 8M*8M 
    end
end