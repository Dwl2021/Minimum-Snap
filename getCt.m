function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    Ct = zeros(8*n_seg,4*n_seg+4);
    Ct(1:4,1:4) = eye(4);
    for k = 1:n_seg
        if k ~= n_seg
            Ctk = zeros(4,4*n_seg+4);
            % 第k个轨迹的终点
            Ctk(1, 4+k) = 1;
            Ctk(2, n_seg+7+3*(k-1)+1) = 1;
            Ctk(3, n_seg+7+3*(k-1)+2) = 1;
            Ctk(4, n_seg+7+3*(k-1)+3) = 1;
            Ctk = [Ctk;Ctk];
            Ct(5+8*(k-1):5+8*(k-1)+7,:) = Ctk;
        else
            Ctk = zeros(4,4*n_seg+4);
            Ctk(1,n_seg+4) = 1;
            Ctk(1,n_seg+5) = 1;
            Ctk(1,n_seg+6) = 1;
            Ctk(1,n_seg+7) = 1;
            Ct(8*n_seg-3:8*n_seg,:) = Ctk;
        end
        
    end

end
