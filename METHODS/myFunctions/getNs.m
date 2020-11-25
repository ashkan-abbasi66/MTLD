function Ns = getNs(imn)
    % Get the number of decomposition levels
    %
    % Borrowed from FusedKSVD package
    
    [R,C] = size(imn);
    lmin = min(R,C);
    if lmin>550
        if lmin > 1200
            Ns = 3;
        else
            Ns=2;
        end
    else
        Ns=1;
    end

end