function [dict, label] = PolarCodeBook(N, d, lambda, beta, rho_min, rho_max)
    c = 3e8;s = 1;
    D = s * N;
    theta = -1 + 2/D : 2/D : 1;
    dict = cell( D, 1 );
    label = cell( D, 1);
    Zmax = (N * d)^2 / 2 / lambda / beta^2;
    kmax = floor(Zmax/rho_min);
    for idx = 1:D
        Z = (N * d)^2 * ( 1 - theta(idx)^2) / 2 / lambda / beta^2;
        kmax = floor(Z/rho_min);
        kmin = floor(Z/rho_max) + 1;
        
        r = zeros(1, kmax - kmin + 2);
        r(:,1) = (N * d)^2 * 2 / lambda; 
        r(:,2:end) = Z./(kmin:kmax);
        
        dict{idx} = zeros(N, kmax + 1);
        label{idx} = zeros(2, kmax + 1);
        for t = 1 : kmax - kmin + 2
            dict{idx}(:, t) = polar_domain_manifold( N, d, c/lambda, r(t), asin(theta(idx)) );
            label{idx}(:, t) = [theta(idx), r(t)]';
        end
    end
    dict = merge(dict, D, N); % converted to matrix  
    label = merge(label, D, 2); % converted to matrix
end

function B = merge(A, D, Q)
    S = zeros(1, D);
    for idx = 1 : D
        S(idx) = size(A{idx}, 2);
    end
    B = zeros(Q, sum(S));
    for idx = 1 : D
        B(:, sum(S(1:idx)) - S(idx) + 1: sum(S(1:idx))) = A{idx};
    end
end