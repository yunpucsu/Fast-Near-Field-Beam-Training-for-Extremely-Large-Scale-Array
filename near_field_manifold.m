function  at = near_field_manifold( Nt, d, f, r0, theta0 )
    c = 3e8;
    nn = [-(Nt-1)/2:1:(Nt-1)/2];
    r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sin(theta0));
    at = exp(-1j*2*pi*f*(r - r0)/c)/sqrt(Nt);
end
