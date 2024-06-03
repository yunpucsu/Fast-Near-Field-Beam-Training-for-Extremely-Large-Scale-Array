function at = far_field_mainfold(Nt,d,f,theta0)
    c = 3e8;
    lambda = c/f;
    at = exp(-1i*2*pi/lambda*[0:Nt-1]'*d*sin(theta0));
    %n = 0:1:Nt-1;
    %at = exp(-1j*2*pi/lambda.*n*d*sin(theta0));
    at = at / norm(at);
end