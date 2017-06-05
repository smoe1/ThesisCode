function x=besselzero(n,k,kind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% besselzero.m
%
% Find first k positive zeros of the Bessel function J(n,x) or Y(n,x) 
% using Halley's method.
%
% Written by: Greg von Winckel - 01/25/05
% Contact: gregvw(at)chtm(dot)unm(dot)edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k3=3*k;

    x=zeros(k3,1);

    for j=1:k3
        
        % Initial guess of zeros 
        x0=2.40482555769*j;
        x0=1+sqrt(2)+(j-1)*pi+n+n^0.4;
        
        % Do Halley's method
        x(j)=findzero(n,x0,kind);

        if x(j)==inf
            error('Bad guess.');
        end
        
    end

    x=sort(x);
    dx=[1;abs(diff(x))];
    x=x(dx>1e-8);

    x=x(1:k);
end

function x=findzero(n,x0,kind)

    n1=n+1;     n2=n*n;

    % Tolerance
    tol=1e-15;

    % Maximum number of times to iterate
    MAXIT=10000;

    % Initial error
    err=1;

    iter=0;

    while abs(err)>tol & iter<MAXIT
        
        switch kind
            case 1
                a=besselj(n,x0);    
                b=besselj(n1,x0);   
            case 2
                a=bessely(n,x0);
                b=bessely(n1,x0);
        end
                
        x02=x0*x0;
        
        err=2*a*x0*(n*a-b*x0)/(2*b*b*x02-a*b*x0*(4*n+1)+(n*n1+x02)*a*a);
        
        x=x0-err;
        x0=x;
        iter=iter+1;
        
    end

    if iter>MAXIT-1
        warning('Failed to converge to within tolerance. ',...
                'Try a different initial guess');
        x=inf;    
    end

end
