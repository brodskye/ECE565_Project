options = optimoptions(@fminunc, 'Display','off');


%Gaussian:
gaussian = 1;

if gaussian
    vals1 = []; %for MLE
    vals2 = []; %for MAP
    vals3 = []; %for CRLB

    mses = [];
    
    A_vec=exp(linspace(log(1),log(100),10));
    for i=1:length(A_vec)
        MSE1=0; %for MLE
        MSE2=0; %for MAP
         A = A_vec(i);
         R=20000;
        for j=1:R
            %
            % rng('default');
           
            n=5000;
            sigma=10;
            y = sign(repmat(A,n,1)+ normrnd(0,sigma,[n,1]));
%qq(j)=sum(y);
            %MLE:
            if (1/n) * sum((y)) == 1
                A_mle = -sigma * sqrt(2)*erfinv(-.9999999);
            elseif (1/n) * sum((y)) == -1
                A_mle = -sigma * sqrt(2)*erfinv(.9999999);
            else
                A_mle = -sigma * sqrt(2)*erfinv(-(1/n) * sum((y)));
            end

            %MAP:
            obj = @(a) -((sum((1+y)/2))*log(1-PHI(-a/sigma)) + (sum((1-y)/2)) * log(PHI(-a/sigma)) -a^2/((2*sigma)^2));
            %[A_map, fval] = fminunc(obj, 0, options);
            
            MSE1 = MSE1 + (A-A_mle)^2;
            %MSE2 = MSE2 + (A-A_map)^2;
        end
 %       hist(qq),pause
        Ey = 1*(1-PHI(-A/sigma)) + (-1) * (PHI(-A/sigma));
        %FIM = ((sum((1-Ey)/2))*(PHI(-A/sigma)*(1/sigma^2)*phip(-A/sigma) - (1/sigma^2)*phi(-A/sigma)^2)/(PHI(-A/sigma)^2)) - ...
            %((sum((1+Ey)/2))*(1-PHI(-A/sigma)*(1/sigma^2)*phip(-A/sigma) - (1/sigma^2)*phi(-A/sigma)^2)/(1-PHI(-A/sigma)^2))
        
        FIM = (1/(sigma^2))*(phi(-A/sigma))^2 * (1/((PHI(-A/sigma)*(1-PHI(-A/sigma)))));
        CRLB = 1/FIM;
        CRLB = CRLB/n;
        ratio = A/sigma
        mses = [mses; MSE1];
        MSE1 = MSE1/R;
        MSE2 = MSE2/R;
        vals1 = [vals1; MSE1];
        vals2 = [vals2; MSE2];
        vals3 = [vals3; CRLB];
        
    end
    figure(1)
    loglog(A_vec/sigma, vals1)
    hold on;
    %figure(2)
    loglog(A_vec/sigma, vals2)
    loglog(A_vec/sigma, vals3)

    %loglog([0.1:0.1:10], mses)


end


%Uniform:
uniform = 0; 

if uniform
    vals1 = []; %for MLE
    vals2 = []; %for MAP
    vals3 = []; %for CRLB
    %mses = [];
    for i=1:100
        MSE1=0; %for MLE
        MSE2=0; %for MAP
        for j=1:200
            A = i;
            n = 5000;
            sigma = 10;
            pd = makedist('Uniform', 'lower',-sqrt(3)*sigma, 'upper',sqrt(3)*sigma);

            y = sign(repmat(A,n,1)+ random(pd, [n, 1]));

            %MLE:
            A_mle = (sqrt(3)*sigma * sum(y))/n;
            
            %MAP:
            obj = @(a) -((sum((1+y)/2))*log(1-((-a+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + ...
                (sum((1-y)/2))* log(((-a+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + log(-a+sqrt(3)*sigma) - log(2*sqrt(3)*sigma));

            [A_map, fval] = fminunc(obj, 0, options);
            
            MSE1 = MSE1 + (A-A_mle)^2;
            MSE2 = MSE2 + (A-A_map)^2;
        end
        %CRLB:
        Ey = 1*(1-((-A+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + (-1) * (((-A+sqrt(3)*sigma)/(2*sqrt(3)*sigma)));
        FIM = (-(-1*sum((1-Ey)/2)/((sqrt(3)*sigma-A)^2)) + sum((1+Ey)/2)/((sqrt(3)*sigma+A)^2));
        CRLB = 1/FIM
        CRLB = CRLB/n;
        ratio = A/sigma
        MSE1 = MSE1/200;
        MSE2 = MSE2/200;
        vals1 = [vals1; MSE1];
        vals2 = [vals2; MSE2];
        vals3 = [vals3; CRLB];
        %mses = [mses; MSE1];
    end
    figure(1)
    loglog([0.1:0.1:10], vals1)
    %figure(2)
    hold on;
    loglog([0.1:0.1:10], vals2)
    %figure(3)
    loglog([0.1:0.1:10], vals3)
    %figure(2)
    %loglog([0.1:0.1:10], mses);
end

function s = PHI(x)
    s = 1/2 + 1/2 * erf(x/sqrt(2));
end

function s = phi(x)
    s = (exp(-(x.^2)/2))/(sqrt(2*pi));
end

function s = phip(x)
    s = (-x).*phi(x);
end

