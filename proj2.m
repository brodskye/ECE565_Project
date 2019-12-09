options = optimoptions(@fminunc, 'Display','off');


%Gaussian:
gaussian = 0;

if gaussian
    vals1 = []; %for MLE
    vals2 = []; %for MAP
    mses = [];
    for i=1:100
        MSE1=0; %for MLE
        MSE2=0; %for MAP
        for j=1:200
            rng('default');
            A = i;
            n=5000;
            sigma=10;
            y = sign(repmat(A,n,1)+ normrnd(0,sigma,[n,1]));

            %MLE:
            if (1/n) * sum((y)) == 1
                A_mle = -sigma * sqrt(2)*erfinv(-.9999999);
            elseif (1/n) * sum((y)) == -1
                A_mle = -sigma * sqrt(2)*erfinv(.9999999);
            else
                A_mle = -sigma * sqrt(2)*erfinv(-(1/n) * sum((y)));
            end

            %MAP:
            obj = @(a) -((sum((1+y)/2))*log(1-phi(-a/sigma)) + (sum((1-y)/2)) * log(phi(-a/sigma)) -a^2/((2*sigma)^2));
            %[A_map, fval] = fminunc(obj, 0, options);
            
            MSE1 = MSE1 + (A-A_mle)^2;
            MSE2 = MSE2 + (A-A_map)^2;
        end
        Ey = 1*(1-phi(-A/sigma)) + (-1) * (phi(-A\sigma));
        
        
        
        
        %ratio = A/sigma
        mses = [mses; MSE1];
        MSE1 = MSE1/200;
        MSE2 = MSE2/200;
        vals1 = [vals1; MSE1];
        vals2 = [vals2; MSE2];
    end
    figure(1)
    loglog([0.1:0.1:10], vals1)
    figure(2)
    %loglog([0.1:0.1:10], vals2)
    loglog([0.1:0.1:10], mses)


end


%Uniform:
uniform = 1; 

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
            %obj = @(a) -((sum((1+y)/2))*log(1-((-a+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + (sum((1-y)/2))* log(((-a+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + log(-a+sqrt(3)*sigma) - log(2*sqrt(3)*sigma));
            obj = @(a) -((sum((1+y)/2))*log(1-((-a+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + (sum((1-y)/2))* log(((-a+sqrt(3)*sigma)/(2*sqrt(3)*sigma))) + log(-a+sqrt(3)*sigma));

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
        %MSE2 = MSE2/200;
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

function s = phi(x)
    s = 1/2 + 1/2 * erf(x/sqrt(2));
end