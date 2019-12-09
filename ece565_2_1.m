%sources
rng(0,'twister');

noisy = 1;

s1 = [1;1];
s2 = [1;-1];
s3 = [-1;-1];
%s4 = [-1;1];
pos = [0.3;-0.2]; %true location

theta = [0;0];

estimates = [theta];

while i < 19
if noisy
    y = [h(s1, pos); h(s2, pos); h(s3, pos)] + noise() - ...
        [h(s1, theta); h(s2, theta); h(s3, theta)];
    
    %y = [h(s1, pos); h(s2, pos); h(s3, pos); h(s4, pos)] + noise() - ...
        %[h(s1, theta); h(s2, theta); h(s3, theta); h(s4, theta)];

else
    y = [h(s1, pos); h(s2, pos); h(s3, pos)] - ...
        [h(s1, theta); h(s2, theta); h(s3, theta)];
    %y = [h(s1, pos); h(s2, pos); h(s3, pos); h(s4, pos)] - ...
        %[h(s1, theta); h(s2, theta); h(s3, theta); h(s4, theta)];
end
H = [dh(s1, theta);dh(s2, theta);dh(s3, theta)];
%H = [dh(s1, theta);dh(s2, theta);dh(s3, theta);dh(s4, theta)];

theta = theta + H\y;
estimates = [estimates, theta];
i = i+1;
end

function dist = h(s, theta) % s: source node, theta: opt variable.
    dist = sqrt((theta(1)-s(1))^2 + (theta(2)-s(2))^2);
end

function ddist = dh(s, theta)
    ddist = [(theta(1)-s(1))/(sqrt((theta(1)-s(1))^2 + (theta(2)-s(2))^2)), ...
        (theta(2)-s(2))/(sqrt((theta(1)-s(1))^2 + (theta(2)-s(2))^2))];
end

function n = noise()
    n = random('normal', 0, 0.02, 3, 1);
end

%{
PART C:
4 nodes, noiseless:
Columns 1 through 10

         0    0.2982    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000
         0   -0.2012   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000

  Columns 11 through 20

    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000
   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000

4 nodes, noise:
olumns 1 through 10

         0    0.2715    0.3051    0.2715    0.3013    0.3142    0.3221    0.2954    0.3031    0.2949
         0   -0.2141   -0.2164   -0.2389   -0.2013   -0.1894   -0.2204   -0.1913   -0.1959   -0.2218

  Columns 11 through 20

    0.3035    0.3073    0.3081    0.2795    0.3248    0.3080    0.3189    0.2717    0.2904    0.2785
   -0.1895   -0.2009   -0.1974   -0.2071   -0.2069   -0.1987   -0.1940   -0.2167   -0.1828   -0.1876

PART D:
3 nodes, noiseless:
  Columns 1 through 10

         0    0.2968    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000
         0   -0.1998   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000

  Columns 11 through 20

    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000    0.3000
   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000   -0.2000

3 nodes, noise:
  Columns 1 through 10

         0    0.2511    0.2780    0.3357    0.3390    0.3061    0.3190    0.2831    0.3138    0.2923
         0   -0.1937   -0.2099   -0.1767   -0.2317   -0.2057   -0.1937   -0.1933   -0.1725   -0.1885

  Columns 11 through 20

    0.2890    0.2991    0.3563    0.3239    0.3060    0.2849    0.3090    0.2714    0.3033    0.2722
   -0.1965   -0.2275   -0.2228   -0.2078   -0.1853   -0.2024   -0.1992   -0.2001   -0.2208   -0.1838
%}
