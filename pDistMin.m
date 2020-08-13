function t = pDistMin(x, muA, muB, sigmaA, sigmaB, rho)
% From https://www.gwern.net/docs/conscientiousness/2008-nadarajah.pdf
t1 = normpdf(x, muA, sigmaA);
tt = rho*(x - muA)/(sigmaA*sqrt(1 - rho^2));
tt = tt - (x - muB)/(sigmaB*sqrt(1 - rho^2));
t1 = t1.*normcdf(tt);
t2 = normpdf(x, muB, sigmaB);
tt = rho*(x-muB)/(sigmaB*sqrt(1 - rho^2));
tt = tt - (x - muA)/(sigmaA*sqrt(1 - rho^2));
t2 = t2.*normcdf(tt);
t = t1 + t2;