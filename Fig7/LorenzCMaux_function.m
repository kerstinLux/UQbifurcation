function Y = LorenzCMaux_function(X)
Y = X(:,1).*(1+X(:,1)).^(-1);
% Y = X(:,1);