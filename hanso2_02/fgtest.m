function [ f , g ] = fgtest( b , pars )
X=pars.X;
Y=pars.Y;
lambda=pars.lambda;
f = (norm(Y-X*b,2)+lambda*(norm(b,1)*norm(b,2))^(1/2));
g = ((-1)*(X'*(Y-X*b))/(norm(Y-X*b,2))+lambda*(1/2)*sign(b)*sqrt(norm(b,2)/norm(b,1))+lambda/2*b/norm(b,2)*sqrt(norm(b,1)/norm(b,2)));
end

