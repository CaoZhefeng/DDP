function dmudx = grad_gaussian(x,u,gprMd)
tic;
mu = predict(gprMd,[x,u]);

A = gprMd.Alpha;
n = gprMd.NumObservations;
kernel_p = exp(gprMd.KernelInformation.KernelParameters);
% compute the Kernal Matrix
Kernel_vector = exp(-0.5*(ones(n,1)*[x,u]-gprMd.X).^2*(1./kernel_p(1:end-1,1).^2))';
x = [x,u];

% for i =1:size(x,2)
%     for j = 1:n
%         grad_Mat(i,j) = Kernel_vector(j)*(gprMd.X(j,i)-x(i))/kernel_p(i)^2;
%     end
% end
grad_Mat=((Kernel_vector'*ones(1,size(x,2))).*(gprMd.X-ones(n,1)*x)./(ones(n,1)*kernel_p(1:end-1)').^2)';
dmudx = grad_Mat*A;
if strcmp(gprMd.BasisFunction,'Linear')
    dmudx = dmudx+gprMd.Beta(2:end);
end
t2=toc;
end