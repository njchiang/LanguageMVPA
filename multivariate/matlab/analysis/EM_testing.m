%EM testing
x0 = [1/3 1/3 1/3 0 1 0 1 0 1 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
A = -[1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];
b = [0;0;0;0];
Aeq = [ 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];
beq = [1;1;1;1];
% initialize theta to uniform and standard normal
% calculate tau's
ts = intVox(1:32);

% synthetic ts
ts = repmat(1:8, 1, 4);
label = repmat(1:8, 1, 4);
c = repmat(eye(max(label)), 4, 1);
tau = zeros(3, length(ts));
prevTau = tau+1;
ts = ts - mean(ts);
x = x0;
thresh = 1;
nIter = 0;
% need to include indicator for z and c.
% z is whether the voxel over or undershoots
% c is the category label
while abs(sum(sum(tau - prevTau))) > thresh
    nIter = nIter+1;
    if nIter > 5000
        disp('Did not converge')
        break;
    end
    prevTau = tau;
    tmpLabs=x(10:17);
    tau(1,:) = eps+x(1).*tmpLabs(label).*normpdf(ts, x(4), x(5)); % the x(10) is wrong...
    tmpLabs=x(18:25);
    tau(2,:) = eps+x(2).*tmpLabs(label).*normpdf(ts, x(6), x(7));
    tmpLabs=x(26:33);
    tau(3,:) = eps+x(3).*tmpLabs(label).*normpdf(ts, x(8), x(9));
    sumTau=tau(1,:)+tau(2,:)+tau(3,:);
    tau(1,:) = tau(1,:)./sumTau;
    tau(2,:) = tau(2,:)./sumTau;
    tau(3,:) = tau(3,:)./sumTau;
    
    fun = @(x)sum(...
        tau(1,:).*(log(x(1))+ log(x(10:17))*c'+log(normpdf(ts,x(4),x(5))))+ ...
        tau(2,:).*(log(x(2))+ log(x(18:25))*c'+log(normpdf(ts,x(6),x(7))))+ ...
        tau(3,:).*(log(x(3))+ log(x(26:33))*c'+log(normpdf(ts,x(8),x(9)))))
%         tau(1,:).*(log(x(1))+   log(x(10))+log(x(11))+log(x(12))+log(x(13))+log(x(14))+log(x(15))+log(x(16))+log(x(17)))+log(normpdf(ts,x(4),x(5)))...
%        +tau(2,:).*(log(x(2))+log(x(18))+log(x(19))+log(x(20))+log(x(21))+log(x(22))+log(x(23))+log(x(24))+log(x(25)))+log(normpdf(ts,x(6),x(7)))...
%        +tau(3,:).*(log(x(3))+log(x(26))+log(x(27))+log(x(28))+log(x(29))+log(x(30))+log(x(31))+log(x(32))+log(x(33)))+log(normpdf(ts,x(8),x(9))));
    
    x = fmincon(fun, x, A, b, Aeq, beq);
    fun(x)
end
