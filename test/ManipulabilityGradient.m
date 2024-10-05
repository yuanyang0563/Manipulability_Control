function f_grad = ManipulabilityGradient(q)

[~,J] = forwardKinematics(q);
pJ = partialJacobian(q);
J_pinv = J'/(J*J');
f = sqrt(det(J*J'));
f_grad = zeros(3,1);
for i=1:3
    f_grad(i) = trace(pJ(:,:,i)*J_pinv);
end
f_grad = f.*f_grad;

end