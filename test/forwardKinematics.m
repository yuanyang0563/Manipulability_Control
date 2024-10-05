function [x, J] = forwardKinematics(q)

armLength = 4;

x = zeros(2,1);
x(1) = cos(q(1))+cos(q(1)+q(2))+cos(q(1)+q(2)+q(3));
x(2) = sin(q(1))+sin(q(1)+q(2))+sin(q(1)+q(2)+q(3));
x = armLength.*x;

J = zeros(2,3);
J(1,1) = -sin(q(1))-sin(q(1)+q(2))-sin(q(1)+q(2)+q(3));
J(1,2) = -sin(q(1)+q(2))-sin(q(1)+q(2)+q(3));
J(1,3) = -sin(q(1)+q(2)+q(3));
J(2,1) =  cos(q(1))+cos(q(1)+q(2))+cos(q(1)+q(2)+q(3));
J(2,2) =  cos(q(1)+q(2))+cos(q(1)+q(2)+q(3));
J(2,3) =  cos(q(1)+q(2)+q(3));
J = armLength.*J;

end