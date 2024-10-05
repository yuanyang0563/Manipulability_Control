function f = Manipulability(q)

[~,J] = forwardKinematics(q);
f = sqrt(det(J*J'));

end