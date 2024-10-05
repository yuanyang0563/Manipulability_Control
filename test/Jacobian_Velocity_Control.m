function Jacobian_Velocity_Control

% this function optimizes the velocity manipulability of a velocity-controlled dual-arm system by using the
% Jacobian-based approach in the nullspace

addpath('../fcts/');

dt = 1e-3;
nbIter = 2000;
Kt = 10;
Kp = 100;

nbDOFs = 3;
armLength = 4;

L1 = Link('d', 0, 'a', armLength, 'alpha', 0);
robot1 = SerialLink(repmat(L1,nbDOFs,1));

L2 = Link('d', 0, 'a', armLength, 'alpha', 0);
robot2 = SerialLink(repmat(L2,nbDOFs,1));

q1_0 = [5*pi/6; -2*pi/3; -pi/4]; 
q2_0 = [1*pi/6;  2*pi/3;  pi/4];

Me_d = [77.3167 -51.8251; -51.8251 57.7210];
xd = [0;6];

q1t = q1_0;
q2t = q2_0;
it = 1;
p1 = [];
p2 = [];

figure('position',[10 10 1000 450],'color',[1 1 1]);
axis equal;
xlim([-10,10]);
ylim([-2,8]);
xlabel('$x_1$','fontsize',20,'Interpreter','latex');
ylabel('$x_2$','fontsize',20,'Interpreter','latex');

while(it<=nbIter)
    
    delete(p1);
    delete(p2);
	
	J1t = robot1.jacob0(q1t);
    J1t_full = J1t;
    J1t = J1t(1:2,:);
	x1t = robot1.fkine(q1t).t(1:2)-[3.8637;0];
    
	J2t = robot2.jacob0(q2t);
    J2t_full = J2t;
    J2t = J2t(1:2,:);
	x2t = robot2.fkine(q2t).t(1:2)+[3.8637;0];
    
    Me1_ct = J1t*J1t';
    Me2_ct = J2t*J2t';

    Jm1_t = compute_red_manipulability_Jacobian(J1t_full,1:2);
    Jm2_t = compute_red_manipulability_Jacobian(J2t_full,1:2);
	
    dxr1 = Kt*(xd-x1t)+Kp*(x2t-x1t);
    dq_t1 = pinv(J1t)*dxr1;
    M1_diff = logmap(Me_d,Me1_ct);
    dq_t1m = pinv(Jm1_t)*symmat2vec(M1_diff);
    Ns1 = eye(nbDOFs)-pinv(J1t)*J1t;
    dq_ns1 = Ns1*dq_t1m;
    
    dxr2 = Kp*(x1t-x2t);
    dq_t2 = pinv(J2t)*dxr2;
    M2_diff = logmap(Me_d,Me2_ct);
    dq_t2m = pinv(Jm2_t)*symmat2vec(M2_diff);
    Ns2 = eye(nbDOFs)-pinv(J2t)*J2t;
    dq_ns2 = Ns2*dq_t2m; 

	colTmp = [1,1,1] - [.8,.8,.8];
	p1 = plotArm(q1t, ones(nbDOFs,1)*armLength, [-3.8637; 0; 1], .2, colTmp);
    p2 = plotArm(q2t, ones(nbDOFs,1)*armLength, [+3.8637; 0; 1], .2, colTmp);
	
	drawnow;
	
	q1t = q1t + (dq_t1 + dq_ns1)*dt;
    q2t = q2t + (dq_t2 + dq_ns2)*dt;
	it = it + 1;

end

end