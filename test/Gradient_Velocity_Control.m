function Gradient_Velocity_Control

% This function optimizes the velocity manipulability of a velocity-controlled dual-arm system using the 
% gradient-based approach in the nullspace.

addpath('../fcts/');

dt = 1e-3;
nbIter = 1000;
Kt = 10;
Kp = 100;

nbDOFs = 3;
armLength = 4;

q10 = [5*pi/6; -2*pi/3; -pi/4]; 
q20 = [1*pi/6;  2*pi/3;  pi/4];

xd = [0;6];

q1t = q10;
q2t = q20;
it = 1;

figure('position',[5 50 1000 567],'color','w');
hold on;
xlim([-10,10]);
ylim([-2,8]);
xlabel('$x$','fontsize',20,'Interpreter','latex');
ylabel('$y$','fontsize',20,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);
set(gca,'Position',[0.1 0.2 0.85 0.75]);
set(gca,'LineWidth',1);

while(it<=nbIter)

    q1t_track(:,it) = q1t;
    q2t_track(:,it) = q2t;
    f1t_track(it) = Manipulability(q1t);
    f2t_track(it) = Manipulability(q2t);

    [x1t,J1t] = forwardKinematics(q1t);
    x1t = x1t-[3.8637;0];
    [x2t,J2t] = forwardKinematics(q2t);
    x2t = x2t+[3.8637;0];

    dx_track(:,it) = x1t-x2t;

    dxr1 = Kt*(xd-x1t)+Kp*(x2t-x1t);
    J1t_pinv = J1t'/(J1t*J1t');
    dq_t1 = J1t_pinv*dxr1;
    f1_grad = ManipulabilityGradient(q1t);
    Ns1 = eye(nbDOFs)-J1t_pinv*J1t;
    dq_ns1 = Ns1*f1_grad;

    dxr2 = Kp*(x1t-x2t);
    J2t_pinv = J2t'/(J2t*J2t');
    dq_t2 = J2t_pinv*dxr2;
    f2_grad = ManipulabilityGradient(q2t);
    Ns2 = eye(nbDOFs)-J2t_pinv*J2t;
    dq_ns2 = Ns2*f2_grad;

	p1 = plotArm(q1t, ones(nbDOFs,1)*armLength, [-3.8637; 0; 1], .2, [0,0,0]);
    p2 = plotArm(q2t, ones(nbDOFs,1)*armLength, [+3.8637; 0; 1], .2, [0,0,0]);
    h1 = plotGMM(x1t, 1e-2*J1t*J1t', [1 0 0], .1, '-', 2, 1);
    h2 = plotGMM(x2t, 1e-2*J2t*J2t', [0 0 1], .1, '-', 2, 1);
    if(it==1)
		plotGMM(x1t, 1e-2*J1t*J1t', [1 0 0], .1, '--', 2, 1);
        plotGMM(x2t, 1e-2*J2t*J2t', [0 0 1], .1, '--', 2, 1);
    end
    drawnow;
	
	q1t = q1t+(dq_t1+dq_ns1)*dt;
    q2t = q2t+(dq_t2+dq_ns2)*dt;
	it = it+1;
    
    if it<=nbIter
        delete(p1);
        delete(p2);
        delete(h1);
        delete(h2);
    end

end

figure('position',[105 150 1000 567],'color','w');
list = [1, 5, 10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 500, 750, 1000];
for i = 1:length(list)
    it = list(i);
	colTmp = [1,1,1]*(1-it/nbIter);
	plotArm(q1t_track(:,it), ones(nbDOFs,1)*armLength, [-3.8637; 0; it/nbIter], .2, colTmp);
    plotArm(q2t_track(:,it), ones(nbDOFs,1)*armLength, [+3.8637; 0; it/nbIter], .2, colTmp);
end
hold on;
xlim([-10,10]);
ylim([-2,8]);
xlabel('$x$','fontsize',20,'Interpreter','latex');
ylabel('$y$','fontsize',20,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);
set(gca,'Position',[0.1 0.2 0.85 0.75]);
set(gca,'LineWidth',1);

figure('position',[205 250 1000 567],'color','w');
plot_f1 = plot(1:nbIter,f1t_track,'Color','r','LineStyle','-','LineWidth',2);
hold on;
plot_f2 = plot(1:nbIter,f2t_track,'Color','b','LineStyle','-','LineWidth',2);
hold on;
leg_f = legend([plot_f1 plot_f2],{'$f_1$','$f_2$'},'Location','southeast','FontSize',20,'Orientation','Horizontal','NumColumns',1);
set(leg_f, 'interpreter', 'latex','color','none');
xlim([0,nbIter]);
ylim([25,45]);
xlabel('Iterations','fontsize',20,'Interpreter','latex');
ylabel('Manipulability','fontsize',20,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);
set(gca,'Position',[0.1 0.2 0.85 0.75]);
set(gca,'LineWidth',1);

figure('position',[305 350 1000 567],'color','w');
plot_dx = plot(1:nbIter,dx_track(1,:),'Color','r','LineStyle','-','LineWidth',2);
hold on;
plot_dy = plot(1:nbIter,dx_track(2,:),'Color','b','LineStyle','-','LineWidth',2);
hold on;
leg_dxy = legend([plot_dx plot_dy],{'$x_1-x_2$','$y_1-y_2$'},'Location','northeast','FontSize',20,'Orientation','Horizontal','NumColumns',1);
set(leg_dxy, 'interpreter', 'latex','color','none');
xlim([0,nbIter]);
xlabel('Iterations','fontsize',20,'Interpreter','latex');
ylabel('$$\textbf{x}_1-\textbf{x}_2$$','fontsize',20,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',20);
set(gca,'Position',[0.1 0.2 0.85 0.75]);
set(gca,'LineWidth',1);

end