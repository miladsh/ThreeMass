function gait = doMPC(constants, dynamics, robot, gait)
axRows = size(dynamics.Ax{1},1);
axCols = size(dynamics.Ax{1},2);
bxRows = size(dynamics.Bx{1},1);
bxCols = size(dynamics.Bx{1},2);


ayRows = size(dynamics.Ay{1},1);
ayCols = size(dynamics.Ay{1},2);
byRows = size(dynamics.By{1},1);
byCols = size(dynamics.By{1},2);

dynamics.PcomxX = dynamics.Px(1:axRows:end,1:axRows);
dynamics.PcomxU = dynamics.Px(1:axRows:end,axRows+1:end);

dynamics.PcomyX = dynamics.Py(1:axRows:end,1:axRows);
dynamics.PcomyU = dynamics.Py(1:ayRows:end,ayRows+1:end);



dynamics.dcmxX = dynamics.Px(2:axRows:end,1:axRows);
dynamics.dcmxU = dynamics.Px(2:axRows:end,axRows+1:end);


dynamics.dcmyX = dynamics.Py(2:ayRows:end,1:ayRows);
dynamics.dcmyU = dynamics.Py(2:ayRows:end,ayRows+1:end);

dynamics.PcopxX = dynamics.Px(3:axRows:end,1:axRows);
dynamics.PcopxU = dynamics.Px(3:axRows:end,axRows+1:end);

dynamics.PcopyX = dynamics.Py(3:ayRows:end,1:ayRows);
dynamics.PcopyU = dynamics.Py(3:ayRows:end,ayRows+1:end);

dynamics.LdotxX = dynamics.Px(4:axRows:end,1:axRows);
dynamics.LdotxU = dynamics.Px(4:axRows:end,axRows+1:end);

dynamics.LdotyX = dynamics.Py(4:ayRows:end,1:ayRows);
dynamics.LdotyU = dynamics.Py(4:ayRows:end,ayRows+1:end);


dynamics.FxX = dynamics.Px(5:axRows:end,1:axRows);
dynamics.FxU = dynamics.Px(5:axRows:end,axRows+1:end);

dynamics.FyX = dynamics.Py(5:ayRows:end,1:ayRows);
dynamics.FyU = dynamics.Py(5:ayRows:end,ayRows+1:end);

[H,f] = computeCost(dynamics,gait, constants);

u0 = -H\f';

Aiq = [];
biq = [];

Aeq = [];
beq = [];
A = [];
b = [];

 [Acop,bcop,gait] = getCopConstraint(gait,dynamics,robot,constants);
 [Aiq,biq] = add_constraint(Aiq,biq,Acop,bcop);

 % [Acmp,bcmp] = getCmpConstraint(gait,dynamics,robot,constants);
 %[Aiq,biq] = add_constraint(Aiq,biq,Acmp,bcmp);
 
% Final COM Position Constraint %%%
AcomEndx = dynamics.PcomxU(end,:);
bcomEndx = gait.COPs{end}(1) - dynamics.PcomxX(end,:)*dynamics.initialConditionX;
AcomEndy = dynamics.PcomyU(end,:);
bcomEndy = gait.COPs{end}(2) - dynamics.PcomyX(end,:)*dynamics.initialConditionY;
AcomEnd = blkdiag(AcomEndx,AcomEndy);
bcomEnd = [bcomEndx;bcomEndy];
[Aeq,beq] = add_constraint(Aeq,beq,AcomEnd,bcomEnd);

%%% Final dcm Constraint %%%
AdcmEndx=dynamics.dcmxU(end,:);
bdcmEndx=gait.COPs{end}(1)-dynamics.dcmxX(end,:)*dynamics.initialConditionX;
AdcmEndy=dynamics.dcmyU(end,:);
bdcmEndy=gait.COPs{end}(2)-dynamics.dcmyX(end,:)*dynamics.initialConditionY;
AdcmEnd=blkdiag(AdcmEndx,AdcmEndy);
bdcmEnd=[bdcmEndx;bdcmEndy];
[Aeq,beq] = add_constraint(Aeq,beq,AdcmEnd,bdcmEnd);

%%% Final COP  Constraint %%%
AcopEndx=dynamics.PcopxU(end,:) ;
bcopEndx= gait.COPs{end}(1)-dynamics.PcopxX(end,:)*dynamics.initialConditionX;
AcopEndy=dynamics.PcopyU(end,:);
bcopEndy= gait.COPs{end}(2)-dynamics.PcopyX(end,:)*dynamics.initialConditionY;
AcopEnd=blkdiag(AcopEndx,AcopEndy);
bcopEnd=[bcopEndx;bcopEndy];

[Aeq,beq] = add_constraint(Aeq,beq,AcopEnd,bcopEnd);

%%note: the ldotxX is for y direction equation so use y initial conditions
%%% Final Ldot Constraint %%%
ALdotEndx = dynamics.LdotxU(end,:);
ALdotEndx =[ALdotEndx];
bLdotEndx = -dynamics.LdotxX(end,:)*dynamics.initialConditionY;
ALdotEndy = dynamics.LdotyU(end,:);
ALdotEndy =[ALdotEndy];
bLdotEndy = -dynamics.LdotyX(end,:)*dynamics.initialConditionX;
ALdotEnd  = blkdiag(ALdotEndy,ALdotEndx);
bLdotEnd  = [bLdotEndy;bLdotEndx];
[Aeq,beq] = add_constraint(Aeq,beq,ALdotEnd,bLdotEnd);

if(constants.runningInMatlab)
options = optimset('Algorithm','interior-point-convex','Display','final','TolFun',1e-6,'TolX',1e-6,'MaxIter',1000);
tic
x = quadprog(H,f,Aiq,biq,Aeq,beq,[],[],u0,options);
toc
else
options = optimset('Display','final','TolFun',1e-6,'TolX',1e-6,'MaxIter',1000);
[x, OBJ, INFO, LAMBDA] = qp(u0,H,f',Aeq,beq);
end

ux = x(1:bxCols*constants.N);
uy = x(byCols*constants.N+1:end);

px = ux(1:2:end);
py = uy(1:2:end);



Lddoty = ux(2:2:end);
Lddotx = uy(2:2:end);

gait.comX = [dynamics.initialConditionX(1);dynamics.PcomxX * dynamics.initialConditionX + dynamics.PcomxU * ux];

gait.comY = [dynamics.initialConditionY(1);dynamics.PcomyX * dynamics.initialConditionY + dynamics.PcomyU * uy];
gait.comZ = dynamics.z;
gait.com = [gait.comX,gait.comY,gait.comZ];

gait.dcmX = [dynamics.initialConditionX(2);dynamics.dcmxX * dynamics.initialConditionX + dynamics.dcmxU * ux];
gait.dcmY = [dynamics.initialConditionY(2);dynamics.dcmyX * dynamics.initialConditionY + dynamics.dcmyU * uy];


gait.copX = [dynamics.initialConditionX(3);dynamics.PcopxX * dynamics.initialConditionX + dynamics.PcopxU * ux];
gait.copY = [dynamics.initialConditionY(3);dynamics.PcopyX * dynamics.initialConditionY + dynamics.PcopyU * uy];

gait.comddotX=((constants.gravity)./(constants.initialCenterOfMassHeight))*(gait.comX-gait.copX);
gait.comddotY=((constants.gravity)./(constants.initialCenterOfMassHeight))*(gait.comY-gait.copY);

gait.LdotX = [dynamics.initialConditionY(4);dynamics.LdotyX*dynamics.initialConditionY + dynamics.LdotyU * uy];
gait.LdotY = [dynamics.initialConditionX(4);dynamics.LdotxX*dynamics.initialConditionX + dynamics.LdotxU * ux];

gait.cmpX = gait.copX + gait.LdotY./(constants.mass.*(dynamics.zddot + constants.gravity));
gait.cmpY = gait.copY - gait.LdotX./(constants.mass.*(dynamics.zddot + constants.gravity));

gait.comdotX=(gait.dcmX-gait.comX)/sqrt(0.7/9.8) ;
gait.comdotY=(gait.dcmY-gait.comY)/sqrt(0.7/9.8) ;

% gait.comY
% gait.comX 
% gait.comdotY
% gait.comdotX
gait.ICPX=[dynamics.initialConditionX(2);dynamics.dcmxX * dynamics.initialConditionX + dynamics.dcmxU * ux]; ;
gait.ICPY=[dynamics.initialConditionY(2);dynamics.dcmyX * dynamics.initialConditionY + dynamics.dcmyU * uy]; ;

x=[0;gait.comX];
y=[0;gait.comY];
z=[0;gait.comZ];


time=[0;gait.t'];
save('xcom.mat','x')
save('zcom.mat','z')
save('ycom.mat','y')

vx=[0;gait.comdotX];
vy=[0;gait.comdotY];
vz=zeros(78,1);

save('vxcom.mat','vx')
save('vzcom.mat','vz')
save('vycom.mat','vy')

ax=[0;gait.comddotX];
ay=[0;gait.comddotY];
az=zeros(78,1);

save('axcom.mat','ax')
save('azcom.mat','az')
save('aycom.mat','ay')


% gait.ICPX
% gait.ICPX
end
