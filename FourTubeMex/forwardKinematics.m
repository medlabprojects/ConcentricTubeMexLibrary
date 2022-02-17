function [s, pStar, qStar, Jh, stability, Ch] = forwardKinematics(tube1, tube2, tube3, tube4, psiL, beta)
kinRet = FourTubeMex(tube1, tube2, tube3, tube4, psiL, beta); 
%p_tip: Position of the tip in the s = 0 frame
%q_tip: Quaternion representation of tip in the s=0 global frame
%J_tip: The body frame jacobian at the tip of the robot
%C_tip: The body frame compliance matrix at the tip of the robot
%stability: The stability metric
%s: Arc lengths - backwards integrated
%p: - backwards (tip position is [0 0 0]')
%q: - backwards (tip rotation is identity)

N = length(kinRet.s);
pStar = zeros(N,3);
qStar = zeros(N,4);
g = zeros(4,4); g(4,4) = 1.0;
s = kinRet.s;

zeroIndices = (s == 0);
k = find(zeroIndices);

if (isempty(k))
     error('Cannot find the zero position');
else
     zeroIndex = k(1);
end

qZero = kinRet.q(zeroIndex,:);
pZero = kinRet.p(zeroIndex,:);
gZero = [quat2rotm(qZero) pZero';
         zeros(1,3) 1];
     
gStarZero = [eye(3) zeros(3,1); zeros(1,3) 1];
gStarL = gStarZero * inverse_g(gZero);
R_tip = quat2rotm(kinRet.q_tip');

for i=1:N
    g(1:3,1:3) = quat2rotm(kinRet.q(i,:));
    g(1:3,4) = kinRet.p(i,:)';
    gStar = gStarL * g;
    qStar(i,:) = rotm2quat(gStar(1:3,1:3));
    pStar(i,:) = gStar(1:3,4)';
end

Jh = [R_tip zeros(3,3); 
      zeros(3,3) R_tip] * kinRet.J_tip;
Ch = [R_tip zeros(3,3); 
      zeros(3,3) R_tip] * kinRet.C_tip;
stability = kinRet.stability;