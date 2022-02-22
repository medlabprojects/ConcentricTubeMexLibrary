function [p_tip,s,pStar,qStar,Jh,JStar, kin] = fwkin_TwoTube(robot,psiL,beta)

    kin = TwoTubeMexWithPsi(robot.tube1,robot.tube2,psiL,beta);
    
    % kinRet contains:
    % p_tip (3x1) - in base frame
    % q_tip (4x1) - in base frame
    % J_tip (6x4) - pose values in base frame wrt joint values in tip frame
    % s (Nx1)
    % p (Nx3)
    % q (Nx4)
    % J (Nx1 cell containing 6x4 matrices)

    s = kin.s;
    pret = kin.p;  % in tip frame
    qret = kin.q;  % in tip frame
    p_tip = kin.p_tip;

    N = length(s); % number of points returned from kinematics

    % These are what we want to fill (in base frame):
    pStar = zeros(N,3);
    qStar = zeros(N,4);

    % We'll use these to help fill out pStar and qStar values for each s:
    g = zeros(4,4); g(4,4) = 1.0;

    % Find the base frame s value:
    [~,zeroIndex] = min(abs(s));
    zeroIndex = zeroIndex(end); % in case more than one zero index is found

    % Find the base frame transformation with respect to the tip frame
    qZero = qret(zeroIndex,:);
    pZero = pret(zeroIndex,:);
    rot = quat2rotm(qZero);

    gZero = [rot, pZero'; % in tip frame
             0 0 0 1];

    % The base frame wrt the base frame will be identity:
    gStarZero = eye(4);

    % The tip frame expressed in the base frame:
    Rgstar = gStarZero(1:3,1:3);
    Pgstar = gStarZero(1:3,4);
    invg = [Rgstar' Rgstar'*Pgstar;
            0 0 0 1];

    gStarL = gStarZero*invg;
    R_tip = quat2rotm(kin.q_tip');

    % Run through the points and find them wrt base frame instead of tip:
    for i = 1:N
        g(1:3,1:3) = quat2rotm(qret(i,:));
        g(1:3,4) = pret(i,:)';
        gStar = robot.base*gStarL*g;

        % Now pull p & q out of the transformation:
        pStar(i,:) = gStar(1:3,4)';
        qStar(i,:) = rotm2quat(gStar(1:3,1:3));
    end

    % Find the hybrid Jacobian for tip motion:
    Jh = [R_tip zeros(3,3);
          zeros(3,3) R_tip] * kin.J_tip;

    % Find the hybrid Jacobian for each point along the backbone:
    JStar{N} = [];
    qstar = 0;
    for i = 1:N
        JStar{i} = [R_tip zeros(3,3);
        zeros(3,3) R_tip]*kin.J{i};
    end

end










