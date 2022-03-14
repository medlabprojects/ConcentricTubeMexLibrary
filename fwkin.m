%% FWKIN
% INPUTS
%   robot: struct containing the n tube structs
%   q: array consisting of values [base rotation, base translation], where
%       rotation alpha or psi [1xn] and translation is beta [1xn]
% OUTPUTS
%   p_tip: [3x1]   position of tip in base frame
%   s:     [nx1]   arc length positions of sampled points
%   ps:    [3xn]   positions of backbone points at position s in base frame
%   qs:    [4xn]   quaternions of backbone points at position s in base frame
%   Jh:    [6x6]   Jacobian at tip
%   Js:    {6x6xn} Jacobian at backbone points
%   kin:           raw struct output of kinematics
% TubeN: Struct defining tube parameters with following elements
%       OD: Outer diameter
%       ID: Inner diameter
%       k:  Tube precurvature (const)
%       L:  Total tube length
%       Lt: Transmission length
%       E:  Young's Modulus
%       G:  Shear Modulus

function [p_tip,s,ps,qs,Jh,Js,kin] = fwkin(robot,q)
    n = length(q);
    nTubes = n/2;
    alphas = q(1:nTubes);
    betas = q(nTubes+1:end);

    % flip if not row vector
    if size(q,1)~= 1
        alphas = alphas';
        betas = betas';
    end
    
    switch nTubes
        case 2
            psis = find_psiTwoTube(robot, alphas, betas);
            kin = TwoTubeMexWithPsi(robot.tube1, robot.tube2, psis, betas);
        case 3
            psis = find_psiThreeTube(robot, alphas, betas);
            kin = ThreeTubeMexWithPsi(robot.tube1, robot.tube2, robot.tube3, psis, betas);
        case 4
            psis = find_psiFourTube(robot, alphas, betas);
            kin = FourTubeMexWithPsi(robot.tube1, robot.tube2, robot.tube3, robot.tube4, psis, betas);
        otherwise
            disp('wrong number of tubes')
    end

    s = kin.s;
    pret = kin.p;  % in tip frame
    qret = kin.q;  % in tip frame
    p_tip = kin.p_tip;

    N = length(s); % number of points returned from kinematics

    % These are what we want to fill (in base frame):
    ps = zeros(N,3);
    qs = zeros(N,4);

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
        ps(i,:) = gStar(1:3,4)';
        qs(i,:) = rotm2quat(gStar(1:3,1:3));
    end

    % Find the hybrid Jacobian for tip motion:
    Jh = [R_tip zeros(3,3);
          zeros(3,3) R_tip] * kin.J_tip;
    kin.Jh = Jh;

    % Find the hybrid Jacobian for each point along the backbone:
    Js{N} = [];
    qstar = 0;
    for i = 1:N
        Js{i} = [R_tip zeros(3,3);
        zeros(3,3) R_tip]*kin.J{i};
    end
end
