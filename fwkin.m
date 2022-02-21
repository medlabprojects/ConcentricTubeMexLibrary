%% FWKIN
% INPUTS
%   robot: struct containing the n tube structs
%   q: array consisting of values [base rotation, base translation], where
%       rotation alpha or psi [1xn] and translation is beta [1xn]
%   optional
% OUTPUTS
%   p_tip:
%   q_tip:
%   s:
%   p:
%   Jh:
%   J:
% TubeN: Struct defining tube parameters with following elements
%       OD: Outer diameter
%       ID: Inner diameter
%       k:  Tube precurvature (const)
%       L:  Total tube length
%       Lt: Transmission length
%       E:  Young's Modulus
%       G:  Shear Modulus

function [p_tip,s,p,q,Jh,J,kin] = fwkin(robot,q)
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
            [p_tip,s,p,q,Jh,J,kin] = fwkin_TwoTube(robot, psis, betas);
        case 3
            psis = find_psiThreeTube(robot, alphas, betas);
            [p_tip,s,p,q,Jh,J,kin] = fwkin_ThreeTube(robot, psis, betas);
        case 4
            psis = find_psiFourTube(robot, alphas, betas);
            [p_tip,s,p,q,Jh,J,kin] = fwkin_FourTube(robot, psis, betas);
        otherwise
            disp('wrong number of tubes')
    end
end
