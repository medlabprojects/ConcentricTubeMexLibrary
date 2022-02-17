function kinematicsInfo = getKinematics(tubes, psiL, beta, NInterp)

%Run the forward kinematics
[s, pStar, qStar, Jh, stability, Ch] = ...
    forwardKinematics(tubes.tube1, tubes.tube2, tubes.tube3, tubes.tube4, psiL, beta);

%Interpolate the solution
[sInterp, gStarInterp] = interpolateBackbone(s,pStar,qStar,NInterp);

%Collect into a structure
kinematicsInfo = collectIntoStruct(Jh, Ch, stability, sInterp, gStarInterp);

end