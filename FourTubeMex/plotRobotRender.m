function h = plotRobotRender(tubes, kinematicsInfo, beta, NTheta)

L = [tubes.tube1.L tubes.tube2.L tubes.tube3.L]';
tubeEndIndices = getTubeEndIndices(beta, L, kinematicsInfo.s);
tubeBackbone = getTubeBackbone(tubes, kinematicsInfo, tubeEndIndices, NTheta);

h1 = surf(tubeBackbone.tube1.x, tubeBackbone.tube1.y, tubeBackbone.tube1.z); hold on;
h2 = surf(tubeBackbone.tube2.x, tubeBackbone.tube2.y, tubeBackbone.tube2.z); hold on;
h3 = surf(tubeBackbone.tube3.x, tubeBackbone.tube3.y, tubeBackbone.tube3.z); hold on;
h4 = surf(tubeBackbone.TwoToOne.x, tubeBackbone.TwoToOne.y, tubeBackbone.TwoToOne.z); hold on;
h5 = surf(tubeBackbone.ThreeToTwo.x, tubeBackbone.ThreeToTwo.y, tubeBackbone.ThreeToTwo.z); hold on;
h6 = surf(tubeBackbone.OneToGrip.x, tubeBackbone.OneToGrip.y, tubeBackbone.OneToGrip.z); hold on;


h = [h1 h2 h3 h4 h5 h6];
set(h, 'FaceColor', rgb('Silver')); %type: "rgb chart" in command window to see choices
set(h, 'EdgeAlpha', 0);
