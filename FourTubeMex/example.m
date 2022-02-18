clf;

tubes = defineTubes(); %all geometry
frontPlate = setFrontPlateSpecs(); %size, thickness, color

theta1 = linspace(0,360,60)'; %rotational vector for plotting
t = linspace(0,4,60)'; %time vector for plotting

NInterp = 50; %number of arc length points
NTheta = 10; %number of points on circle around each arc length

for i=1:length(theta1)
    
    psiL = [0  *deg2rad(theta1(i)),...
            1  *deg2rad(theta1(i)),...
            0 *deg2rad(theta1(i)),...
            0 *deg2rad(theta1(i))]; %rotational joints
        
    beta = [-95e-3+ (10e-3)*sin(0.2*pi*t(i)),...
            -60e-3+ (5e-3)*sin(0.2*pi*t(i)),...
            -25e-3+ (5e-3)*cos(0.2*pi*t(i)),...
             -10e-3+ (5e-3)*cos(0.2*pi*t(i))]; %translational joints
    
    if (i~=1)    
        delete(allHandles);
    end

    kinematicsInfo = getKinematics(tubes, psiL, beta, NInterp);
    allHandles = plotRobotRender(tubes, kinematicsInfo, beta, NTheta);
    
    if (i==1)
        setPlotEnvironment(frontPlate);
    end
    
    lighting phong;
    drawnow;
end
