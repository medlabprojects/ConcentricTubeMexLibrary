function psiL = find_psiThreeTube(robot,alpha,beta)
    options = optimoptions('fsolve');
    options.Display = 'off';    
            
    psiL = fsolve(@(psiL) base_angle_err3(psiL,alpha,beta,robot),alpha,options);
    end
    
function e = base_angle_err3(psiL,alpha,beta,robot)
    kinRet = ThreeTubeMexWithPsi(robot.tube1,robot.tube2,robot.tube3,psiL,beta);
    psibase = kinRet.psi(end,1:3);
    e = (psibase - alpha);
end