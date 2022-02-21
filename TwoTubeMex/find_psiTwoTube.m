function psiL = find_psiTwoTube(robot,alpha,beta)
    options = optimoptions('fsolve');
    options.Display = 'off';    
            
    psiL = fsolve(@(psiL) base_angle_err2(psiL,alpha,beta,robot),alpha,options);
    end
    
function e = base_angle_err2(psiL,alpha,beta,robot)
    kinRet = TwoTubeMexWithPsi(robot.tube1,robot.tube2,psiL,beta);
    psibase = kinRet.psi(end,1:2);
    e = (psibase - alpha);
end