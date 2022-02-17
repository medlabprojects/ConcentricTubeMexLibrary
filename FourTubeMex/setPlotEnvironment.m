function setPlotEnvironment(frontPlate)

        axis equal;
        xlim([-.040 .040]);
        ylim([-.040 .040]);
        zlim([-frontPlate.thickness .080]);

        hLight1 = camlight('headlight');
        set(hLight1,'Color',rgb('lightYellow'));
        set(hLight1,'Position',[-.010 -.010 .020]);

        hLight2 = camlight('headlight');
        set(hLight2,'Color',rgb('lightYellow'));
        set(hLight2,'Position',[-.010 -.010 .050]);
        set(gca,'visible','off');
        
        plotFrontPlate(frontPlate.xlims, frontPlate.ylims, frontPlate.thickness, frontPlate.color);
end