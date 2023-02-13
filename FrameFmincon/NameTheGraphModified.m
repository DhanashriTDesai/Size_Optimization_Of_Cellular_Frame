function NameTheGraphModified(xLabel,yLabel,Title,nLegends,Leg1,Leg2,Leg3,Leg4,Leg5,Leg6,Leg7,Leg8,LegLocation)

xlabel(xLabel,'Interpreter','latex','fontWeight','bold'); 
ylabel(yLabel,'Interpreter','latex','fontWeight','bold');
title(Title,'FontSize', 32,'Interpreter','latex');
set(gca,'FontSize',30); set(gca,'FontWeight','bold');
h=gca; h.LineWidth=3;

switch nLegends
    case 2
        legend(Leg1,Leg2,'FontSize', 20,'Location',LegLocation,'Interpreter','latex');
    case 3
        legend(Leg1,Leg2,Leg3,'FontSize', 20,'Location',LegLocation,'Interpreter','latex');
    case 4
        legend(Leg1,Leg2,Leg3,Leg4,'FontSize', 20,'Location',LegLocation,'Interpreter','latex');
    case 6
        legend(Leg1,Leg2,Leg3,Leg4,Leg5,Leg6,'FontSize',20,'Location',LegLocation,'Interpreter','latex');
    case 7
        legend(Leg1,Leg2,Leg3,Leg4,Leg5,Leg6,Leg7,'FontSize',20,'Location',LegLocation,'Interpreter','latex');
    case 8
        legend(Leg1,Leg2,Leg3,Leg4,Leg5,Leg6,Leg7,Leg8,'FontSize',20,'Location',LegLocation,'Interpreter','latex');
end

end