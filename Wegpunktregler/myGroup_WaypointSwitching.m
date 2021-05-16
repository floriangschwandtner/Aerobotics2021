function [CurrentWP_ID, PreviousWP] = myGroup_WaypointSwitching(r, VK, WayPoints, CurrentWP_ID, PreviousWP)
% This function contains the algorithm for switching the Waypoints.
%
% Inputs:
%   - r             The current position of the Aircraft -> [North, East, Down]
%   - VK            Ground Speed of the Aircraft
%   - WayPoints     A matrix containing all WayPoints (one Waypoint(WP) per row: [North, East, Down]) 
%   - Current_WP_ID The ID/row number of the current WP
%   - PreviousWP    The last WP [North, East, Down]
%
% Outputs:
%   - Current_WP_ID The ID/row number of the current WP - Must only be changed if the switch condition is fullfilled 
%   - PreviousWP    The last WP [North, East, Down] - Must only be changed if the switch condition is fullfilled

    % Distanz der momentanen Position zum nächsten Wegpunkt als Kugel. Kann
    % problematisch werden, wenn Höhenregler schlecht. Dann wäre ein
    % Zylinder sinnvoller.
    Dist = sqrt((r(1)-WayPoints(1,CurrentWP_ID))^2+(r(2)-WayPoints(2,CurrentWP_ID))^2+(r(3)-WayPoints(3,CurrentWP_ID))^2);

    basic = true;
    
    if(basic)
       minDist = 40.0;
    else
        %Abhängigkeit des minimalen Abstands von Bahnparametern, kann durch
        %"basic" flag umgeschaltet werden
        %-PC    Vektor von vorherigem Punkt zum aktuellen
        PC = WayPoints(:,CurrentWP_ID)-PreviousWP;
        
        %-CN    Vektor von aktuellem zum nächsten Punkt
        if(CurrentWP_ID<size(WayPoints,2))
            CN = WayPoints(:,CurrentWP_ID+1)-WayPoints(:,CurrentWP_ID);
        else
            CN = WayPoints(:,1)-WayPoints(:,CurrentWP_ID);
        end
        
        %Winkel zwischen den Vektoren
        ang_PN = atan2(norm(cross(PC,CN)),dot(PC,CN));
        
        %Berechnung der minimalen Distanz, Parameter geschätzt. Hier noch
        %Verbesserungsbedarf!
        if(CN(3)>0)
            minDist = 3*VK+abs((ang_PN-pi/2))*55-0.75*CN(3)
        else
            minDist = 3*VK+abs((ang_PN-pi/2))*55
        end
    end
    
    %Umschalten der Wegpunkte, wenn Distanz klein genug.
    if(Dist<minDist)
        
        PreviousWP = WayPoints(:,CurrentWP_ID);
        CurrentWP_ID = CurrentWP_ID+1;
        
        %Springt zum Anfang der Liste, ermöglicht das Kreisen des
        %Flugzeugs, verhindert Absturz der Simulation.
        if(CurrentWP_ID > size(WayPoints,2))
            CurrentWP_ID = 1;
        end
    end
end

