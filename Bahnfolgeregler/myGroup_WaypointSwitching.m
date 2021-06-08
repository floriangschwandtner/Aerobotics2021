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
    
    %Rollwinkel

    beta  = 30;
    Kugel = true;
    basic = false;
    
    if(Kugel)
        Dist = sqrt((r(1)-WayPoints(1,CurrentWP_ID))^2+(r(2)-WayPoints(2,CurrentWP_ID))^2+(r(3)-WayPoints(3,CurrentWP_ID))^2);
    else
        Dist = sqrt((r(1)-WayPoints(1,CurrentWP_ID))^2+(r(2)-WayPoints(2,CurrentWP_ID))^2);
    end
        
    if(basic)
       minDist = 40.0; %Zylinder oder Kugel als Umschaltbedingung
    else        
        %Abhängigkeit des minimalen Abstands von Bahnparametern, kann durch
        %"basic" flag umgeschaltet werden
        % Kurvenradius als Umschaltbedingung
        
        if CurrentWP_ID + 1 > size(WayPoints, 2)
            nextWP_ID = 1;
        else
            nextWP_ID = CurrentWP_ID + 1;
        end
        CurrentWP = WayPoints(:,CurrentWP_ID);
        NextWP    = WayPoints(:,nextWP_ID);
        
        %-PC    Vektor von vorherigem Punkt zum aktuellen
        
        PC     = CurrentWP-PreviousWP;
        
        %-CN    Vektor von aktuellem zum nächsten Punkt
        
        CN     = NextWP-CurrentWP;
        
        %Winkel zwischen den Vektoren
        
        cos_path_angle = PC(1:2)'*CN(1:2)/(norm(PC(1:2))*norm(CN(1:2)));
        cos_path_angle = min(max(cos_path_angle,-1),1);
        path_angle     = acos(cos_path_angle);
        
        %Limitierung der Winkel ziwschen den Wegpunkten auf sinnvolle Werte
        
        if path_angle >= 150*pi/180
            path_angle = single(150*pi/180);
        elseif path_angle <= 20*pi/180
            path_angle = single(20*pi/180);
        end
        
        % Berechnung des Kurvenradius und des Umschaltradius
        
        phiTurn = beta*pi/180;

        KurvenRadius = VK^2 / (tan(phiTurn)*9.81);
        switchRadius = KurvenRadius*tan(path_angle/2);
        
        minDist      = switchRadius;
    end
    
    %Umschalten der Wegpunkte, wenn Distanz klein genug.
    if(Dist<minDist)
        
        PreviousWP   = WayPoints(:,CurrentWP_ID);
        
        %Springt zum Anfang der Liste, ermöglicht das Kreisen des
        %Flugzeugs, verhindert Absturz der Simulation.
        
        if CurrentWP_ID >= size(WayPoints,2)
            CurrentWP_ID = 1;
        else
            CurrentWP_ID = CurrentWP_ID + 1;
        end
    end
end



