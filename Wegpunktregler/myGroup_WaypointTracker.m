function [hCmd, chiCmd] = myGroup_WaypointTracker(r, WayPoints, CurrentWP_ID)
% This function contains the algorithm for the waypoint-tracker ("Wegpunktregler").
%
% Inputs:
%   - r             The current position of the Aircraft -> [North, East, Down]
%   - V_K           Ground Speed of the Aircraft
%   - WayPoints     A matrix containing all WayPoints (one Waypoint(WP) per row: [North, East, Down]) 
%   - Current_WP_ID The ID/row number of the current WP
%   - PreviousWP    The last WP [North, East, Down]
%
% Outputs:
%   - hCmd         The computed "height"-command (see note above)
%   - chiCmd       The commanded flight-path angle
% 
% Empfohlene Vorgehensweise:
%   1.) Berechnen Sie das Heading von der aktuellen Flugzeugposition zum nächsten Wegpunkt
%   2.) Übergeben Sie das Heading an chiSoll und die Höhe des nächsten Zielpunkts an hSoll (Achtung Vorzeichen!)
    
    % Kommando nach Formel in Folien
    chiCmd = atan2(WayPoints(2,CurrentWP_ID)-r(2),WayPoints(1,CurrentWP_ID)-r(1));
    
    % Höhenkommando, vorzeichen umgekehrt.
    hCmd = -single(WayPoints(3,CurrentWP_ID));
    
end