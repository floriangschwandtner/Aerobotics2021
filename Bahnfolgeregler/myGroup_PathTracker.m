function [gammaCmd, chiCmd] = myGroup_PathTracker(r, VK, WayPoints, CurrentWP_ID, PreviousWP)
% This function contains the algorithm for the path-tracker ("Bahnregler").
%
% Inputs:
%   - r             The current position of the Aircraft -> [North, East, Down]
%   - VK            Ground Speed of the Aircraft
%   - WayPoints     A matrix containing all WayPoints (one Waypoint(WP) per row: [North, East, Down]) 
%   - CurrentWP_ID  The ID/row number of the current WP
%   - PreviousWP    The last WP [North, East, Down]
%
% Outputs:
%   - gammaCmd      The commanded flight-path angle
%   - chiCmd        The commanded pitch angle
% 
% Empfohlene Vorgehensweise:
%   1.) Berechnen Sie die Sollbahnneigung (gamma0) und den Sollbahnazimuth (chi0)
%   2.) Bauen Sie die Transfomationsmatrix TrE auf
%   3.) Transfomieren Sie die aktuelle Flugzeugposition in Sollbahnkoordinaten und speichern Sie das Ergebnis als rrFZ
%   4.) Berechnen Sie den kommandierten Bahnwinkel (gammaCmd) und das kommandierte Heading (chiCmd)
%   5.) Limitieren Sie gammaSoll und chiCmd auf sinnvolle Wertebereiche

gammaCmd = single(0);
chiCmd   = single(0);

end