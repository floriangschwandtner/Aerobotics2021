function [gammaCmd, chiCmd] = Fabrizio_PathTracker(r, VK, WayPoints, CurrentWP_ID, PreviousWP)
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

%% Parameter
ky =  4.4928;
kz =  2.0406;
gammaMax = pi/4;
chiMax = pi/2;

%% Berechnung der Sollbahnwerte
% Sollbahnneigung
gamma0 = atan2(Waypoints(CurrentWP_ID,3)-PreviousWP(3),Waypoints(CurrentWP_ID,2)-PreviousWP(2));
% Sollbahnazimuth
chi0 = atan2(Waypoints(CurrentWP_ID,2)-PreviousWP(2),Waypoints(CurrentWP_ID,1)-PreviousWP(1));

%% Transformation der aktuellen Flugzeugposition in Sollbahnkoordinaten
% Transformationsmatrix
TrE = [cos(gamma0)*cos(chi0) cos(gamma0)*sin(chi0) -sin(gamma0);...
        -sin(chi0) cos(chi0) 0;...
        sin(gamma0)*cos(chi0) sin(gamma0)*sin(chi0) cos(gamma0)];

rrFZ = Tre*(r-PreviousWP);

%% Berechnung der Kommandos
% Pseudosteuerungen
uH = -ky*rrFZ(2);
uV = -kz*rrFZ(3);
% Begrenzung der Pseudosteuerungen
uH = uH;
uV = uV;
% Heading
chiCmd = chi0+asin(uH/(cos(gamma0*VK)));
% Bahnwinkel
a = -cos(gamma0);
b = sin(gamma0)*cos(chiCmd-chi0);
c = uV/VK;
delta = atan2(b,-a);
gammaCmd = asin(-c/sqrt(a^2+b^2)) + delta;

%% Limitierung der Kommandos auf sinnvolle Bereiche
chiCmd = min(max(chiCmd,-chiMax),chiMax);
gammaCmd = min(max(gammaCmd,-gammaMax),gammaMax);
end