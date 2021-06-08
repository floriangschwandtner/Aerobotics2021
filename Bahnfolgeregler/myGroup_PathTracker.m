function [gammaCmd, chiCmd] = myGroup_PathTracker(r, VK, WayPoints, CurrentWP_ID, PreviousWP)
% This function contains the algorithm for the path-tracker ("Bahnregler").
%
% Inputs:
%   - r             The current position of the Aircraft -> [North, East, Down]
%   - VK            Ground Speed of the Aircraft
%   - WayPoints     A matrix containing all WayPoints (one Waypoint(WP) per row: [North, East, Down]) 
%   - Current_WP_ID The ID/row number of the current WP
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

kH          = 1.4966/4;
kV          = 1.8195/4;
gammaMax    = 10/180*pi;

%   1.) %%%%%%%%%

gamma0   = -atan2(WayPoints(3,CurrentWP_ID)-PreviousWP(3),sqrt((WayPoints(2,CurrentWP_ID)-PreviousWP(2))^2+(WayPoints(1,CurrentWP_ID)-PreviousWP(1))^2));
chi0     = atan2(WayPoints(2,CurrentWP_ID)-PreviousWP(2),WayPoints(1,CurrentWP_ID)-PreviousWP(1));

%   2.) %%%%%%%%%

TrE = [cos(gamma0)*cos(chi0) cos(gamma0)*sin(chi0) -sin(gamma0);...
-sin(chi0) cos(chi0) 0;...
sin(gamma0)*cos(chi0) sin(gamma0)*sin(chi0) cos(gamma0)];

%   3.) %%%%%%%%%

rrFZ = TrE *(r-PreviousWP);

%   4.) %%%%%%%%%

% Lateral Inversion
uH       = -kH * rrFZ(2);
if abs(uH) > VK * cos(gamma0)
    uH = sign(-kH * rrFZ(2))*VK*cos(gamma0);
end

chiCmd   = chi0+asin(uH/(VK*cos(gamma0)));

% Vertical Inversion
a        = -cos(gamma0);
b        = sin(gamma0)*cos(chiCmd-chi0);

uV       = -kV * rrFZ(3);
if abs(uV) > VK*sqrt(a^2+b^2)
    uV = sign(uV) * VK*sqrt(a^2+b^2);
end

gammaCmd = atan2(b,-a)+asin(-uV/(VK*sqrt(a^2+b^2)));

%   5.) %%%%%%%%%

gammaCmd = max(-gammaMax, gammaCmd);
gammaCmd = min(+gammaMax, gammaCmd);

end