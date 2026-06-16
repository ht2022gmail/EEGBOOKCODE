function drdt = excitatoryinhinitorynetwork(t, r, params)

drdt = zeros(2, 1);

WEE = params.WEE;
WIE = params.WIE;
WEI = params.WEI;
WII = params.WII;
IE = params.IE;
II = params.II;
tauE = params.tauE;
tauI = params.tauI;

drdt(1) = (-r(1) + max(WEE*r(1)+WEI*r(2)+IE,0))/tauE;
drdt(2) = (-r(2) + max(WIE*r(1)+WII*r(2)+II,0))/tauI;