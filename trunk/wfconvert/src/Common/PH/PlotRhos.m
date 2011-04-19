% From plane waves
PlaneRho = load('rho.dat'); 

% From radial equation
RadialRho = load('RadialRho.dat'); 

plot (PlaneRho(:,1), 1.133*PlaneRho(:,2), RadialRho(:,1), RadialRho(:,2), ...
      'r');

xlabel ('x (bohr)');
ylabel ('\rho(x)');
title ('Electron density comparison for Na atom');
legend ('From plane waves', 'From radial eq.');
