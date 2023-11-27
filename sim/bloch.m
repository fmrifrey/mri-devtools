function dMdt = bloch(t,M,M0,B1,G,T1,T2,r)
    gamma = 4.258; % Gyromagnetic ratio (Hz/mG)
    
    % Compute effective magnetic field
    B = B1(t) + [0; 0; dot(G(t),r)];
    
    % Compute dM/dt using the bloch equation
    dMdt = cross(M,gamma*B) - ... % Precession
        1/T2*[M(1); M(2); 0] - ... % Transverse relaxation
        1/T1*[0; 0; M(3)-M0(3)]; % Longitudinal relaxation
end