 function[Bx_temp,By_temp,Bz_temp] = Magnetic_Field(Coil_Input,x,y,z)
 
 %% Magnetic Field in Ring

% Magnetic Field Formulas : Solving by Biot-Savart Law will yield Bx,By and
% Bz upto an integral equation similar to elliptic integrals. The numerically solved
% and simplified form of those integrals can be found in the paper "Simple 
% Analytic Expressions for the Magnetic Field of a Circular Current Loop" by James Simpson et al. 
% Link : https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010038494.pdf

% Formula in above paper is given assuming ring centre is at origin. Shift
% origin and compute input variables for B for both rings separately.

if (x==0)
    x = 1e-7;
end

if (y==0)
    y = 1e-7;
end

if (z==0)
    z = 1e-7;
end
% B due to Ring 1 ( Centre at (0,0,-d12) => (x,y,z) -> (x,y,z + d12) )
z0 = -1 * Coil_Input(2,1);
C = 4e-7*Coil_Input(3,1); % mu*I1/pi
rho_sq = x^2 + y^2;
rho = sqrt(rho_sq);

r_sq = x^2 + y^2 + (z - z0)^2;
alpha_sq = Coil_Input(1,1)^2 + r_sq - (2*Coil_Input(1,1)*rho);
beta_sq  = Coil_Input(1,1)^2 + r_sq + (2*Coil_Input(1,1)*rho);
beta = sqrt(beta_sq);
              
k_sq = 1 - (alpha_sq/beta_sq);
[K,E] = ellipke(k_sq);
constant = C/(2*alpha_sq*beta);
              
Bz1 = constant * (((Coil_Input(1,1)^2 - r_sq)*E) + (alpha_sq*K));
Bx1 = constant * ((x*(z - z0))/(rho_sq)) * (((Coil_Input(1,1)^2 + r_sq)*E) - (alpha_sq*K));
By1 = constant * ((y*(z - z0))/(rho_sq)) * (((Coil_Input(1,1)^2 + r_sq)*E) - (alpha_sq*K));

% B due to Ring 2 ( Centre at (0,0,0) => (x,y,z) -> (x,y,z) )

z0 = Coil_Input(2,1);
C = 4e-7*Coil_Input(3,2); % mu*I2/pi
rho_sq = x^2 + y^2;
rho = sqrt(rho_sq);

r_sq = x^2 + y^2 + (z - z0)^2;
alpha_sq = Coil_Input(1,2)^2 + r_sq - (2*Coil_Input(1,2)*rho);
beta_sq  = Coil_Input(1,2)^2 + r_sq + (2*Coil_Input(1,2)*rho);
beta = sqrt(beta_sq);
              
k_sq = 1 - (alpha_sq/beta_sq);
[K,E] = ellipke(k_sq);
constant = C/(2*alpha_sq*beta);
              
Bz2 = constant * (((Coil_Input(1,2)^2 - r_sq)*E) + (alpha_sq*K));
Bx2 = constant * ((x*(z - z0))/(rho_sq)) * (((Coil_Input(1,2)^2 + r_sq)*E) - (alpha_sq*K));
By2 = constant * ((y*(z - z0))/(rho_sq)) * (((Coil_Input(1,2)^2 + r_sq)*E) - (alpha_sq*K));

% B due to Ring 3 ( Centre at (0,0,+d23) => (x,y,z) -> (x,y,z - d23) )

z0 = 1 * Coil_Input(2,3);
C = 4e-7*Coil_Input(3,3); % mu*I3/pi
rho_sq = x^2 + y^2;
rho = sqrt(rho_sq);

r_sq = x^2 + y^2 + (z - z0)^2;
alpha_sq = Coil_Input(1,3)^2 + r_sq - (2*Coil_Input(1,3)*rho);
beta_sq  = Coil_Input(1,3)^2 + r_sq + (2*Coil_Input(1,3)*rho);
beta = sqrt(beta_sq);
              
k_sq = 1 - (alpha_sq/beta_sq);
[K,E] = ellipke(k_sq);
constant = C/(2*alpha_sq*beta);
              
Bz3 = constant * (((Coil_Input(1,3)^2 - r_sq)*E) + (alpha_sq*K));
Bx3 = constant * ((x*(z - z0))/(rho_sq)) * (((Coil_Input(1,3)^2 + r_sq)*E) - (alpha_sq*K));
By3 = constant * ((y*(z - z0))/(rho_sq)) * (((Coil_Input(1,3)^2 + r_sq)*E) - (alpha_sq*K));

% B due to Ring 4 ( Centre at (0,0,+d23+d34) => (x,y,z) -> (x,y,z - d23 - d34) )

z0 = 1 * (Coil_Input(2,3)+Coil_Input(2,4));
C = 4e-7*Coil_Input(3,4); % mu*I3/pi
rho_sq = x^2 + y^2;
rho = sqrt(rho_sq);

r_sq = x^2 + y^2 + (z - z0)^2;
alpha_sq = Coil_Input(1,4)^2 + r_sq - (2*Coil_Input(1,4)*rho);
beta_sq  = Coil_Input(1,4)^2 + r_sq + (2*Coil_Input(1,4)*rho);
beta = sqrt(beta_sq);
              
k_sq = 1 - (alpha_sq/beta_sq);
[K,E] = ellipke(k_sq);
constant = C/(2*alpha_sq*beta);
              
Bz4 = constant * (((Coil_Input(1,4)^2 - r_sq)*E) + (alpha_sq*K));
Bx4 = constant * ((x*(z - z0))/(rho_sq)) * (((Coil_Input(1,4)^2 + r_sq)*E) - (alpha_sq*K));
By4 = constant * ((y*(z - z0))/(rho_sq)) * (((Coil_Input(1,4)^2 + r_sq)*E) - (alpha_sq*K));

Bx_temp = Bx1+Bx2+Bx3+Bx4; % in [T]
By_temp = By1+By2+By3+By4; % in [T]
Bz_temp = Bz1+Bz2+Bz3+Bz4; % in [T]
 
 end

