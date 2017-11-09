%% EOM Numerical Integrator:
% ====================================================
% change this for different scenarios

function [ t , z , tfinal ] = SimEoM_VC_D_v0_1( qf , uf , lambdaf , sf , qf0 , par )
t0 = par.t0 ; dt = par.dt ; stepT = par.stepT ;

[ ~ , nq ] = size( qf ) ;
[ ~ , nlambda ] = size( lambdaf ) ;

eps = 0 ; % for singularity prevention
z0 = [ double( qf0 ) eps*ones( 1 , nlambda ) eps*ones( 1 , nq ) eps*ones( 1 , nlambda ) ] ; % initial condition

% Standard ODE solver:
par.nq = nq ; par.nlambda = nlambda ;
options = odeset () ;  %,'abstol',1*1e-6,'reltol',1*1e-6) ;
tspan = t0 : stepT: t0 + dt ;
[ t , z , tfinal ] = ode113( @EOM , tspan , z0 , options , par ) ;


function dz = EOM ( t , z , par )
if mod( t , 1 ) < 1e-6 ; t % report every second
end

lambdal = 1 ;
% lambdal = z(1) / par.l_b ;

% pressure parameters
d_p = par.d_p / sqrt( lambdal ) ; a_p = par.a_p ; n_seg = par.n_seg ;

p_i = [ par.p_rho( floor( t * 1/par.stepT + 1 ) , 1 ) ...
    par.p_rho( floor( t * 1/par.stepT + 1 ) , 2 ) ...
    par.p_rho( floor( t * 1/par.stepT + 1 ) , 3 ) ] - [ 1 1 1 ] * 1e5 ;
fj1 = [ ( p_i(1) + p_i(2) + p_i(3) )*a_p + ...
            -( p_i(1) + p_i(2) * cos( 2 * pi / 3 ) + p_i(3) * cos( - 2 * pi / 3  ) )*a_p*d_p/2*z(2) + ...
            -( p_i(2) * sin( 2 * pi / 3 ) + p_i(3) * sin( - 2 * pi / 3  ) )*a_p*d_p/2*z(3) ;
        -( p_i(1) + p_i(2) * cos( 2 * pi / 3 ) + p_i(3) * cos( - 2 * pi / 3  ) )*a_p*d_p/2*z(1) ;
        -( p_i(2) * sin( 2 * pi / 3 ) + p_i(3) * sin( - 2 * pi / 3  ) )*a_p*d_p/2*z(1) ] ;
% n_seg = 1 ;
fj = [] ; for i = 1 : n_seg ; fj = [ fj ; fj1 ] ; end % multi segment case

fex = [ 0 0 0 ] ;

u = z( par.nq + par.nlambda + 1 : end ) ;

% for rigid elements
M = double( MF( lambdal , z.' ) ) ;
T = double( TF( lambdal , z.' ) ) ;
D = double( DdF( lambdal , z.' ) ) ;
fg = double( fgF( lambdal , z.' ) ) ;
fj = fj + double( fjF( lambdal , z.' ) ) ;
Tex = double( TefF( lambdal , z.' ) ) ;

TMTr = T.' * M * T ;
Tdr = T.' * ( fg - M * D * u ) + fj + Tex.' * fex.' ;

% for non-rigid elements
TMTc = 0 ; Tdc = 0 ;
intlim = double( intlimF( lambdal , z.' ) ) ;
if par.nC ~= 0
%     for nC = 1 : par.nC
        [ tmp1 , tmp2 ] = crvint( lambdal , z , u , intlim(1) ) ;
        TMTc = TMTc + tmp1 ; Tdc = Tdc + tmp2 ;
%     end
end

% EOM
TMT = TMTr + TMTc ; Td = Tdr + Tdc ;
dzt = TMT \ Td ;
% dzt = pinv( TMT ) * Td ;
dz = [ u ; dzt ];


function [ TMTc , Tdc ] = crvint( var , z , u , intlim )
% trapezoidal numerical integration

n = 10 ; % integration step
TMTc = 0 ; Tdc = 0 ;
for i = 0 : 1/n : 1 
    
    s = i * intlim ;
    
    MC = double( MCF( var , z.' , s ) ) ;
    TC = double( TCF( var , z.' , s ) ) ;
    DC = double( DdCF( var , z.' , s ) ) ;
    fgC = double( fgCF( var , z.' , s ) ) ;
    
    if i == 0 || i == 1
        TMTc = TMTc + TC.' * MC * TC / 2 ;
        Tdc = Tdc + TC.' * ( fgC - MC * DC * u ) / 2 ;
    else
        TMTc = TMTc + TC.' * MC * TC ;
        Tdc = Tdc + TC.' * ( fgC - MC * DC * u ) ;
    end
    
end

TMTc = TMTc * intlim / n ; Tdc = Tdc * intlim / n ;

 
