%% Animation:
% ====================================================

function AnimEOM ( t , z , rj , q , u )

% Data:

sk = 20; % step skip
p = 1e-1; % pause time
qu = [ q , u ];

[ n , m ] = size ( rj ); % number of links
ss = length ( t ); % simulation steps
sa = floor ( ( ss - 1 ) / sk ) + 2; % animation steps
n = n + 1;
rjt = zeros ( n , 4 , sa );

parfor i = 1 : ss
%     rjt(:,:,i) = subs (  rj , qu , z(i,:) );
    rjt(:,:,i) = [ 0 0 0 0 ; rjF( z(i,:) ) ];
end

rtip = zeros(1,3);
parfor i = 1 : ss
    rtip(i,:) = [rjt(end,1,i) rjt(end,2,i) rjt(end,3,i)+0.0075]; % corrected
end
% rjt(:,:,end) = subs (  rj , qu , z(end,:) );
rjt(:,:,end) = [ 0 0 0 0 ; rjF( z(end,:) ) ];
rjt = double ( rjt );

% tip position
figure;
plot(t,rtip(:,1),t,rtip(:,2),t,rtip(:,3)); title('sim tip pos');

% Windows set:
figure;
clf

window_xmin = min ( min ( rjt(:,1,:) ) );
window_xmax = max ( max ( rjt(:,1,:) ) );
window_ymin = min ( min ( rjt(:,2,:) ) );
window_ymax = max ( max ( rjt(:,2,:) ) );
window_zmin = min ( min ( rjt(:,3,:) ) );
window_zmax = max ( max ( rjt(:,3,:) ) );

ad = max( abs ( [ window_xmin window_xmax window_ymin window_ymax window_zmin window_zmax ] ) );

axis('equal')
axis ( [ window_xmin window_xmax window_ymin window_ymax window_zmin window_zmax ] ...
    + [ -ad ad -ad ad -ad ad ] );
xlabel('x [m]') , ylabel('y [m]') , zlabel('z [m]')

set(gcf,'Color',[1,1,1])


% Animation:

clr = { 'blue' , 'red' , 'green' , 'magenta' , 'cyan' , 'yellow' , 'black' };

for i = 1 : n - 1

	if rjt(i+1,4,1) == 0
    
		xs = [ rjt(i,1,1) , rjt(i+1,1,1) ];
		ys = [ rjt(i,2,1) , rjt(i+1,2,1) ];
		zs = [ rjt(i,3,1) , rjt(i+1,3,1) ];
    
	else
	
		ip = rjt(i+1,4,1);
		xs = [ rjt(ip,1,1) , rjt(i+1,1,1) ];
		ys = [ rjt(ip,2,1) , rjt(i+1,2,1) ];
		zs = [ rjt(ip,3,1) , rjt(i+1,3,1) ];	
	
        
	end
	
    if xs(1) - xs(2) ~= 0 || ys(1) - ys(2) ~= 0 || zs(1) - zs(2) ~= 0
      
        ln(i) = line ( xs , ys , zs , ...
           'linewidth' , 2 , 'erase' , 'xor' ,'color' , clr{i-floor(i/8)*7} );
    
    end

end

for i1 = 1 : sk : ss
    
    title( t(i1) );
    
    for i = 1 : n - 1
        
		if rjt(i+1,4,i1) == 0
    
			xs = [ rjt(i,1,i1) , rjt(i+1,1,i1) ];
			ys = [ rjt(i,2,i1) , rjt(i+1,2,i1) ];
			zs = [ rjt(i,3,i1) , rjt(i+1,3,i1) ];
    
		else
	
			ip = rjt(i+1,4,i1);
			xs = [ rjt(ip,1,i1) , rjt(i+1,1,i1) ];
			ys = [ rjt(ip,2,i1) , rjt(i+1,2,i1) ];
			zs = [ rjt(ip,3,i1) , rjt(i+1,3,i1) ];
            
		end
		
        if xs(1) - xs(2) ~= 0 || ys(1) - ys(2) ~= 0 || zs(1) - zs(2) ~= 0
            
            set( ln(i) , 'XData' , xs , 'YData' , ys , 'ZData' , zs );

        end
        
    end
            
%     if floor((i1-1)/10) == 0; hold on; end
    drawnow 
    pause( p );

end