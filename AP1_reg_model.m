function dydt = AP1_reg_model( t, y0, kinaseSignals )
% Calculates reaction rates vector at time t
% 
% Input:
%	t  - time
%	y0 - initial condition
%	kinaseSignals
%
% output:
%	dydt - change in activity levels
%

switch nargin
    case 2
	FRK_t = 0;
	ERK_t = 0;
	JNK_t = 0;
    case 3
	FRK_t = interpn( kinaseSignals(:,1), kinaseSignals(:,2), t );
	ERK_t = interpn( kinaseSignals(:,3), kinaseSignals(:,4), t );
	JNK_t = interpn( kinaseSignals(:,5), kinaseSignals(:,6), t );
    otherwise
	error('Insufficient Arguments');
end


k(1) = 0.154;
k(2) = 8.4053;
k(3) = 0.2741;
k(4) = 5.1573;
k(5) = 3.8786;

k(6)  = 1.1317;
k(7)  = 0.073;
k(8)  = 0.2705;
k(9)  = 4.5337;
k(10) = 6.5615;

k(11) = 69.5154;
k(12) = 0.0214;
k(13) = 1.3943;
k(14) = 11.9303;
k(15) = 0.0053;

k(16) = 0.00090396;
k(17) = 0.4058;
k(18) = 17.7672;
k(19) = 1.3369;
k(20) = 8.4851;

k(21) = 0.2867;
k(22) = 0.0774;
k(23) = 55.9578;
k(24) = 0.0546;
k(25) = 0.0498;

k(26) = 105.164;
k(27) = 0.9128;
k(28) = 125.366;
k(29) = 3.4048;
k(30) = 0.4247;

k(32) = 0.3316;
k(33) = 0.1557;
k(34) = 16.6871;
k(35) = 2.6297;

k(36) = 0.2352;
k(37) = 0.4267;
k(39) = 0.6621;
k(40) = 0.068;

k(41) = 59.9722;
k(42) = 2.5531;
k(43) = 0.0224;
k(44) = 0.0042;
k(45) = 1.1333;

k(46) = 0.0042;
k(47) = 0.06;
k(48) = 0.0031;

kb(21) = 1.2753;
kb(24) = 428.219;
kb(27) = 1.7441;
kb(29) = 3.1303;
kb(36) = 14.6428;
kb(43) = 4.3273;
kb(45) = 0.6137;

K(1)  = 1.6189;
K(2)  = 6.8804;
K(4)  = 3.0176;
K(5)  = 15.0773;
K(8)  = 65.2625;
K(9)  = 9.534;
K(10) = 27.1651;
K(11) = 112.334;
K(13) = 2.069;
K(14) = 110.623;
K(17) = 18.5934;
K(18) = 1.7933;
K(19) = 21.5017;
K(20) = 3.3207;
K(23) = 8.947;
K(26) = 60.7342;
K(28) = 56.3667;

R(31) = 0.000082285;
R(38) = 0.00081524;


V_N = 14.1*10^(-9);  % nucleus volume in microliters
V_C = 65.3*10^(-9);  % cell volume in microliters

v = Calculate_v( y0, k, kb, K, R, FRK_t, ERK_t, JNK_t);

% dydt or d[activity of species] in time interval dt
dydt(1) = v(2) + v(23) + (V_C/V_N)*v(35) - v(1) - v(3);
dydt(2) = v(1) + v(5) - v(2) - v(4) - v(6);
dydt(3) = v(4) - v(5) - v(7) - v(21);
dydt(4) = v(9) - v(8);
dydt(5) = v(8) - v(9) - v(29);

dydt(6) = v(11) + v(23) + 2*v(26) + v(28) + (V_C/V_N)*v(42) - v(10) - v(12);
dydt(7) = v(10) + v(14) - v(11) - v(13) - v(15);
dydt(8) = v(13) - v(14) - v(16) - v(21) - 2*v(24) - v(27);
dydt(9) = v(18) + v(28) - v(17);
dydt(10) = v(17) + v(20) - v(18) - v(19);

dydt(11) = v(19) - v(20) - v(27);
dydt(12) = v(21) - v(22) - v(23) - v(45);
dydt(13) = v(24) - v(25) - v(26) - v(43);
dydt(14) = v(27) - v(28) - v(36);
dydt(15) = -v(29);

dydt(16) = v(29);
dydt(17) = v(31) + v(30) - (V_C/V_N)*v(32);
dydt(18) = v(32) - v(33);
dydt(19) = v(34) - v(35);
dydt(20) = -v(36);

dydt(21) = v(36);
dydt(22) = v(37) + v(38) - (V_C/V_N)*v(39);
dydt(23) = v(39) - v(40);
dydt(24) = v(41) - v(42);
dydt(25) = -v(43) - v(45);

dydt(26) = v(43);
dydt(27) = v(45);
dydt(28) = v(44) + v(46) - (V_C/V_N)*v(47);
dydt(29) = v(47) - v(48);
dydt(30) = dydt(12) + dydt(13);

clear v

% matlab output format. Else will give error.
dydt = dydt(:);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YY = interpn( X, Y, XX )
% custom interpolation function

  I = find( X >= XX);
  i = I(1)-1;
  if i==0
  i=1;
  end
  YY = Y(i) + ( ( XX - X(i) )*Y(i+1) - ( XX - X(i) )*Y(i)  ) / ( X(i+1) - X(i) );
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = Calculate_v( y0, k, kb, K, R, FRK_t, ERK_t, JNK_t)
% Calculates reaction rates vector at time t
% t - time
% y0 - species conc at time t
% k, kb, K, R - parse rate constants
% FRK_t, ERK_t, JNK_t - parse the signal activity levels


V_N = 14.1*10^(-9); % microliters
V_C = 65.3*10^(-9); % microliters

% Calculate Reaction Rates
v(1) = michaelismenten( k(1), K(1), [y0(1),FRK_t] );
v(2) = michaelismenten( k(2), K(2), y0(2) );
v(3) = reactionrate( k(3), y0(1) );
v(4) = michaelismenten( k(4), K(4), [y0(2),FRK_t] );
v(5) = michaelismenten( k(5), K(5), y0(3) );

v(6) = reactionrate( k(6), y0(2) );
v(7) = reactionrate( k(7), y0(3) );
v(8) = michaelismenten( k(8), K(8), [y0(4),ERK_t] );
v(9) = michaelismenten( k(9), K(9), y0(5) );
v(10) = michaelismenten( k(10), K(10), [y0(6),JNK_t] );

v(11) = michaelismenten( k(11), K(11), y0(7) );
v(12) = reactionrate( k(12), y0(6) );
v(13) = michaelismenten( k(13), K(13), [y0(7),JNK_t] );
v(14) = michaelismenten( k(14), K(14), y0(8) );
v(15) = reactionrate( k(15), y0(7) );

v(16) = reactionrate( k(16), y0(8) );
v(17) = michaelismenten( k(17), K(17), [y0(9),JNK_t] );
v(18) = michaelismenten( k(18), K(18), y0(10) );
v(19) = michaelismenten( k(19), K(19), [y0(10),JNK_t] );
v(20) = michaelismenten( k(20), K(20), y0(11) );

v(21) = reactionrate( [k(21),kb(21)], [y0(3),y0(8),y0(12)] );
v(22) = reactionrate( k(22), y0(12) );
v(23) = michaelismenten( k(23), K(23), y0(12) );
v(24) = reactionrate( [k(24),kb(24)], [y0(8),y0(8),y0(13)] );
v(25) = reactionrate( k(25), y0(13) );

v(26) = michaelismenten( k(26), K(26), y0(13) );
v(27) = reactionrate( [k(27),kb(27)], [y0(11),y0(8),y0(14)] );
v(28) = michaelismenten( k(28), K(28), y0(14) );
v(29) = reactionrate( [k(29),kb(29)], [y0(15),y0(5),y0(16)] );
v(30) = reactionrate( k(30), y0(16) );

v(31) = R(31);
v(32) = reactionrate( (k(32) * (V_N / V_C) ), y0(17) );
v(33) = reactionrate( k(33), y0(18) );
v(34) = reactionrate( k(34), y0(18) );
v(35) = reactionrate( k(35), y0(19) );

v(36) = reactionrate( [k(36),kb(36)], [y0(20),y0(14),y0(21)] );
v(37) = reactionrate( k(37), y0(21) );
v(38) = R(38);
v(39) = reactionrate( (k(39) * (V_N / V_C) ), y0(22) );
v(40) = reactionrate( k(40), y0(23) );

v(41) = reactionrate( k(41), y0(23) );
v(42) = reactionrate( k(42), y0(24) );
v(43) = reactionrate( [k(43),kb(43)], [y0(25),y0(13),y0(26)] );
v(44) = reactionrate( k(44), y0(26) );
v(45) = reactionrate( [k(45),kb(45)], [y0(25),y0(12),y0(27)] );

v(46) = reactionrate( k(46),y0(27) );
v(47) = reactionrate( ( k(47) * (V_N / V_C) ), y0(28) );
v(48) = reactionrate( k(48), y0(29) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = michaelismenten( k,K,conc )
%Calculates Michaelis Menten Reaction rate

n = length( conc );

if n==1
    v = ( k*conc(1) ) / ( K + conc(1) );
elseif n==2
    v=( k*conc(1)*conc(2) ) / ( K + conc(1) );
else 
    error( 'Wrong number of input arguments' )
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = reactionrate( k,conc )
% Calculates first order rection rate

n = length( conc );
if n==1
    v = k*conc(1) ;
elseif n==3
    v = k(1)*conc(1)*conc(2) - k(2)*conc(3) ;
else 
    error( 'Wrong number of input arguments' )
end

end
