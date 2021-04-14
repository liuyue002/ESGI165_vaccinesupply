function [position,isterminal,direction] = event_negative(t,X,M)


position = X(1:M);
isterminal = ones(M,1);  
direction = 0; 