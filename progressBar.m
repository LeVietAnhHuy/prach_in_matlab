imax = 100;
prog = 0;
fprintf('Computation Progress: %3d%%\n',prog);
for k = 1:1:imax
	prog = ( 100*(k/imax) );
	fprintf('\b\b\b\b%3.0f%%',prog); pause(0.1); % Deleting 4 characters (The three digits and the % symbol)
end
fprintf('\n'); % To go to a new line after reaching 100% progress