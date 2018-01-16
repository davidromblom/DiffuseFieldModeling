function [ virMicDir numLdspkr ] = getLoudspeakerConfiguration(ldspkrSystem)

switch ldspkrSystem
    case 2,
        % 16.0 setup
        virMicDir   = [
                (0.0/180)*pi,  0.0 ; ...
                (22.5/180)*pi,  0.0 ; ...
                (45/180)*pi, 0.0 ; ... 
                (67.5/180)*pi,  0.0 ; ...
                (90/180)*pi, 0.0 ; ...
                (112.5/180)*pi, 0.0 ; ...
                (135/180)*pi, 0.0 ; ... 
                (157.5/180)*pi, 0.0; ...
                (180/180)*pi, 0.0 ; ...
                -(157.5/180)*pi, 0.0 ; ...
                -(135/180)*pi, 0.0 ; 
                -(112.5/180)*pi, 0.0 ; ...
                -(90/180)*pi, 0.0 ; ...
                -(67.5/180)*pi, 0.0 ; ...
                -(45/180)*pi, 0.0; ...
                -(22.5/180)*pi, 0.0 ...
        ];
        numLdspkr = 16;
        
    case 3,
        % 9.1 cube
        virMicDir   = [ 
                (30/180)*pi,  0.0 ; ...     % Left
                -(30/180)*pi, 0.0 ; ...     % Right
                (110/180)*pi,  0.0 ; ...    % Left Surround
                -(110/180)*pi, 0.0 ; ...    % Right Surround
                 (0/180)*pi, 0.0 ; ...      % Center 
                (30/180)*pi, (40/180)*pi ; ...     % Left Height
                -(30/180)*pi,(40/180)*pi ; ...     % Right Height
                (110/180)*pi, (40/180)*pi ; ...    % LS Height
                -(110/180)*pi, (40/180)*pi ; ...   % RS Height
        ];
        numLdspkr = 4;
        
     case 4,
        % Faux Binaural
        virMicDir = [ -(100/180)*pi,  0.0 ; (100/180)*pi, 0.0 ];
        numLdspkr = 2;
        
    otherwise,
        disp('That loudspeaker spacing was not found');

end;

return;


