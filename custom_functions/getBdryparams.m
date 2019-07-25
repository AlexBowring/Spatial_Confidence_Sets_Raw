function [bdry_params] = getBdryparams(field, c, connectivity, a)

% Finds edges and weights needed to linearly interpolate a discretized field 
% to find the locations where the true continuous field = c
%
% Input:
% field:     random field over a domain in R^2, it is an (2+1)-dimensional array,
%            where the last dimension enumerates realisations
% c:         threshold value for excursion
% connectivity: A binary value, if 0 then 6cc is used in 3D case, if 1 then 26cc used. 
%
%Output:
% bdry_params is a struct containing the edge locations either side of the
% true continuous boundary where the field = c, as well as the weights that
% the adjA_cent voxels need to be multiplied by to give the precise location where the 
% field = c assuming that the gradient of the signal is linear between adjA_cent
% voxels. 
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Alex Bowring (alex.bowring@bdi.ox.A_c.uk)
% Last changes: 10/25/2018
%__________________________________________________________________________


A_c = (field >= c);
dim = size(A_c);
D   = length(dim);

switch D
    case 2
        %%%%%%%%%%%%%%%% Case 2D random field with 4-connectivity %%%%%%%%%%%%%%%%%
        %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the second component might be called vertical)
        horz = A_c(:,2:end) | A_c(:,1:end-1);
        %%% Compute the left shifted horizontal edges
        lshift            = A_c; % initialize
        lshift(:,1:end-1) = horz;
        lshift            = lshift & ~A_c;
        % Signal locations to be used with lshift edges
        lshift_sig        = lshift(:,[dim(2) 1:dim(2)-1]);
        %%% Compute the right shifted horizontal edges
        rshift          = A_c; % initialize
        rshift(:,2:end) = horz;
        rshift          = rshift & ~A_c;
        % Signal locations to be used with rshift edges
        rshift_sig      = rshift(:,[2:dim(2) 1]);
        %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the first component might be called horizontal)
        vert = A_c(1:end-1,:) | A_c(2:end,:);
        %%% Compute the right shifted horizontal edges
        ushift = A_c;
        ushift(1:end-1,:) = vert;
        ushift = ushift & ~A_c;
        % Signal locations to be used with ushift edges
        ushift_sig      = ushift([dim(1) 1:dim(1)-1],:);
        %%% Compute the down shifted vertical edges
        dshift = A_c;
        dshift(2:end,:)   = vert;
        dshift = dshift & ~A_c;
        % Signal locations to be used with dshift edges
        dshift_sig      = dshift([2:dim(1) 1],:);
        
        % Computing the weights for the weighted linear boundary of the
        % field
        lshift_w1 = abs(field(lshift_sig) - c)./abs(field(lshift) - field(lshift_sig));
        lshift_w2 = abs(field(lshift) - c)./abs(field(lshift) - field(lshift_sig));

        rshift_w1 = abs(field(rshift_sig) - c)./abs(field(rshift) - field(rshift_sig));
        rshift_w2 = abs(field(rshift) - c)./abs(field(rshift) - field(rshift_sig));

        ushift_w1 = abs(field(ushift_sig) - c)./abs(field(ushift) - field(ushift_sig));
        ushift_w2 = abs(field(ushift) - c)./abs(field(ushift) - field(ushift_sig));

        dshift_w1 = abs(field(dshift_sig) - c)./abs(field(dshift) - field(dshift_sig));
        dshift_w2 = abs(field(dshift) - c)./abs(field(dshift) - field(dshift_sig));
        
        switch a 
            case 0
                % Compute the length of the boundary
                len    = length(lshift_w1) + length(rshift_w1) + length(ushift_w1) + length(dshift_w1);
                
                % Storing parameters in a structure
                bdry_params = struct('length', len, ...
                                     'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                                     'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                                     'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                                     'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2));

            case 1
                lshift_w2(lshift_w1 < 0.5) = 0; 
                lshift_w1(lshift_w1 < 0.5) = 1; 
                
                rshift_w2(rshift_w1 < 0.5) = 0; 
                rshift_w1(rshift_w1 < 0.5) = 1; 
                
                ushift_w2(ushift_w1 < 0.5) = 0; 
                ushift_w1(ushift_w1 < 0.5) = 1; 
                
                dshift_w2(dshift_w1 < 0.5) = 0; 
                dshift_w1(dshift_w1 < 0.5) = 1; 
                                
                % Compute the length of the boundary
                len    = length(lshift_w1) + length(rshift_w1) + length(ushift_w1) + length(dshift_w1);
                
                % Storing parameters in a structure
                bdry_params = struct('length', len, ...
                                     'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                                     'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                                     'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                                     'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2));
        end          
    case 3

        switch connectivity
            case 0 
            %%%%%%%%%%%%%%%% Case 3D random field with 6-connectivity %%%%%%%%%%%%%%%%%
            %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
            %%%%%% usually the second component might be called vertical)
            horz = A_c(:,2:end,:) | A_c(:,1:end-1,:);
            %%% Compute the left shifted horizontal edges
            lshift              = A_c; % initialize
            lshift(:,1:end-1,:) = horz;
            lshift              = lshift & ~A_c;
            % Signal locations to be used with lshift edges
            lshift_sig          = lshift(:,[dim(2) 1:dim(2)-1],:);
            %%% Compute the right shifted horizontal edges
            rshift            = A_c; % initialize
            rshift(:,2:end,:) = horz;
            rshift            = rshift & ~A_c;
            % Signal locations to be used with rshift edges
            rshift_sig        = rshift(:,[2:dim(2) 1],:);
            %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
            %%%%%% usually the first component might be called horizontal)
            vert = A_c(1:end-1,:,:) | A_c(2:end,:,:);
            %%% Compute the up shifted vertical edges
            ushift = A_c;
            ushift(1:end-1,:,:) = vert;
            ushift = ushift & ~A_c;
            % Signal locations to be used with ushift edges
            ushift_sig = ushift([dim(1) 1:dim(1)-1],:,:);
            %%% Compute the down shifted vertical edges
            dshift = A_c;
            dshift(2:end,:,:)   = vert;
            dshift = dshift & ~A_c;
            % Signal locations to be used with dshift edges
            dshift_sig = dshift([2:dim(1) 1],:,:);
            %%%%%% depth edges
            depth = A_c(:,:,1:end-1) | A_c(:,:,2:end);
            %%% Compute the bA_ck shifted depth edges
            bshift = A_c;
            bshift(:,:,1:end-1) = depth;
            bshift = bshift & ~A_c;
            % Signal locations to be used with bshift edges
            bshift_sig = bshift(:,:,[dim(3) 1:dim(3)-1]);
            %%% Compute the front shifted depth edges
            fshift = A_c;
            fshift(:,:,2:end)   = depth;
            fshift = fshift & ~A_c;
            % Signal locations to be used with fshift edges
            fshift_sig = fshift(:,:,[2:dim(3) 1]);

            % Computing the weights for the weighted linear boundary of the
            % field
            lshift_w1 = abs(field(lshift_sig) - c)./abs(field(lshift) - field(lshift_sig));
            lshift_w2 = abs(field(lshift) - c)./abs(field(lshift) - field(lshift_sig));

            rshift_w1 = abs(field(rshift_sig) - c)./abs(field(rshift) - field(rshift_sig));
            rshift_w2 = abs(field(rshift) - c)./abs(field(rshift) - field(rshift_sig));

            ushift_w1 = abs(field(ushift_sig) - c)./abs(field(ushift) - field(ushift_sig));
            ushift_w2 = abs(field(ushift) - c)./abs(field(ushift) - field(ushift_sig));

            dshift_w1 = abs(field(dshift_sig) - c)./abs(field(dshift) - field(dshift_sig));
            dshift_w2 = abs(field(dshift) - c)./abs(field(dshift) - field(dshift_sig));

            bshift_w1 = abs(field(bshift_sig) - c)./abs(field(bshift) - field(bshift_sig));
            bshift_w2 = abs(field(bshift) - c)./abs(field(bshift) - field(bshift_sig));

            fshift_w1 = abs(field(fshift_sig) - c)./abs(field(fshift) - field(fshift_sig));
            fshift_w2 = abs(field(fshift) - c)./abs(field(fshift) - field(fshift_sig));
            
            
            switch a 
                case 0
                    % Compute the length of the boundary
                    len    = length(lshift_w1) + length(rshift_w1) + length(ushift_w1) + length(dshift_w1) ...
                               + length(bshift_w1) + length(fshift_w1);


                    % Create structure for storing parameters
                    bdry_params = struct('length', len, ...
                                         'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                                         'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                                         'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                                         'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2), ...
                                         'bshift', struct('edges', bshift, 'sig_edges', bshift_sig, 'w1', bshift_w1,'w2', bshift_w2), ...
                                         'fshift', struct('edges', fshift, 'sig_edges', fshift_sig, 'w1', fshift_w1,'w2', fshift_w2)); 

                case 1
                    lshift_w2(lshift_w1 < 0.5) = 0; 
                    lshift_w1(lshift_w1 < 0.5) = 1; 

                    rshift_w2(rshift_w1 < 0.5) = 0; 
                    rshift_w1(rshift_w1 < 0.5) = 1; 

                    ushift_w2(ushift_w1 < 0.5) = 0; 
                    ushift_w1(ushift_w1 < 0.5) = 1; 

                    dshift_w2(dshift_w1 < 0.5) = 0; 
                    dshift_w1(dshift_w1 < 0.5) = 1; 
                    
                    bshift_w2(bshift_w1 < 0.5) = 0; 
                    bshift_w1(bshift_w1 < 0.5) = 1; 
                    
                    fshift_w2(fshift_w1 < 0.5) = 0; 
                    fshift_w1(fshift_w1 < 0.5) = 1;
                    
                    % Compute the length of the boundary
                    len    = length(lshift_w1) + length(rshift_w1) + length(ushift_w1) + length(dshift_w1) ...
                               + length(bshift_w1) + length(fshift_w1);


                    % Create structure for storing parameters
                    bdry_params = struct('length', len, ...
                                         'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                                         'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                                         'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                                         'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2), ...
                                         'bshift', struct('edges', bshift, 'sig_edges', bshift_sig, 'w1', bshift_w1,'w2', bshift_w2), ...
                                         'fshift', struct('edges', fshift, 'sig_edges', fshift_sig, 'w1', fshift_w1,'w2', fshift_w2)); 
            end
            
            case 1
            %%%%%%%%%%%%%%%% Case 3D random field with 26-connectivity %%%%%%%%%%%%%%%%%
            %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
            %%%%%% usually the second component might be called vertical)
            horz = A_c(:,2:end,:) | A_c(:,1:end-1,:);
            %%% Compute the left shifted horizontal edges
            lshift                                = A_c; % initialize
            lshift(:,1:end-1,:)                   = horz;
            lshift                                = lshift & ~A_c;
            % Signal locations to be used with lshift edges
            lshift_sig                            = lshift(:,[dim(2) 1:dim(2)-1],:);
            %% Compute left and vertically shifted edges
            lshift_vert                           = lshift(1:end-1,:,:) | lshift(2:end,:,:);
            % Left and up shifted
            lshift_ushift                         = lshift;
            lshift_ushift(1:end-1,:,:)            = lshift_vert;
            lshift_ushift                         = lshift_ushift & ~A_c & ~lshift;
            % Signal locations to be used with lshift_ushift edges
            lshift_ushift_sig                     = lshift_ushift([dim(1) 1:dim(1)-1],[dim(2) 1:dim(2)-1],:); 
            % Left and down shifted
            lshift_dshift                         = lshift;
            lshift_dshift(2:end,:,:)              = lshift_vert;
            lshift_dshift                         = lshift_dshift & ~A_c & ~lshift;
            % Signal locations to be used with lshift_dshift edges
            lshift_dshift_sig                     = lshift_dshift([2:dim(1) 1],[dim(2) 1:dim(2)-1],:);
            %% Compute left and depth shifted edges
            lshift_depth                          = lshift(:,:,1:end-1) | lshift(:,:,2:end);
            % Left and bA_ck shifted
            lshift_bshift                         = lshift;
            lshift_bshift(:,:,1:end-1)            = lshift_depth;
            lshift_bshift                         = lshift_bshift & ~A_c & ~lshift;
            % Signal locations to be used with lshift_bshift edges
            lshift_bshift_sig                     = lshift_bshift(:,[dim(2) 1:dim(2)-1],[dim(3) 1:dim(3)-1]);
            % Left and front shifted 
            lshift_fshift                         = lshift;
            lshift_fshift(:,:,2:end)              = lshift_depth;
            lshift_fshift                         = lshift_fshift & ~A_c & ~lshift;
            % Signal locations to be used with lshift_fshift edges
            lshift_fshift_sig                     = lshift_fshift(:,[dim(2) 1:dim(2)-1],[2:dim(3) 1]);
            %% Compute left, bA_ck, and vertically shifted edges
            lshift_bshift_vert                    = lshift_bshift(1:end-1,:,:) | lshift_bshift(2:end,:,:);
            % Left, bA_ck, and up shifted
            lshift_bshift_ushift                  = lshift_bshift;
            lshift_bshift_ushift(1:end-1,:,:)     = lshift_bshift_vert;
            lshift_bshift_ushift                  = lshift_bshift_ushift & ~A_c & ~lshift_bshift;
            % Signal locations to be used with lshift_bshift_ushift edges
            lshift_bshift_ushift_sig              = lshift_bshift_ushift([dim(1) 1:dim(1)-1],[dim(2) 1:dim(2)-1],[dim(3) 1:dim(3)-1]);
            % Left, bA_ck, and down shifted
            lshift_bshift_dshift                  = lshift_bshift;
            lshift_bshift_dshift(2:end,:,:)       = lshift_bshift_vert;
            lshift_bshift_dshift                  = lshift_bshift_dshift & ~A_c & ~lshift_bshift;
            % Signal locations to be used with lshift_bshift_dshift edges
            lshift_bshift_dshift_sig              = lshift_bshift_dshift([2:dim(1) 1],[dim(2) 1:dim(2)-1],[dim(3) 1:dim(3)-1]); 
            %% Compute left, front, and vertically shifted edges
            lshift_fshift_vert                    = lshift_fshift(1:end-1,:,:) | lshift_fshift(2:end,:,:);
            % Left, front, and up shifted
            lshift_fshift_ushift                  = lshift_fshift;
            lshift_fshift_ushift(1:end-1,:,:)     = lshift_fshift_vert;
            lshift_fshift_ushift                  = lshift_fshift_ushift & ~A_c & ~lshift_fshift;
            % Signal locations to be used with lshift_fshift_ushift edges
            lshift_fshift_ushift_sig              = lshift_fshift_ushift([dim(1) 1:dim(1)-1],[dim(2) 1:dim(2)-1],[2:dim(3) 1]);
            % Left, front, and down shifted
            lshift_fshift_dshift                  = lshift_fshift;
            lshift_fshift_dshift(2:end,:,:)       = lshift_fshift_vert;
            lshift_fshift_dshift                  = lshift_fshift_dshift & ~A_c & ~lshift_fshift;
            % Signal locations to be used with lshift_fshift_dshift edges
            lshift_fshift_dshift_sig              = lshift_fshift_dshift([2:dim(1) 1],[dim(2) 1:dim(2)-1],[2:dim(3) 1]);  
            
            %%% Compute the right shifted horizontal edges
            rshift            = A_c; % initialize
            rshift(:,2:end,:) = horz;
            rshift            = rshift & ~A_c;
            % Signal locations to be used with rshift edges
            rshift_sig        = rshift(:,[2:dim(2) 1],:);
            %% Compute right and vertically shifted edges
            rshift_vert                           = rshift(1:end-1,:,:) | rshift(2:end,:,:);
            % Right and up shifted
            rshift_ushift                         = rshift;
            rshift_ushift(1:end-1,:,:)            = rshift_vert;
            rshift_ushift                         = rshift_ushift & ~A_c & ~rshift;
            % Signal locations to be used with rshift_ushift edges
            rshift_ushift_sig                     = rshift_ushift([dim(1) 1:dim(1)-1],[2:dim(2) 1],:); 
            % Right and down shifted
            rshift_dshift                         = rshift;
            rshift_dshift(2:end,:,:)              = rshift_vert;
            rshift_dshift                         = rshift_dshift & ~A_c & ~rshift;
            % Signal locations to be used with rshift_dshift edges
            rshift_dshift_sig                     = rshift_dshift([2:dim(1) 1],[2:dim(2) 1],:);
            %% Compute right and depth shifted edges
            rshift_depth                          = rshift(:,:,1:end-1) | rshift(:,:,2:end);
            % Right and bA_ck shifted
            rshift_bshift                         = rshift;
            rshift_bshift(:,:,1:end-1)            = rshift_depth;
            rshift_bshift                         = rshift_bshift & ~A_c & ~rshift;
            % Signal locations to be used with rshift_bshift edges
            rshift_bshift_sig                     = rshift_bshift(:,[2:dim(2) 1],[dim(3) 1:dim(3)-1]);
            % Right and front shifted 
            rshift_fshift                         = rshift;
            rshift_fshift(:,:,2:end)              = rshift_depth;
            rshift_fshift                         = rshift_fshift & ~A_c & ~rshift;
            % Signal locations to be used with rshift_fshift edges
            rshift_fshift_sig                     = rshift_fshift(:,[2:dim(2) 1],[2:dim(3) 1]); 
            %% Compute right, bA_ck, and vertically shifted edges
            rshift_bshift_vert                    = rshift_bshift(1:end-1,:,:) | rshift_bshift(2:end,:,:);
            % right, bA_ck, and up shifted
            rshift_bshift_ushift                  = rshift_bshift;
            rshift_bshift_ushift(1:end-1,:,:)     = rshift_bshift_vert;
            rshift_bshift_ushift                  = rshift_bshift_ushift & ~A_c & ~rshift_bshift;
            % Signal locations to be used with rshift_bshift_ushift edges
            rshift_bshift_ushift_sig              = rshift_bshift_ushift([dim(1) 1:dim(1)-1],[2:dim(2) 1],[dim(3) 1:dim(3)-1]);
            % right, bA_ck, and down shifted
            rshift_bshift_dshift                  = rshift_bshift;
            rshift_bshift_dshift(2:end,:,:)       = rshift_bshift_vert;
            rshift_bshift_dshift                  = rshift_bshift_dshift & ~A_c & ~rshift_bshift;
            % Signal locations to be used with rshift_bshift_dshift edges
            rshift_bshift_dshift_sig              = rshift_bshift_dshift([2:dim(1) 1],[2:dim(2) 1],[dim(3) 1:dim(3)-1]); 
            %% Compute right, front, and vertically shifted edges
            rshift_fshift_vert                    = rshift_fshift(1:end-1,:,:) | rshift_fshift(2:end,:,:);
            % Right, front, and up shifted
            rshift_fshift_ushift                  = rshift_fshift;
            rshift_fshift_ushift(1:end-1,:,:)     = rshift_fshift_vert;
            rshift_fshift_ushift                  = rshift_fshift_ushift & ~A_c & ~rshift_fshift;
            % Signal locations to be used with rshift_fshift_ushift edges
            rshift_fshift_ushift_sig              = rshift_fshift_ushift([dim(1) 1:dim(1)-1],[2:dim(2) 1],[2:dim(3) 1]);
            % Right, front, and down shifted
            rshift_fshift_dshift                  = rshift_fshift;
            rshift_fshift_dshift(2:end,:,:)       = rshift_fshift_vert;
            rshift_fshift_dshift                  = rshift_fshift_dshift & ~A_c & ~rshift_fshift;
            % Signal locations to be used with rshift_fshift_dshift edges
            rshift_fshift_dshift_sig              = rshift_fshift_dshift([2:dim(1) 1],[2:dim(2) 1],[2:dim(3) 1]); 

            %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
            %%%%%% usually the first component might be called horizontal)
            vert = A_c(1:end-1,:,:) | A_c(2:end,:,:);
            %%% Compute the up shifted horizontal edges
            ushift = A_c;
            ushift(1:end-1,:,:) = vert;
            ushift = ushift & ~A_c;
            % Signal locations to be used with ushift edges
            ushift_sig = ushift([dim(1) 1:dim(1)-1],:,:);
            %% Compute up and depth shifted edges
            ushift_depth                          = ushift(:,:,1:end-1) | ushift(:,:,2:end);
            % Up and bA_ck shifted
            ushift_bshift                         = ushift;
            ushift_bshift(:,:,1:end-1)            = ushift_depth;
            ushift_bshift                         = ushift_bshift & ~A_c & ~ushift;
            % Signal locations to be used with ushift_bshift edges
            ushift_bshift_sig                     = ushift_bshift([dim(1) 1:dim(1)-1],:,[dim(3) 1:dim(3)-1]);
            % Up and front shifted 
            ushift_fshift                         = ushift;
            ushift_fshift(:,:,2:end)              = ushift_depth;
            ushift_fshift                         = ushift_fshift & ~A_c & ~ushift;
            % Signal locations to be used with ushift_fshift edges
            ushift_fshift_sig                     = ushift_fshift([dim(1) 1:dim(1)-1],:,[2:dim(3) 1]); 

            %%% Compute the down shifted vertical edges
            dshift = A_c;
            dshift(2:end,:,:)   = vert;
            dshift = dshift & ~A_c;
            % Signal locations to be used with dshift edges 
            dshift_sig = dshift([2:dim(1) 1],:,:);
            %% Compute down and depth shifted edges
            dshift_depth                          = dshift(:,:,1:end-1) | dshift(:,:,2:end);
            % Down and bA_ck shifted
            dshift_bshift                         = dshift;
            dshift_bshift(:,:,1:end-1)            = dshift_depth;
            dshift_bshift                         = dshift_bshift & ~A_c & ~dshift;
            % Signal locations to be used with dshift_bshift edges
            dshift_bshift_sig                     = dshift_bshift([2:dim(1) 1],:,[dim(3) 1:dim(3)-1]);
            % Down and front shifted 
            dshift_fshift                         = dshift;
            dshift_fshift(:,:,2:end)              = dshift_depth;
            dshift_fshift                         = dshift_fshift & ~A_c & ~dshift;
            % Signal locations to be used with dshift_fshift edges
            dshift_fshift_sig                     = dshift_fshift([2:dim(1) 1],:,[2:dim(3) 1]); 

            %%%%%% depth edges
            depth = A_c(:,:,1:end-1) | A_c(:,:,2:end);
            %%% Compute the bA_ck shifted depth edges
            bshift = A_c;
            bshift(:,:,1:end-1) = depth;
            bshift = bshift & ~A_c;
            % Signal locations to be used with bshift edges
            bshift_sig = bshift(:,:,[dim(3) 1:dim(3)-1]);
            %%% Compute the front shifted depth edges
            fshift = A_c;
            fshift(:,:,2:end)   = depth;
            fshift = fshift & ~A_c;
            % Signal locations to be used with fshift edges
            fshift_sig = fshift(:,:,[2:dim(3) 1]);

            % Computing the weights for the weighted linear boundary of the
            % field
            lshift_w1 = abs(field(lshift_sig) - c)./abs(field(lshift) - field(lshift_sig));
            lshift_w2 = abs(field(lshift) - c)./abs(field(lshift) - field(lshift_sig));

            lshift_ushift_w1 = abs(field(lshift_ushift_sig) - c)./abs(field(lshift_ushift) - field(lshift_ushift_sig));
            lshift_ushift_w2 = abs(field(lshift_ushift) - c)./abs(field(lshift_ushift) - field(lshift_ushift_sig));

            lshift_dshift_w1 = abs(field(lshift_dshift_sig) - c)./abs(field(lshift_dshift) - field(lshift_dshift_sig));
            lshift_dshift_w2 = abs(field(lshift_dshift) - c)./abs(field(lshift_dshift) - field(lshift_dshift_sig));

            lshift_bshift_w1 = abs(field(lshift_bshift_sig) - c)./abs(field(lshift_bshift) - field(lshift_bshift_sig));
            lshift_bshift_w2 = abs(field(lshift_bshift) - c)./abs(field(lshift_bshift) - field(lshift_bshift_sig));

            lshift_fshift_w1 = abs(field(lshift_fshift_sig) - c)./abs(field(lshift_fshift) - field(lshift_fshift_sig));
            lshift_fshift_w2 = abs(field(lshift_fshift) - c)./abs(field(lshift_fshift) - field(lshift_fshift_sig));

            lshift_bshift_ushift_w1 = abs(field(lshift_bshift_ushift_sig) - c)./abs(field(lshift_bshift_ushift) - field(lshift_bshift_ushift_sig));
            lshift_bshift_ushift_w2 = abs(field(lshift_bshift_ushift) - c)./abs(field(lshift_bshift_ushift) - field(lshift_bshift_ushift_sig));

            lshift_bshift_dshift_w1 = abs(field(lshift_bshift_dshift_sig) - c)./abs(field(lshift_bshift_dshift) - field(lshift_bshift_dshift_sig));
            lshift_bshift_dshift_w2 = abs(field(lshift_bshift_dshift) - c)./abs(field(lshift_bshift_dshift) - field(lshift_bshift_dshift_sig));

            lshift_fshift_ushift_w1 = abs(field(lshift_fshift_ushift_sig) - c)./abs(field(lshift_fshift_ushift) - field(lshift_fshift_ushift_sig));
            lshift_fshift_ushift_w2 = abs(field(lshift_fshift_ushift) - c)./abs(field(lshift_fshift_ushift) - field(lshift_fshift_ushift_sig));

            lshift_fshift_dshift_w1 = abs(field(lshift_fshift_dshift_sig) - c)./abs(field(lshift_fshift_dshift) - field(lshift_fshift_dshift_sig));
            lshift_fshift_dshift_w2 = abs(field(lshift_fshift_dshift) - c)./abs(field(lshift_fshift_dshift) - field(lshift_fshift_dshift_sig));

            rshift_w1 = abs(field(rshift_sig) - c)./abs(field(rshift) - field(rshift_sig));
            rshift_w2 = abs(field(rshift) - c)./abs(field(rshift) - field(rshift_sig));

            rshift_ushift_w1 = abs(field(rshift_ushift_sig) - c)./abs(field(rshift_ushift) - field(rshift_ushift_sig));
            rshift_ushift_w2 = abs(field(rshift_ushift) - c)./abs(field(rshift_ushift) - field(rshift_ushift_sig));

            rshift_dshift_w1 = abs(field(rshift_dshift_sig) - c)./abs(field(rshift_dshift) - field(rshift_dshift_sig));
            rshift_dshift_w2 = abs(field(rshift_dshift) - c)./abs(field(rshift_dshift) - field(rshift_dshift_sig));

            rshift_bshift_w1 = abs(field(rshift_bshift_sig) - c)./abs(field(rshift_bshift) - field(rshift_bshift_sig));
            rshift_bshift_w2 = abs(field(rshift_bshift) - c)./abs(field(rshift_bshift) - field(rshift_bshift_sig));

            rshift_fshift_w1 = abs(field(rshift_fshift_sig) - c)./abs(field(rshift_fshift) - field(rshift_fshift_sig));
            rshift_fshift_w2 = abs(field(rshift_fshift) - c)./abs(field(rshift_fshift) - field(rshift_fshift_sig));

            rshift_bshift_ushift_w1 = abs(field(rshift_bshift_ushift_sig) - c)./abs(field(rshift_bshift_ushift) - field(rshift_bshift_ushift_sig));
            rshift_bshift_ushift_w2 = abs(field(rshift_bshift_ushift) - c)./abs(field(rshift_bshift_ushift) - field(rshift_bshift_ushift_sig));

            rshift_bshift_dshift_w1 = abs(field(rshift_bshift_dshift_sig) - c)./abs(field(rshift_bshift_dshift) - field(rshift_bshift_dshift_sig));
            rshift_bshift_dshift_w2 = abs(field(rshift_bshift_dshift) - c)./abs(field(rshift_bshift_dshift) - field(rshift_bshift_dshift_sig));

            rshift_fshift_ushift_w1 = abs(field(rshift_fshift_ushift_sig) - c)./abs(field(rshift_fshift_ushift) - field(rshift_fshift_ushift_sig));
            rshift_fshift_ushift_w2 = abs(field(rshift_fshift_ushift) - c)./abs(field(rshift_fshift_ushift) - field(rshift_fshift_ushift_sig));

            rshift_fshift_dshift_w1 = abs(field(rshift_fshift_dshift_sig) - c)./abs(field(rshift_fshift_dshift) - field(rshift_fshift_dshift_sig));
            rshift_fshift_dshift_w2 = abs(field(rshift_fshift_dshift) - c)./abs(field(rshift_fshift_dshift) - field(rshift_fshift_dshift_sig));

            ushift_w1 = abs(field(ushift_sig) - c)./abs(field(ushift) - field(ushift_sig));
            ushift_w2 = abs(field(ushift) - c)./abs(field(ushift) - field(ushift_sig));

            ushift_bshift_w1 = abs(field(ushift_bshift_sig) - c)./abs(field(ushift_bshift) - field(ushift_bshift_sig));
            ushift_bshift_w2 = abs(field(ushift_bshift) - c)./abs(field(ushift_bshift) - field(ushift_bshift_sig));

            ushift_fshift_w1 = abs(field(ushift_fshift_sig) - c)./abs(field(ushift_fshift) - field(ushift_fshift_sig));
            ushift_fshift_w2 = abs(field(ushift_fshift) - c)./abs(field(ushift_fshift) - field(ushift_fshift_sig));

            dshift_w1 = abs(field(dshift_sig) - c)./abs(field(dshift) - field(dshift_sig));
            dshift_w2 = abs(field(dshift) - c)./abs(field(dshift) - field(dshift_sig));

            dshift_bshift_w1 = abs(field(dshift_bshift_sig) - c)./abs(field(dshift_bshift) - field(dshift_bshift_sig));
            dshift_bshift_w2 = abs(field(dshift_bshift) - c)./abs(field(dshift_bshift) - field(dshift_bshift_sig));

            dshift_fshift_w1 = abs(field(dshift_fshift_sig) - c)./abs(field(dshift_fshift) - field(dshift_fshift_sig));
            dshift_fshift_w2 = abs(field(dshift_fshift) - c)./abs(field(dshift_fshift) - field(dshift_fshift_sig));

            bshift_w1 = abs(field(bshift_sig) - c)./abs(field(bshift) - field(bshift_sig));
            bshift_w2 = abs(field(bshift) - c)./abs(field(bshift) - field(bshift_sig));

            fshift_w1 = abs(field(fshift_sig) - c)./abs(field(fshift) - field(fshift_sig));
            fshift_w2 = abs(field(fshift) - c)./abs(field(fshift) - field(fshift_sig));
            
            % Compute the length of the boundary
            len    = length(lshift_w1) + length(lshift_ushift_w1) + length(lshift_dshift_w1) ...
                     + length(lshift_bshift_w1) + length(lshift_fshift_w1) + length(lshift_bshift_ushift_w1) ...
                     + length(lshift_bshift_dshift_w1) + length(lshift_fshift_ushift_w1) + length(lshift_fshift_dshift_w1) ...
                     + length(rshift_w1) + length(rshift_ushift_w1) + length(rshift_dshift_w1) ...
                     + length(rshift_bshift_w1) + length(rshift_fshift_w1) + length(rshift_bshift_ushift_w1) ...
                     + length(rshift_bshift_dshift_w1) + length(rshift_fshift_ushift_w1) + length(rshift_fshift_dshift_w1) ...
                     + length(ushift_w1) + length(ushift_bshift_w1) + length(ushift_fshift_w1) ...
                     + length(dshift_w1) + length(dshift_bshift_w1) + length(dshift_fshift_w1) ...
                     + length(bshift_w1) + length(bshift_w2);

            
            % Create structure for storing parameters
            bdry_params = struct('length', len, ...
                                 'lshift', struct('edges', lshift, 'sig_edges', lshift_sig, 'w1', lshift_w1,'w2', lshift_w2), ...
                                 'lshift_ushift', struct('edges', lshift_ushift, 'sig_edges', lshift_ushift_sig, 'w1', lshift_ushift_w1,'w2', lshift_ushift_w2), ...
                                 'lshift_dshift', struct('edges', lshift_dshift, 'sig_edges', lshift_dshift_sig, 'w1', lshift_dshift_w1,'w2', lshift_dshift_w2), ...
                                 'lshift_bshift', struct('edges', lshift_bshift, 'sig_edges', lshift_bshift_sig, 'w1', lshift_bshift_w1,'w2', lshift_bshift_w2), ...
                                 'lshift_fshift', struct('edges', lshift_fshift, 'sig_edges', lshift_fshift_sig, 'w1', lshift_fshift_w1,'w2', lshift_fshift_w2), ...
                                 'lshift_bshift_ushift', struct('edges', lshift_bshift_ushift, 'sig_edges', lshift_bshift_ushift_sig, 'w1', lshift_bshift_ushift_w1,'w2', lshift_bshift_ushift_w2), ...
                                 'lshift_bshift_dshift', struct('edges', lshift_bshift_dshift, 'sig_edges', lshift_bshift_dshift_sig, 'w1', lshift_bshift_dshift_w1,'w2', lshift_bshift_dshift_w2), ...
                                 'lshift_fshift_ushift', struct('edges', lshift_fshift_ushift, 'sig_edges', lshift_fshift_ushift_sig, 'w1', lshift_fshift_ushift_w1,'w2', lshift_fshift_ushift_w2), ...
                                 'lshift_fshift_dshift', struct('edges', lshift_fshift_dshift, 'sig_edges', lshift_fshift_dshift_sig, 'w1', lshift_fshift_dshift_w1,'w2', lshift_fshift_dshift_w2), ...
                                 'rshift', struct('edges', rshift, 'sig_edges', rshift_sig, 'w1', rshift_w1,'w2', rshift_w2), ...
                                 'rshift_ushift', struct('edges', rshift_ushift, 'sig_edges', rshift_ushift_sig, 'w1', rshift_ushift_w1,'w2', rshift_ushift_w2), ...
                                 'rshift_dshift', struct('edges', rshift_dshift, 'sig_edges', rshift_dshift_sig, 'w1', rshift_dshift_w1,'w2', rshift_dshift_w2), ...
                                 'rshift_bshift', struct('edges', rshift_bshift, 'sig_edges', rshift_bshift_sig, 'w1', rshift_bshift_w1,'w2', rshift_bshift_w2), ...
                                 'rshift_fshift', struct('edges', rshift_fshift, 'sig_edges', rshift_fshift_sig, 'w1', rshift_fshift_w1,'w2', rshift_fshift_w2), ...
                                 'rshift_bshift_ushift', struct('edges', rshift_bshift_ushift, 'sig_edges', rshift_bshift_ushift_sig, 'w1', rshift_bshift_ushift_w1,'w2', rshift_bshift_ushift_w2), ...
                                 'rshift_bshift_dshift', struct('edges', rshift_bshift_dshift, 'sig_edges', rshift_bshift_dshift_sig, 'w1', rshift_bshift_dshift_w1,'w2', rshift_bshift_dshift_w2), ...
                                 'rshift_fshift_ushift', struct('edges', rshift_fshift_ushift, 'sig_edges', rshift_fshift_ushift_sig, 'w1', rshift_fshift_ushift_w1,'w2', rshift_fshift_ushift_w2), ...
                                 'rshift_fshift_dshift', struct('edges', rshift_fshift_dshift, 'sig_edges', rshift_fshift_dshift_sig, 'w1', rshift_fshift_dshift_w1,'w2', rshift_fshift_dshift_w2), ...
                                 'ushift', struct('edges', ushift, 'sig_edges', ushift_sig, 'w1', ushift_w1,'w2', ushift_w2), ...
                                 'ushift_bshift', struct('edges', ushift_bshift, 'sig_edges', ushift_bshift_sig, 'w1', ushift_bshift_w1,'w2', ushift_bshift_w2), ...
                                 'ushift_fshift', struct('edges', ushift_fshift, 'sig_edges', ushift_fshift_sig, 'w1', ushift_fshift_w1,'w2', ushift_fshift_w2), ...
                                 'dshift', struct('edges', dshift, 'sig_edges', dshift_sig, 'w1', dshift_w1,'w2', dshift_w2), ...
                                 'dshift_bshift', struct('edges', dshift_bshift, 'sig_edges', dshift_bshift_sig, 'w1', dshift_bshift_w1,'w2', dshift_bshift_w2), ...
                                 'dshift_fshift', struct('edges', dshift_fshift, 'sig_edges', dshift_fshift_sig, 'w1', dshift_fshift_w1,'w2', dshift_fshift_w2), ...
                                 'bshift', struct('edges', bshift, 'sig_edges', bshift_sig, 'w1', bshift_w1,'w2', bshift_w2), ...
                                 'fshift', struct('edges', fshift, 'sig_edges', fshift_sig, 'w1', fshift_w1,'w2', fshift_w2)); 
        end
end

