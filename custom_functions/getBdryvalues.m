function [bdry_values] = getBdryvalues(field, bdry_params, connectivity)
% Linear interpolation of a field to a boundary.
% Input:
%  field:       random field over a domain in R^D, it is an (D+1)-dimensional array,
%               where the last dimension enumerates the samples
%  bdry_params: structure generated from getBdryparams
%  connectivity: A binary value, use 0 if bdry_values have been obtained in 6cc, 1 if 26cc
%
% Output:
%  bdry_values are the linear interpolated values of field along the
%  boundary.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Alex Bowring
% Last changes: 10/25/2018
%__________________________________________________________________________

%%%%% Compute parameters of the field input
dim = size(field);
D   = length(dim);

%%%%% Treat cases of 1D and 2D domain of the field differently
switch D
    case 2
        %%%%% Linear Interpolation onto the boundary
        lshift_bdry_values = bdry_params.lshift.w1.*field(bdry_params.lshift.edges) + bdry_params.lshift.w2.*field(bdry_params.lshift.sig_edges);
        rshift_bdry_values = bdry_params.rshift.w1.*field(bdry_params.rshift.edges) + bdry_params.rshift.w2.*field(bdry_params.rshift.sig_edges);
        ushift_bdry_values = bdry_params.ushift.w1.*field(bdry_params.ushift.edges) + bdry_params.ushift.w2.*field(bdry_params.ushift.sig_edges);
        dshift_bdry_values = bdry_params.dshift.w1.*field(bdry_params.dshift.edges) + bdry_params.dshift.w2.*field(bdry_params.dshift.sig_edges);

        bdry_values = [lshift_bdry_values; rshift_bdry_values; ushift_bdry_values; dshift_bdry_values];

    case 3
        switch connectivity
            case 0 
            %%%%% Linear Interpolation onto the boundary
            lshift_bdry_values = bdry_params.lshift.w1.*field(bdry_params.lshift.edges) + bdry_params.lshift.w2.*field(bdry_params.lshift.sig_edges);
            rshift_bdry_values = bdry_params.rshift.w1.*field(bdry_params.rshift.edges) + bdry_params.rshift.w2.*field(bdry_params.rshift.sig_edges);
            ushift_bdry_values = bdry_params.ushift.w1.*field(bdry_params.ushift.edges) + bdry_params.ushift.w2.*field(bdry_params.ushift.sig_edges);
            dshift_bdry_values = bdry_params.dshift.w1.*field(bdry_params.dshift.edges) + bdry_params.dshift.w2.*field(bdry_params.dshift.sig_edges);
            bshift_bdry_values = bdry_params.bshift.w1.*field(bdry_params.bshift.edges) + bdry_params.bshift.w2.*field(bdry_params.bshift.sig_edges);
            fshift_bdry_values = bdry_params.fshift.w1.*field(bdry_params.fshift.edges) + bdry_params.fshift.w2.*field(bdry_params.fshift.sig_edges);
            
            bdry_values = [lshift_bdry_values; rshift_bdry_values; ushift_bdry_values; dshift_bdry_values; bshift_bdry_values; fshift_bdry_values];

            case 1
            %%%%% Linear Interpolation onto the boundary
            lshift_bdry_values        = bdry_params.lshift.w1.*field(bdry_params.lshift.edges) + bdry_params.lshift.w2.*field(bdry_params.lshift.sig_edges);
            lshift_ushift_bdry_values = bdry_params.lshift_ushift.w1.*field(bdry_params.lshift_ushift.edges) + bdry_params.lshift_ushift.w2.*field(bdry_params.lshift_ushift.sig_edges);
            lshift_dshift_bdry_values = bdry_params.lshift_dshift.w1.*field(bdry_params.lshift_dshift.edges) + bdry_params.lshift_dshift.w2.*field(bdry_params.lshift_dshift.sig_edges);
            lshift_bshift_bdry_values = bdry_params.lshift_bshift.w1.*field(bdry_params.lshift_bshift.edges) + bdry_params.lshift_bshift.w2.*field(bdry_params.lshift_bshift.sig_edges);
            lshift_fshift_bdry_values = bdry_params.lshift_fshift.w1.*field(bdry_params.lshift_fshift.edges) + bdry_params.lshift_fshift.w2.*field(bdry_params.lshift_fshift.sig_edges);
            lshift_bshift_ushift_bdry_values = bdry_params.lshift_bshift_ushift.w1.*field(bdry_params.lshift_bshift_ushift.edges) + bdry_params.lshift_bshift_ushift.w2.*field(bdry_params.lshift_bshift_ushift.sig_edges);
            lshift_bshift_dshift_bdry_values = bdry_params.lshift_bshift_dshift.w1.*field(bdry_params.lshift_bshift_dshift.edges) + bdry_params.lshift_bshift_dshift.w2.*field(bdry_params.lshift_bshift_dshift.sig_edges);
            lshift_fshift_ushift_bdry_values = bdry_params.lshift_fshift_ushift.w1.*field(bdry_params.lshift_fshift_ushift.edges) + bdry_params.lshift_fshift_ushift.w2.*field(bdry_params.lshift_fshift_ushift.sig_edges);
            lshift_fshift_dshift_bdry_values = bdry_params.lshift_fshift_dshift.w1.*field(bdry_params.lshift_fshift_dshift.edges) + bdry_params.lshift_fshift_dshift.w2.*field(bdry_params.lshift_fshift_dshift.sig_edges);
            rshift_bdry_values        = bdry_params.rshift.w1.*field(bdry_params.rshift.edges) + bdry_params.rshift.w2.*field(bdry_params.rshift.sig_edges);
            rshift_ushift_bdry_values = bdry_params.rshift_ushift.w1.*field(bdry_params.rshift_ushift.edges) + bdry_params.rshift_ushift.w2.*field(bdry_params.rshift_ushift.sig_edges);
            rshift_dshift_bdry_values = bdry_params.rshift_dshift.w1.*field(bdry_params.rshift_dshift.edges) + bdry_params.rshift_dshift.w2.*field(bdry_params.rshift_dshift.sig_edges);
            rshift_bshift_bdry_values = bdry_params.rshift_bshift.w1.*field(bdry_params.rshift_bshift.edges) + bdry_params.rshift_bshift.w2.*field(bdry_params.rshift_bshift.sig_edges);
            rshift_fshift_bdry_values = bdry_params.rshift_fshift.w1.*field(bdry_params.rshift_fshift.edges) + bdry_params.rshift_fshift.w2.*field(bdry_params.rshift_fshift.sig_edges);
            rshift_bshift_ushift_bdry_values = bdry_params.rshift_bshift_ushift.w1.*field(bdry_params.rshift_bshift_ushift.edges) + bdry_params.rshift_bshift_ushift.w2.*field(bdry_params.rshift_bshift_ushift.sig_edges);
            rshift_bshift_dshift_bdry_values = bdry_params.rshift_bshift_dshift.w1.*field(bdry_params.rshift_bshift_dshift.edges) + bdry_params.rshift_bshift_dshift.w2.*field(bdry_params.rshift_bshift_dshift.sig_edges);
            rshift_fshift_ushift_bdry_values = bdry_params.rshift_fshift_ushift.w1.*field(bdry_params.rshift_fshift_ushift.edges) + bdry_params.rshift_fshift_ushift.w2.*field(bdry_params.rshift_fshift_ushift.sig_edges);
            rshift_fshift_dshift_bdry_values = bdry_params.rshift_fshift_dshift.w1.*field(bdry_params.rshift_fshift_dshift.edges) + bdry_params.rshift_fshift_dshift.w2.*field(bdry_params.rshift_fshift_dshift.sig_edges);
            ushift_bdry_values = bdry_params.ushift.w1.*field(bdry_params.ushift.edges) + bdry_params.ushift.w2.*field(bdry_params.ushift.sig_edges);
            ushift_bshift_bdry_values = bdry_params.ushift_bshift.w1.*field(bdry_params.ushift_bshift.edges) + bdry_params.ushift_bshift.w2.*field(bdry_params.ushift_bshift.sig_edges);
            ushift_fshift_bdry_values = bdry_params.ushift_fshift.w1.*field(bdry_params.ushift_fshift.edges) + bdry_params.ushift_fshift.w2.*field(bdry_params.ushift_fshift.sig_edges);
            dshift_bdry_values = bdry_params.dshift.w1.*field(bdry_params.dshift.edges) + bdry_params.dshift.w2.*field(bdry_params.dshift.sig_edges);
            dshift_bshift_bdry_values = bdry_params.dshift_bshift.w1.*field(bdry_params.dshift_bshift.edges) + bdry_params.dshift_bshift.w2.*field(bdry_params.dshift_bshift.sig_edges);
            dshift_fshift_bdry_values = bdry_params.dshift_fshift.w1.*field(bdry_params.dshift_fshift.edges) + bdry_params.dshift_fshift.w2.*field(bdry_params.dshift_fshift.sig_edges);
            bshift_bdry_values = bdry_params.bshift.w1.*field(bdry_params.bshift.edges) + bdry_params.bshift.w2.*field(bdry_params.bshift.sig_edges);
            fshift_bdry_values = bdry_params.fshift.w1.*field(bdry_params.fshift.edges) + bdry_params.fshift.w2.*field(bdry_params.fshift.sig_edges);

            bdry_values = [lshift_bdry_values; lshift_ushift_bdry_values; lshift_dshift_bdry_values; lshift_bshift_bdry_values; lshift_fshift_bdry_values; ...
                           lshift_bshift_ushift_bdry_values; lshift_bshift_dshift_bdry_values; lshift_fshift_ushift_bdry_values; lshift_fshift_dshift_bdry_values; ...
                           rshift_bdry_values; rshift_ushift_bdry_values; rshift_dshift_bdry_values; rshift_bshift_bdry_values; rshift_fshift_bdry_values; ...
                           rshift_bshift_ushift_bdry_values; rshift_bshift_dshift_bdry_values; rshift_fshift_ushift_bdry_values; rshift_fshift_dshift_bdry_values; ...
                           ushift_bdry_values; ushift_bshift_bdry_values; ushift_fshift_bdry_values; dshift_bdry_values; dshift_bshift_bdry_values; dshift_fshift_bdry_values; ...
                           bshift_bdry_values; fshift_bdry_values];

        end
end 
   
     