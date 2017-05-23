% -------------------------------------------------------------------------
% Load two asteroid populations from .mat files, with the names
% of the .mat files given as strings in the inputs 'population01'
% and 'population02' (both .mat files must contain a directory
% listing cell variable 'lcg_files').
%
% When the populations have been loaded, the function reconstructs
% their p and beta distributions and performs a comparison.
%
% Outputs:
% - dp: the statistical differences for the CDFs of p (L1, L2, L_inf)
% - dB: the statistical differences for the CDFs of beta (as above)
% - results: DFs and CDFs so that they can be plotted
% -------------------------------------------------------------------------

function [dp,dB,results]=ast_comparison_WISE(population01, population02)

close all
drawnow

% Reconstructing population 1
load(population01)
[P1, BETA1, Pmargin1, Bmargin1] = ast_comparison_loop(lcg_files);
save('sample01', 'P1', 'BETA1', 'Pmargin1', 'Bmargin1')

% Reconstructing population 2
load(population02)
[P2, BETA2, Pmargin2, Bmargin2] = ast_comparison_loop(lcg_files);
save('sample02', 'P2', 'BETA2', 'Pmargin2', 'Bmargin2')

% Comparing populations
if ( nargout == 2 )
    [dp,dB]=KS_comparison('sample01', 'sample02');
elseif ( nargout == 3 )
    [dp,dB,results]=KS_comparison('sample01', 'sample02');
else
    KS_comparison('sample01', 'sample02')
end

end



function [P, BETA, Pmargin, Bmargin] = ast_comparison_loop(database)

lcg_files = database;
Nkierroksia=length(lcg_files);
leader_main_WISE

end