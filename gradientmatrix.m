function [nabla,L] =gradientmatrix(R,U)
nabla=R.'/(R*(R.'));
% disp('nabla=');
% disp(nabla);
L=(U*nabla).';
% disp('L=');
% disp(L);

