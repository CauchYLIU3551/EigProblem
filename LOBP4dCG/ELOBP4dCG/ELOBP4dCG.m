## Copyright (C) 2021 root
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} ELOBP4dCG (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-04-23

function [eigenvalue, eigenvector] = ELOBP4dCG (K, M, Z0, E_plus, E_minus, nb=1, m=3, tol=10^-3, opts=1)
%% This function computes the smallest k eigenvalues and corresponding eigenvec,
%% by ELOBP4dCG method.

%% Input:
%       K      array (n-by-n), SPD
%       M      array (n-by-n), SPD
%       E_plus      array (n-by-n), SPD
%       E_minus      array (n-by-n), SPD
%       Z0     array (2n-by-nb) whose columns span an approximate 
%              invariant subspace associated with (k0+1)st to (k0+k)th 
%              smallest positive eigenvalues
%       k,     int , the number of desired eigenvalues
%       m      int >=2, the order of krylov subspace default set as 3
%       tol    tolerence between true value and numerical result
%       opts   precondition options.

if opts==1
  n_cvgd=0;
  A_cvgd=[];
  Z_cvgd=[];
end

%initial approximate
X=Z0(n+1:2*n,:); 
Y=Z0(1:n,:);  

XKX=X'*KX; colXKX=chol(XKX); invXKX=inv(colXKX); X=X*invXKX; KX=KX*invXKX; % X'*shiftedK*X=I
MY=M*Y;
YMY=Y'*MY; colYMY=chol(YMY); invYMY=inv(colYMY); Y=Y*invYMY; MY=MY*invYMY; % Y'*M*Y=I

W=X'*Y; W=0.5*(W+W');
RX=KX*W-Y; RY=MY*W-X;

for i=1:nb
  RY(:,i)=CG(M,RY(:,i));
  RX(:,i)=CG(K,RX(:,i));
end

RY=[Y,RY];
RX=[X,RX];
RY=M_GS(RY,M);
RX=M_GS(RX,K);

[V,D,U]=svd(RY'*E_minus*RX);

[s,idx]=sort(diag(S1),'descend'); 
ns=length(s); %s=s(1:nb);
endfunction
