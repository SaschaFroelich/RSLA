function [ntrans,trans]=transitive(x)
% Input:
%       x:      Matrix
% Output:
%       ntrans: number of non-transitive triples in matrix x.
%       trans:  number of transitive triples in matrix x.
%
% May contain NaN-entries. This script is only useful for TD-matrices that
% are skew-symmetric (i.e. they were computed using xcorr).
[rows,columns]=size(x);

if rows~=columns
   msgbox('Matrix must be square for estimation of transitivity.','Error','error'); 
end

if length(find(isnan(x)))
   for entry=1:rows*columns
      if isnan(x(entry))
         x(entry)=0; 
      end
   end
end


if ~issymmetric(x,'skew')
    msgbox('Matrix must be skew-symmetric for transitivity-estimation.','error','error');
end

nodes = rows;
    
%Compute no. of all possible triples
triples=nchoosek(nodes,3);

trans = 0; %no. of trans triples
ntrans = 0; %no. of non-trans triples

%run over all triples
for i=1:nodes
	for j=i+1:nodes
        for k=j+1:nodes
            scoreseq = [0 0 0];
                %stores the number of edges coming into nodes as in scoreseq=[i_in,j_in,k_in]. In a trans triple this has to be 0,1,2 (or a permutation thereof). 
                %For example, 2_in increases by one if x(1,2)<0 because this means that
              	%node 2 fires after node 1
                %disp(['Check if ' num2str(x(i,j)) ', ' num2str(x(i,k)) ',' num2str(x(k,j)) ' is a trans triple']);
                
                %If a triple contains a single NaN-entry, the triple must
                %be considered non-trans
                if isnan(x(i,j)) || isnan(x(i,k)) || isnan(x(j,k))
                    ntrans = ntrans+1;
                else
                    %(i,j)
                    if x(i,j) < 0
                        scoreseq(2) = scoreseq(2) +1;
                    elseif x(i,j) > 0
                        scoreseq(1) = scoreseq(1) + 1;
                    end

                    %(i,k)
                    if x(i,k) < 0
                        scoreseq(3) = scoreseq(3) +1;
                    elseif x(i,k) > 0
                        scoreseq(1) = scoreseq(1) + 1;
                    end

                    %(j,k)
                    if x(j,k) < 0
                        scoreseq(3) = scoreseq(3) +1;
                    elseif x(j,k) > 0
                        scoreseq(2) = scoreseq(2) + 1;
                    end
                    if sortrows(scoreseq') == [0 1 2]'
                        %disp(['The triple (' num2str(i) ',' num2str(j) '),(' num2str(i) ',' num2str(k) '),(' num2str(j) ',' num2str(k) ') is trans!']);
                        trans = trans +1;
                    else 
                        %disp(['The triple (' num2str(i) ',' num2str(j) '),(' num2str(i) ',' num2str(k) '),(' num2str(j) ',' num2str(k) ') is NOT trans!']);
                        ntrans = ntrans +1;
                    end
                end
                
        end
	end
end



end