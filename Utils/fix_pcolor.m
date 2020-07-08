function [new_horz,new_vert,new_C] = fix_pcolor(horz,vert,C)
% Duplicate the last column and row in matrix C so that pcolor plots all
% data. horz and vert are vectors of the horizontal and vertical positions,
% respectively. Last column and row are just set to the same width of the
% duplicated column and row.

shift =(0.5.*diff(horz));
test = horz(1:end-1) + shift;

if iscolumn(horz)    
%     horz = [horz; horz(end)+horz(end)-horz(end-1)];
    test = [horz(1)-shift(1); test];
else
%     horz = [horz horz(end)+horz(end)-horz(end-1)];
test = [horz(1)-shift(1) test];
end
test(end+1) = horz(end) + shift(end);

if iscolumn(vert)
    vert = [vert; vert(end)+vert(end)-vert(end-1)];
else
    vert = [vert vert(end)+vert(end)-vert(end-1)];
end

C = horzcat(C,C(:,end));

C = vertcat(C,C(end,:));

new_horz = test;
new_vert = vert;
new_C = C;

end