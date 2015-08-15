function [ind] = renumbering(ind1)
% Renumbering neighbours

o1 = ind1(2);
o2 = ind1(4);
o3 = ind1(6);
o4 = ind1(5);
o5 = ind1(3);
o6 = ind1(1);
ind = [o1 o2 o3 o4 o5 o6];

end
