function [dobs, dpred] = get_station_loc_2D(s,dobs,dpred,m)
% get locations of stations in km
% use column locations from the par file, and ksurf values from the model
% file
dobs.y = m.cy(s.allsites)./1000;
dobs.z = m.z(m.z_ind(s.allsites))./1000;

dpred.y = m.cy(s.allsites)./1000;
dpred.z = m.z(m.z_ind(s.allsites))./1000;

end
