% storeexists.m: Determine whether a store is present in the current header.
%   Spencer Torene, BU 2013-10-16
function SE = storeexists(DT,req_tag)

ft = find(strcmp({DT.contents.tag},req_tag),1);

SE = ~isempty(ft);