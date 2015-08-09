% test script for local interpolation of the corrlation matrix

val = 3+randn(10 ,10);
val(t1) = -1;

numIdx = length(t1);
% find consequtive indexs in array
ind_pad = [0 ;t1; inf];
a=diff(ind_pad);
b = find( [a ;inf]>1 );
b_tag = [a; inf] > 1;
b_tag = b_tag & circshift(b_tag, -1); % all non conseq fv indexs
c = diff([0 ;b]); % all non 1 are sequence
c_idx = c ~= 1; % location with the length of the sequence
d=cumsum(c)-1; % endpoints of the sequences corresponding to c
tf = val;
tf(t1(b_tag(1:numIdx))) = mean([tf(t1(b_tag(1:numIdx))-1); tf(t1(b_tag(1:numIdx))+1)]); % means non conseq fv indexs
for  k = find(c_idx)' % means conseq fv indexs
    tf(t1(d(k))+1-c(k):t1(d(k))) = mean(val([t1(d(k))-c(k), t1(d(k))+1]));
end
% d_p = [1 ;d(1:end-1)+1];
% IND = d-d_p ~= 0;
