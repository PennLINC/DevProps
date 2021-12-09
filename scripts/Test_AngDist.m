function Test_AngDist(subj)

% evaluate tru angular distance between PG and Op Fl direction

TrueAngDist(subj)

% evaluate permuted angular distance between PG g and Op Fl directions
for i=1:1000
	PermAngDist(subj,i)
end
