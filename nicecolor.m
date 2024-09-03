function cout = nicecolor(c)

% cvector = nicecolor('r' | 'b' | etc. | 'R' | 'B' | etc. | 'yr' | 'brW' | etc.);
% cvector = nicecolor('q1' ... 'q9');
% cvector = nicecolor(same cvector)
% clist = nicecolor(cell array of the above);
%
% NICECOLOR(a letter) returns the numbers for a standard Matlab color ('r','k',
% etc.) or nicer versions of the same ('R','B',etc.). Blends are allowable too:
% nicecolor('Bkr') is the element-by-element average of 'B' and 'k' and 'r'.
% NICECOLOR('q' followed by a number) is a shade of gray: 'q7' is [.7 .7 .7].
% NICECOLOR('q') = NICECOLOR('q5').
% if a 3-element color vector is passed to NICECOLOR, it passes out the other
% end unchanged.
%
% neil banas, feb 2002
% (neil@ocean.washington.edu)

%-------------------
cnames = 'rgbcmykwRGBCMYKW';
colors = [1  0  0; % r
0  1  0; % g
0  0  1; % b
0  1  1; % c
1  0  1; % m
1  1  0; % y
0  0  0; % k
1  1  1; % w
1 .4 .4; % R: salmony
0 .7  0; % G: a bit darker
0 .4  1; % B: lighter: prints like 'b' appears on screen
.2  1  1; % C: a bit darker
.8  0 .6; % M: purple
.9 .8  0; % Y: a bit darker
0  0  0; % K
1  1  1];% W
%--------------------

cc={};
if ~iscell(c)
	for r = 1:size(c,1), cc = {cc{:} c(r,:)}; end	
else
	cc = c;
end
for r = 1:length(cc)
	if ischar(cc{r})
		cout(r,:) = colorblend(cc{r},colors,cnames);
	else
		cout(r,:) = cc{r};
	end
end


%--------------------
function cout0 = colorblend(c0,colors,cnames)
cl = [];
i = 1;
while i <= length(c0)
	doublelength = 0;
	if c0(i)=='q' & i~=length(c0)
		doublelength = (c0(i+1) >= '0' & c0(i+1) <= '9');
	end
	if doublelength
		cl = [cl; singlecolor(c0(i:i+1),colors,cnames)];
		i = i + 2;
	else
		cl = [cl; singlecolor(c0(i),colors,cnames)];
		i = i + 1;
	end
end
if size(cl,1) > 1
	cout0 = mean(cl);
else
	cout0 = cl;
end


%--------------------
function cout1 = singlecolor(c1,colors,cnames);

if c1(1)=='q'
	if length(c1)==1
		cout1 = [.5 .5 .5];
	else
		cout1 = str2num(c1(2))/10 .* [1 1 1];
	end
elseif length(c1) > 1
	error(['nicecolor.m: bad color name (''' c1 ''')']);
else
	j = find(cnames==c1);
	if isempty(j), error(['nicecolor.m: bad color name (''' c1 ''')']); end
	cout1 = colors(find(cnames==c1),:);
end



