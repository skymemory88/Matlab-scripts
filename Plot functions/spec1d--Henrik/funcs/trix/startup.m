% Matlab libraries ==================================================

% This line indicates the path to Matlab libraries 19/11/98
%libroot = [ matlabroot '/toolbox/local/' ];
libroot = '/home/tas/matroot/matlab';

path(path,libroot);
libroot = [ libroot '/' ];

%update startup.m for next session in case of modification by root
copyfile([ libroot 'startup.m' ],[ pwd '/startup.m' ]);

disp('Install local libs : SpectralTools, FileTools, MiscTools, funcs...');
path(path,[ libroot 'FileTools' ]);
path(path,[ libroot 'MiscTools' ]);
path(path,[ libroot 'SpectralTools' ]);
%path(path,[ libroot 'funcs' ]);

fprintf(1,'Install Matlab Graphic apps : ');
fprintf(1,'Mfit, ');
path(path,[ libroot 'mfit4' ]);
path(path,[ libroot 'load' ]);
path(path,[ libroot 'nllsq' ]);
path(path,[ libroot 'funcs' ]);

fprintf(1,'Rescal, ');
path(path,[ libroot 'rescal5' ]);
path(path,[ libroot 'funcs/trix' ]);

fprintf(1,'Mview, ');
path(path,[ libroot 'mview4' ]);

fprintf(1,'SSym, ');
path(path,[ libroot 'ssym1' ]);

fprintf(1,'Spec1d packages. \n');

fprintf(1,'Install temporary packages : ');

fprintf(1,'Riso 98.\n');

path(path,[ '/home/tas/wildes/matlab/fit' ]);


% Matlab auto-install ===========================================

st_user = getenv('USER');
st_path = getenv('PATH');
[st_a,st_b] = unix('looktxt > /dev/null');

if st_a == 1	% looktxt hasn't been compiled yet...
	disp('WARN : can''t find "looktxt" program (Needed for some load routines).')
	disp('Put it in your bin directory and add the following line at the end of your .login file :')
	disp([ 'set path = ($path ~/bin)' ])
end
	
if exist('ffind') ~= 3
	disp('WARN : can''t find "ffind" mex file (Needed for some load routines).')
end

fprintf(1,'Hello %s ! Your PRINTER is ',st_user);

if isempty(getenv('PRINTER'))
	fprintf(1,'(undefined).\n');
	if strcmp(st_user,'in1') | strcmp(st_user,'in8') | strcmp(st_user,'in14') | strcmp(st_user,'in20')
	disp('WARN : Add the following line at the end of your .cshrc file :')
	disp([ 'setenv PRINTER hp1_' st_user '_ps' ]);
end
else
	fprintf(1,'%s.\n',getenv('PRINTER'))
end



%disp('use "doc" command to access Hypertext OnLine Help (with Netscape).');
disp('TIP : Ctrl-a = go to start of line, Ctrl-e = go to end of line')
disp('      Up ar. = backward history completion. Type "help cedit" for more.')

colordef none

format short e
eval([ 'cd ' getenv('HOME') ]);
if exist('startuser')
	disp('Now reading your $HOME/startuser.m file ...');
	startuser
else
	disp('OK, you don"t have any $HOME/startuser.m file.');
end

fprintf(1,'Current directory is : %s\n',pwd);

set(0,'DefaultFigurePaperUnits','centimeters');
set(0,'DefaultFigurePaperType','A4');

clear libroot st_a st_b st_user st_path
