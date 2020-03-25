/* --------------------- Programme looktxt.c ---------------------- */
 
/* Analyse d'un fichier texte pour rechercher les champs 
 * numeriques sous forme de reels, vecteurs, matrices.
 * Prise en compte des lignes de commentaires (%,#,c...\n) 
 */

/* compile with : gcc looktxt.c -o looktxt */

/* For PC version change '/' into '\\' in lines :
 *
 * if (strrchr(filename, '/') != NULL) filename = strrchr(filename, '/')+1;
 * becomes 
 * if (strrchr(filename, '\\') != NULL) filename = strrchr(filename, '\\')+1;
 *
 * if ( (filename[i] > 122) || ischr(filename[i], "!\"#$%&'()*+,-.:;<=>?@[\\]^_`") ) filename[i] = '_';
 * becomes
 * if ( (filename[i] > 122) || ischr(filename[i], "!\"#$%&'()*+,-.:;<=>?@[/]^_`") ) filename[i] = '_';
 */

#define auteur  "Farhi E. 01/97"
#define date    "04/03/99"
#define version "0.87"
/*
 * content: C language
 * tab = 2 chars
 */


/* Usage :
 * ARGS : filename [ -p={.|,} ]
                   [ -s={\t|\v|,|;} ]
                   [ -c={#|%|c} ]
                   [ -n={\n\r\f} ]
                   [ -o={oct|m|txt} ]
                   [ -v ]
                   [ -g ]
                   [ -r ]
                   [ -a ]
                   [ -h ]
                   [ -l=X ]
                   [ -F=<output> ]
 * OUT : matrix of fields sent to stdout (printed to terminal)
         can be used after to extract particular data.

 * Creates an '.m', '.oct' or simple numeric output file

 * Example :
             looktxt foo.txt -p="." -s='\t\v,;' -c='#'

             'looktxt -h' to get help.
 */

/* remove comments if your compiler does not auto include some libraries */
#include <stdio.h>
/* #include <stdlib.h> */
#include <string.h>
/* #include <ctype.h> */
#include <sys/types.h>
#include <sys/stat.h>


/* -- Declaration section ----------------------------------------- */


#define Bnumber      1
#define Balpha       2
#define Bpoint       4
#define Beol         8
#define Bexp        16
#define Bsign       32
#define Bcomment    64
#define Bseparator 128

#define OUT_OCT      1
#define OUT_M        2
#define OUT_TXT      3

#define filetable "lktmp000.txt"

/* number format : [+-][0-9][.][0-9][eE][+-][0-9] */

/* event types                           {0,a,.,l,e,-,c,s} */
const long     all         =  0xFF;   /* {1,1,1,1,1,1,1,1}; */
const long     none        =  0x00;   /* {0,0,0,0,0,0,0,0}; */

/* events for num search                 {n,a,p,l,e,s,c,s} */
const long needforstartnum =  Bnumber + Bpoint + Bsign;
const long needafternumber =  Bnumber + Bpoint + Bexp + Beol + Bseparator;
const long needaftersign   =  Bnumber + Bpoint;
const long needafterexp    =  Bnumber + Bsign;
const long needafterpoint  =  Bnumber + Bexp   + Beol + Bseparator;
const long needaftereol    =  Bnumber + Bpoint + Bsign + Beol + Bseparator;
const long needaftersep    =  Bnumber + Bpoint + Bsign + Beol + Bseparator; 


char verbose  = 0;
char passascii= 1;
char groupnum = 0;
char noroot   = 0;
char outtype  = OUT_M;
char forcefile= 0;
long line     = 0;
char tablextr = 0;
char username = 0;
char *userfile;
char *filename;

#define number    "0123456789"
#define alpha     " !\"#$&'()*+,-./123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
#define Cpoint     "."
#define Ceol       "\n\f"
#define exp       "eE"
#define sign      "+-"
#define Ccomment   "#%"
#define Cseparator "\t\v\r,; ()[]{}"

char *point, *comment, *separator, *eol;


char *make_esc_char(char *s)
{	/* make escape char such as \[abfnrtv\].... */
	int i=0;	/* in s */
	int j=0;	/* in output */
	char c1 = '\0';
	char c2 = '\0';
	char *sout;

	sout = (char *)malloc(strlen(s)+1);
	do
	{
		c1 = s[i];
		c2 = s[++i];
		if (c1 == '\\')
		{
			switch (c2)
			{
				case 'a' : c1 = '\a'; i++; break;
				case 'b' : c1 = '\b'; i++; break;
				case 'f' : c1 = '\f'; i++; break;
				case 'v' : c1 = '\v'; i++; break;
				case 't' : c1 = '\t'; i++; break;
				case 'n' : c1 = '\n'; i++; break;
				case 'r' : c1 = '\r'; i++; break;
				case '\\': c1 = '\\'; i++; break;
			}
		}
		sout[j++] = c1;
		if (i >255)
		{
			printf("looktxt: make esc char : string too long\n");
			exit(-1);
		}
	}
	while ((c2 != '\0') && (c1 != '\0'));
	sout[j] = '\0';
	return(sout);
}

char *undo_esc_char(char *s)
{	/* undo escape char such as \[abfnrtv\].... */
	int i=0;	/* in s */
	long j=0;
	char c2 = '\0';
	char *sout;

	sout = (char *)malloc(strlen(s)*5);
	sout[0]='\0';
	do
	{
		c2 = s[i++];
		switch (c2)
		{
			case '\a' : sout = (char *)strcat(sout," BEL "); break;
			case '\b' : sout = (char *)strcat(sout," BS "); break;
			case '\f' : sout = (char *)strcat(sout," FF "); break;
			case '\v' : sout = (char *)strcat(sout," VT "); break;
			case '\t' : sout = (char *)strcat(sout," HT "); break;
			case '\n' : sout = (char *)strcat(sout," NL "); break;
			case '\r' : sout = (char *)strcat(sout," CR "); break;
			default   : j=strlen(sout); sout[j] = c2; sout[j+1]='\0'; break;
		}
		if (i > 255)
		{
			printf("looktxt: undo esc char : string too long\n");
			exit(-1);
		}
	}
	while (c2 != '\0');
	return(sout);
}


int ischr(char c, char *category)
{
	return (strchr( category   ,c) != NULL);
}

void printfield(char *filestr, long pos, long filesize, long size)
{
	int i;

	for (i = 0; ((pos+i < filesize) && (i < 20) && (i <= size)); i++)
		if (isprint(filestr[pos+i])) printf("%c",filestr[pos+i]);
		else printf(" ");
}

long scanfile(char *filename)
{
	FILE *fnum;
	FILE *fout;
	FILE *ftmp;
	char c;
	char *fileout;
	char *filestr;
	char *filestrcp;
	char *tablestr;
	struct stat stfile;
	long filesize;

	int  last_is     = 0; /* type of preceeding char */
	long last_eolpos = 0; /* last EOL position */
	long last_seppos = 0; /* Separator */
	int  is          = 0; /* current char type */
	int  found       = 0; /* is that char expected ? */
	long pos         = 0; /* current pos */
	long rows        = 0; /* number of rows in matrix */
	long columns     = 0; /* number of columns in matrix */
	long last_columns= 0; /* to test if number of columns has changed since last line */
	long startcmtpos = 0; /* comment field starting pos */
	long startnumpos = 0; /* num field start pos */
	long startcharpos= 0; /* char field start pos */
	long endcmtpos   = 0; /* comment field end pos */
	long endnumpos   = 0; /* num field end pos */
	long endcharpos  = 0; /* char field end pos */
	int  possiblenum = 0; /* we migth be in a number ... */
	int  possiblecmt = 0; /* we are in a comment */
	int  fieldend    = 0; /* detect end of field */
	long id          = 0; /* field number */
	long i           = 0; /* index */
	int  need        = 0; /* what is to be expected for next char */
	char inpoint     = 0; /* if we are after a point */
	char inexp       = 0; /* if we are after an exp */
	long prod        = 0; /* cols * rows */
	long numextr     = 0;

	long maxfield=0;
	long maxprod=0;


/* file section ----------------------------------------------------- */	

/* open, read then close input file */

	fnum = fopen(filename,"r");
	if (fnum == NULL)
	{
		printf("looktxt: unable to open file %s\n",filename);
		exit(-1);
	}
	stat(filename,&stfile);
	filesize = stfile.st_size;
	filestr = (char *)malloc(filesize+2);
	if (filestr == NULL)
	{
		printf("looktxt: unable to malloc %i bytes for filestr %s\n",filesize,filename);
		exit(-1);
	}
	if (!fread(filestr, filesize,1, fnum))
	{
		printf("looktxt: unable to read file %s\n",filename);
		exit(-1);
	}
	fclose(fnum);
	filestr[filesize] = '\n';
	filestr[filesize+1] = '\0';

/* make a copy of input string, and manage file names */

	filestrcp = (char *)malloc(filesize+2);
	if (filestrcp == NULL)
	{
		printf("looktxt: unable to malloc %i bytes for filestrcp %s\n",filesize,filename);
		exit(-1);
	}
	strcpy(filestrcp, filestr);
	for (i=0; i < strlen(filename); i++)
		if ( (filename[i] > 122) || ischr(filename[i], "!\"#$%&'()*+,-.:;<=>?@[\\]^_`") ) filename[i] = '_';
	fileout = (char *)malloc(strlen(filename)+10);
	if (username)
		strcpy(fileout, userfile);
	else
		strcpy(fileout, filename);

	if (outtype == OUT_OCT) fileout = (char *)strcat(fileout,".oct");
	if (outtype == OUT_M) fileout = (char *)strcat(fileout,".m");
	if (outtype == OUT_TXT) fileout = (char *)strcat(fileout,".txt");

/* open output file and table file */

	if (!stat(fileout,&stfile) && !forcefile)
	{
		printf("looktxt: %s already exists. Use '-f' to force file writing.\n",fileout);
		exit(-1);
	}
	if (verbose) printf("Creating output file %s\n",fileout);
	fout = fopen(fileout,"w+");
	if (fout == NULL)
	{
		printf("looktxt: unable to create file %s\n",fileout);
		exit(-1);
	}
	if (tablextr)
	{
		ftmp = fopen(filetable,"w+");
		if (ftmp == NULL)
		{
			printf("looktxt: unable to create file %s\n",filetable);
			exit(-1);
		}
	}

/* get root name for field names */

	if (strrchr(filename, '/') != NULL) filename = strrchr(filename, '/')+1; /* used for field names */
	if (!noroot && (strlen(filename) > 15)) 
	{	
		sprintf(filename, "var%i", (long)time(NULL));
		printf("Warn : filename too long, using private ID for fields : %s\n",filename);
	}
	if (verbose) 
	{
		printf("Scanning...\n");
		printf("   id type     start        end       rows    columns      line\n");
	}

/* init scan ---- */

rows         = 0;
columns      = 0;
last_columns = 0;
is           = Beol;
startcmtpos  = filesize;
endcmtpos    = -1;


/* scan --------- */

do
{
	last_is = is;

	c=filestr[pos];

	if (pos > filesize) /* end of file reached : real exit */
	{
		last_eolpos = filesize;
		last_seppos = filesize;
		last_is     = Beol;
		possiblecmt = 0;
		need = 0; /* generates end of field : end of line */
		c = '\n';
	}
	if (pos == filesize) /* end of file reached */
	{
		c = '\n';
	}	

	is =   Bnumber    * ischr(c, number   )
	     + Balpha     * ischr(c, alpha    )
	     + Bpoint     * ischr(c, point    )
	     + Beol       * ischr(c, eol      )
	     + Bexp       *(ischr(c, exp      ) && possiblenum) 
	                                        /* must be in a number field */
	     + Bsign      *(ischr(c, sign     ) && (last_is & (Bexp | Bseparator | Beol)))
	                                        /* must be after exponent or not in
	                                           number because starting number */
	     + Bcomment   * ischr(c, comment  )
	                                        /* comment starts if we are waiting for it */
	     + Bseparator * ischr(c, separator);


	if (is & Bseparator) c=' ';

	filestrcp[pos] = c;

	if (!(possiblecmt) && (is & Bcomment))  /* activate comment field */
	{
		possiblecmt = 1;
		startcmtpos = pos;
	}
	
	if (possiblecmt)
	{

		if (is & Beol)                        /* end comment on eol */
		{
			endcmtpos = pos;
			fieldend |= Bcomment;
			possiblecmt = 0;
/*			if (endcharpos < startcharpos)
				endcharpos= startcmtpos - 1;
			if (startcharpos < endcharpos) fieldend |= Balpha;*/ 
		}
		else
		{
			is = last_is;
			filestrcp[pos] = ' ';
		}
	}

	if ( (!possiblecmt) && (!possiblenum) && (is & needforstartnum) && (last_is & (Bseparator | Beol)))    /* activate num search */
	{
		possiblenum = 1;
		startnumpos = pos;
		need = needforstartnum;
		inpoint = 0;
		inexp   = 0;
	}

	if (possiblenum && !(possiblecmt))     /* in num field */
	{
		found = is & need;

/* last column : update when found and (EOL and not groupnum) or (EOL and groupnum and columns (not empty previous num line)) */
/* OK : columns = 0 when EOL */
/* OK : newline : when found after EOL (can be really before -> if columns == 1)*/
/* end of num field : found && columns != last_columns */

		if ((last_is & (Bnumber | Bpoint)) && (is & (Beol | Bseparator)))
		{
			columns++; /* detects num end : one more column */
			if (found && (columns == 1)) rows++;	/* this is a new line starting */
		}

		if (is & Beol)
		{
			if (!groupnum)
				if ((columns != last_columns) && (startnumpos < last_eolpos))
				{	                                  /* change in columns -> end of preceeding num field */
					endnumpos = last_eolpos - 1;
					pos       = last_eolpos;
					is        = Beol;
					need      = needforstartnum;
					fieldend |= Bnumber;
					endcharpos= startnumpos - 1;
					columns   = last_columns;
					if (startcharpos <= endcharpos) fieldend |= Balpha;
				}
				else
				{
					last_columns = columns;
					columns = 0;
				}
			else	/* group num */
			{
				if (!columns && last_columns)
				{
					filestrcp[pos] = ' ';	/* pass on : continue in same field */
				}
				else
				{
				  if (columns && !last_columns)
				  {
					last_columns = columns;	/* first line on matrix */
					columns = 0;
				  }
				  else
				  if (last_columns && (columns != last_columns)  && (startnumpos < last_eolpos))
				  {	                                  /* change in columns -> end of preceeding num field */
					endnumpos = last_eolpos - 1;
					pos       = last_eolpos;
					is        = Beol;
					need      = needforstartnum;
					fieldend |= Bnumber;
					endcharpos= startnumpos - 1;
					columns   = last_columns;
					if (startcharpos <= endcharpos) fieldend |= Balpha;
				  }
				  else
					columns = 0;
				}
			}
		}

		if (!found)
		{
			if (last_is & (Beol | Bseparator)) /* end of num field ok */
			{
				endnumpos  = pos - 2;
				endcharpos = startnumpos - 1;
				if (startcharpos <= endcharpos) fieldend |= Balpha;
				fieldend  |= Bnumber;
			}
			else /* anomalous end of num */
			{
				if (startnumpos >= last_seppos) /* first possible number is not a number */
				{
					columns   = last_columns;
					possiblenum = 0; /* abort and pass */
				}
				else
				if ((columns > 0) && (startnumpos >= last_eolpos)) /* only a line */
				{
					endnumpos = last_seppos - 1;
					pos       = last_seppos;
					is        = Bseparator;
					need      = needforstartnum;
					fieldend |= Bnumber;
					endcharpos = startnumpos - 1;
					if (startcharpos <= endcharpos) fieldend |= Balpha;
				}
				else /* already passed more than one line */
				{
					endnumpos = last_eolpos - 1;
					pos       = last_eolpos;
					is        = Beol;
					need      = needforstartnum;
					fieldend |= Bnumber;
					endcharpos= startnumpos - 1;
					columns   = last_columns;
					if (startcharpos <= endcharpos) fieldend |= Balpha;
				}
			}
		}
		else	                              /* still in num */
		{
			if (is & Bpoint)
			{
				if (inpoint || inexp) 
				{
					need = 0; 
				}
				else
				{
					need = needafterpoint;
					inpoint = 1; 
				}
			}
			else
			if (is & Bsign)        need = needaftersign;
			else
			if (is & Bexp)       { need = needafterexp; inpoint = 0; inexp = 1; }
			else
			if (is & Bseparator) { need = needaftersep; inpoint = 0; inexp = 0; }
			else
			if (is & Bnumber)      need = needafternumber;
			else
			if (is & Beol)       { need = needaftereol; inpoint = 0; inexp = 0; }
			else
				need = needafternumber;
		}
	}

	if (fieldend)
	{
		if ((fieldend & Balpha))
		{
			if (verbose)
			{
				printf("%5i %s %10i %10i %10i %10i : ",id,"chr",startcharpos,endcharpos,rows, columns);
				printfield(filestr, startcharpos, filesize, endcharpos - startcharpos);
				printf("\n");
			}
			if (tablextr)
				fprintf(ftmp, "%5i %5i %10i %10i %10i %10i\n",id,Balpha,startcharpos,endcharpos,rows, columns);
			if ((!passascii) && (outtype == OUT_OCT))
			{
				fprintf(fout,"# name: ");
				if (!noroot) fprintf(fout,"%s_",filename);
				fprintf(fout,"s%i\n",id);
				fprintf(fout,"# type: string\n");
				fprintf(fout,"# length: %i\n", endcharpos - startcharpos + 1);
				fwrite(filestr+startcharpos, endcharpos - startcharpos + 1, 1, fout);
				fprintf(fout,"\n");
			}
			else
			if ((!passascii) && (outtype == OUT_M))
			{
				if (!noroot) fprintf(fout,"%s_",filename);
				fprintf(fout,"s%i = '",id);
				for (i=startcharpos; i <= endcharpos; i++)
				{
					if (filestr[i] == '\'')
						fprintf(fout,"\"");
					else
					{
					if (isprint(filestr[i]))
						fprintf(fout,"%c", filestr[i]);
					else
						fprintf(fout," ");
					}
				}
/*				fwrite(filestr+startcharpos, endcharpos - startcharpos + 1, 1, fout); */
				fprintf(fout,"';\n");
			}
			
			startcharpos = pos;
			fieldend -= Balpha;
			id++;
		}	

		if (fieldend & Bnumber)
		{
			columns = last_columns;
			if (columns && (startnumpos <= endnumpos))
			{
				if (rows <= 0) rows = 1;
				if (verbose)
				{
					printf("%5i %s %10i %10i %10i %10i : ",id,"num",startnumpos,endnumpos,rows, columns);
					printfield(filestr, startnumpos, filesize, endnumpos - startnumpos);
					printf("\n");
				}
				if (tablextr)
					fprintf(ftmp, "%5i %5i %10i %10i %10i %10i\n",id,Bnumber,startcharpos,endcharpos,rows, columns);
				prod = rows*columns;
				if (prod > maxprod) { maxprod = prod; maxfield = id; }
				if ((outtype == OUT_OCT) && (prod >= line))
				{
					numextr++;
					fprintf(fout,"# name: ");
					if (!noroot) fprintf(fout,"%s_",filename);
					fprintf(fout,"n%i\n",id);
					if ((rows == 1) && (columns == 1))
						fprintf(fout,"# type: scalar\n");
					else
					{
						fprintf(fout,"# type: matrix\n");
						fprintf(fout,"# rows: %i\n", rows);
						fprintf(fout,"# columns: %i\n", columns);
					}
					for (i=startnumpos; i <= endnumpos; i++)
						if ( (i < startcmtpos) || (i > endcmtpos) )
							fprintf(fout,"%c", filestrcp[i]);
					fprintf(fout,"\n");
				}
				else
				if ((outtype == OUT_M) && (prod >= line))
				{
					numextr++;
					if (!noroot) fprintf(fout,"%s_",filename);
					if ((rows == 1) && (columns == 1)) 	
						fprintf(fout,"n%i =  ",id);
					else		
						fprintf(fout,"n%i = [ ",id);
					for (i=startnumpos; i <= endnumpos; i++)
						if ( (i < startcmtpos) || (i > endcmtpos) )
							fprintf(fout,"%c", filestrcp[i]);
					if ((rows == 1) && (columns == 1)) 	
						fprintf(fout,";\n");
					else
						fprintf(fout," ];\n");
				}
				else
				if ((outtype == OUT_TXT) && (prod >= line))
				{
					numextr++;
					for (i=startnumpos; i <= endnumpos; i++)
						if ( (i < startcmtpos) || (i > endcmtpos) )
							fprintf(fout,"%c", filestrcp[i]);
					fprintf(fout,"\n");
				}
			}
			possiblenum = 0;
			last_eolpos = endnumpos;
			startcharpos = endnumpos+1;
			fieldend -= Bnumber;
			columns = 0;
			rows = 0;
			last_columns = 0;
			pos = endnumpos+1;
			id++;
			inpoint = 0;
		}

		if (fieldend & Bcomment)
		{
			if (verbose)
			{
				printf("%5i %s %10i %10i %10i %10i : ",id,"cmt",startcmtpos,endcmtpos,rows, columns);
				printfield(filestr, startcmtpos, filesize, endcmtpos - startcmtpos);
				printf("\n");
			}
			if (tablextr)
				fprintf(ftmp, "%5i %5i %10i %10i %10i %10i\n",id,Bcomment,startcharpos,endcharpos,rows, columns);
			possiblecmt = 0;
			fieldend -= Bcomment;
			startcharpos = pos+1;
			id++;
		}
	}

	if (is & Beol)
	{
		last_eolpos = pos;
		last_seppos = pos;
	}
		
	if (is & Bseparator)
		last_seppos = pos;

	pos++;
}
	while (pos <= filesize+1);

/* end scan */

	printf("Importation into %s done. %i fields, %i numerics extracted.\n", fileout,id, numextr);
	if (numextr)
	{
		printf("Main numeric field is : ");
		if (!noroot) printf("%s_",filename);
		printf("n%i (%i elements)\n",maxfield, maxprod);
	}
	fflush(NULL);

	free(filestr); 
	free(fileout);
	free(filestrcp);

	if (tablextr)
	{
		fclose(ftmp); 
		ftmp = fopen(filetable,"r");
		if (ftmp == NULL)
		{
			printf("looktxt: unable to open file %s\n",filetable);
			exit(-1);
		} 
		stat(filetable,&stfile);
		filesize = stfile.st_size;
		tablestr = (char *)malloc(filesize+10); 
		if (tablestr == NULL)
		{
			printf("looktxt: unable to malloc %i bytes for filestr %s\n",filesize,filetable);
			exit(-1);
		}
		if (!fread(tablestr, filesize,1, ftmp))
		{
			printf("looktxt: unable to read file %s\n",filetable);
			exit(-1);
		}
		fclose(ftmp);
		tablestr[filesize] = '\0';

		if (outtype == OUT_M)
		{
			if (!noroot) fprintf(fout,"%s_",filename);
			fprintf(fout,"table = [ ");
		}
		if (outtype == OUT_OCT)
		{	
			fprintf(fout,"# name: ");
			if (!noroot) fprintf(fout,"%s_",filename);
			fprintf(fout,"table\n");
			fprintf(fout,"# type: matrix\n");
			fprintf(fout,"# rows: %i\n", id);
			fprintf(fout,"# columns: %i\n", 6);
		}
		fprintf(fout,"%s", tablestr);
		if (outtype == OUT_M) fprintf(fout," ];\n");
		system("rm " filetable); 
		free(tablestr);
		if (verbose)
		{
			printf("Table of fields is : ");
			if (!noroot) printf("%s_",filename);
			printf("table\n");
		}
	}
	fclose(fout);
	return (id);
}

/* MAIN ====================================================== */

int main(int argc, char *argv[])
{
	int i;
	char *arg;

	if (argc == 1)
	{
		printf("Usage : looktxt [options] filename\n");
		printf("Action: Search and export numerics in a text/ascii file.\n");
		printf("Type looktxt -h for extended help, looktxt -v to get default parameters.\n");
		printf("Example : looktxt foo.txt \n");
		printf(auteur " version " version " (" date ")\n");
		exit(-1);
	}
else
	{
		arg = (char *)malloc(1024);
		point = (char *)malloc(256);
		separator = (char *)malloc(256);
		eol = (char *)malloc(256);
		comment = (char *)malloc(256);

		strcpy(point,Cpoint);
		strcpy(separator,Cseparator);
		strcpy(eol,Ceol);
		strcpy(comment,Ccomment);

		userfile = (char *)malloc(1024);
		filename = (char *)malloc(1024);
		filename[0] = '\0';

		for (i=1; i < argc; i++)
		{
			strcpy(arg,argv[i]);
/*			printf("arg = %s\n",arg); */
			if (strlen(arg) > 255)
			{
				printf("looktxt: arg %i too long\n",i-1);
				exit(-1);
			}
			if (arg[0] == '-')	/* option ? */
			{
				if (strncmp(arg,"-v",2) == 0)
					verbose = 1;
				else
				if (strncmp(arg,"-p",2) == 0)
				{
					strcpy(point,arg+3);
				}
				else
				if (strncmp(arg,"-c",2) == 0)
				{
					strcpy(comment,arg+3);
				}
				else
				if (strncmp(arg,"-s",2) == 0)
				{		
					free(separator);
					separator = (char *)make_esc_char((char *)(arg + 3));		
				}
				else
				if (strncmp(arg,"-n",2) == 0)
				{
					free(eol);
					eol = (char *)make_esc_char((char *)(arg+3));
				}
				else
				if (strncmp(arg,"-F",2) == 0)
				{
					username = 1;
					strcpy(userfile,arg+3);;
				}
				else
				if (strncmp(arg,"-o",2) == 0)
				{
					if (arg[3] == 'o')
						outtype = OUT_OCT;
					else
					if (arg[3] == 'm')
						outtype = OUT_M;
					else
					if (arg[3] == 't')
						outtype = OUT_TXT;
					else
						outtype = OUT_M;
				}
				else
				if (strncmp(arg,"-a",2) == 0)
					passascii = 0;
				else
				if (strncmp(arg,"-g",2) == 0)
					groupnum = 1;
				else
				if (strncmp(arg,"-f",2) == 0)
					forcefile = 1;
				else
				if (strncmp(arg,"-r",2) == 0)
					noroot = 1;
				else
				if (strncmp(arg,"-l",2) == 0)   /* -l=X */
        {
					line = atoi(arg+3);
          if (line < 0) line = 0;
          if (line > 60000) line = 60000; 
        }
				else
				if (strncmp(arg,"-t",2) == 0)
					tablextr = 1;
				if (strncmp(arg,"-h",2) == 0)
				{
					printf("Usage : looktxt [options] filename\n\n");
					printf("Action: Search and export numerics in a text/ascii file.\n");
					printf("        The programs looks into your file some numeric fields (and optionally characters).\n");
					printf("        Each identified numeric field is then named and exported into an output filename.\n");
					printf("        The output file can be a .m (Matlab), .oct (Octave) or .txt file, with automatic or specified user name.\n");
					printf("        For a special type of data, you can specify point, separator, end of line and comment conventions.\n");
					printf("        A 'table' field can be exported, containing informations about all other found fields (see end of this help).\n");
					printf("Examples : \n");
					printf("  'looktxt foo.txt'\n");
					printf("           will create 'foo_txt.m' file with all numerics found in 'foo.txt'. Fields have names such as 'foo_txt_n12')\n");
					printf("  'looktxt -f -F=\"mynums\" -l=100 -v foo.txt' \n");
					printf("           will create/overwrite the file 'mynums.m' with only vectors/matrices containing at least 100 elements, and showing all the process and parameters.\n");
					printf("  'looktxt -a -g -r -o=\"oct\" foo.txt'\n");
					printf("           will search character and numeric fields, grouping them if possible, naming them with short names (such as 'n12' or 's10'), and creating 'foo_txt.oct' file\n");
					printf("Options : [ -p=\"{.|,}\" ]           to set point equivalents\n");
					printf("          [ -s=\"{\\t|\\v|\\r|,| |;}\" ]to set numeric separators\n");
					printf("          [ -c=\"{#|%%|c}\" ]         to set comment start character (till end of line)\n");
					printf("          [ -n=\"{\\n|\\f}\" ]         to set end of line equivalent character\n");
					printf("          [ -o=\"{oct|m|txt}\" ]     to set output type (extension) .oct, .m (default) or .txt file.\n");
					printf("          [ -F=\"filename\" ]        to force output filename to 'filename.extension'\n");
					printf("          [ -v ]                   to set verbose mode (show process and parameters)\n");
					printf("          [ -a ]                   to look also for ascii (character) fields.\n");
					printf("          [ -g ]                   to group numeric fields if possible\n");
					printf("          [ -r ]                   to use short field names (such as 'nXX' instead of 'filename_nXX')\n");
					printf("          [ -h ]                   this help\n");
					printf("          [ -f ]                   to force creation/overwriting of output file\n");
					printf("          [ -l=X ]                 to extract only num fields with at least X elements\n");
					printf("          [ -t ]                   to get field info 'table' matrix\n");
					printf("The 'table' field (with -t option) is a matrix of rows :\n");
					printf("   [ id type start_pos end_pos rows columns ]\n");
					printf("describing all other exported fields (identificator, type, start, end, dimensions)\n");
					printf("with type = 1 (num) | 2 (text) | 64 (comment)\n\n");
					printf("Type : looktxt -v to get default parameters.\n");
					printf(auteur " version " version " (" date ")\n");
					free(point);
					free(separator);
					free(eol);
					free(comment);

					free(filename);
					free(userfile);

					free(arg);
					exit(-1);
				}
			} /* if option */
			else
				strcpy(filename,arg);
		}
		if (verbose)
		{
			printf("* looktxt " version ": verbose mode\n");
			printf("point = \"%s\"\n",point);
			printf("separator = \"%s\"\n",undo_esc_char(separator));
			printf("eol = \"%s\"\n",undo_esc_char(eol));
			printf("comment = \"%s\"\n",comment);
			printf("filename = \"%s\"\n",filename);
			printf("Extract num fields with more than %i elements\n",line);
			if (passascii) printf("Pass ascii mode\n");
			if (groupnum) printf("Group numerics mode\n");
			if (noroot) printf("No root for fields names mode\n");
			if (outtype == OUT_OCT) printf("Output is .oct file\n");
			if (outtype == OUT_M) printf("Output is .m file\n");
			if (outtype == OUT_TXT) printf("Output is .txt file\n");
			if (forcefile) printf("Force creation of output file\n");
			if (tablextr) printf("Extract field info in table 'matrix'\n");
			if (username) printf("User output root filename \"%s\"\n",userfile);
		}	
		if (strlen(filename) > 0)
			scanfile(filename);
		else
			printf("no file to process\n");

		free(point);
		free(separator);
		free(eol);
		free(comment);

		free(filename);
		free(userfile);

		free(arg);
	} /* for arg */
	return(0);
}