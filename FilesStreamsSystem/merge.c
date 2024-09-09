/* merge two files: one line from one, then one line from another.
(These go in one line of output file, separated by tabs.)
If they don't have same no. of lines, complain and
quit when you hit first EOF */

/*    usage:   merge file1 file2 > output    */


#include <stdio.h>

main(argc,argv) int argc; char **argv; {
int i,j,c;
FILE *fopen(),*fp1,*fp2;
int flag1,flag2;

    if(argc != 3){fprintf(stderr,"Wrong no. of arguments\n");exit(-1);}
    if((fp1=fopen(argv[1],"r"))==NULL){
	fprintf(stderr,"Can't open file %s\n",argv[1]);
    }
    if((fp2=fopen(argv[2],"r"))==NULL){
	fprintf(stderr,"Can't open file %s\n",argv[2]);
    }
    flag1=flag2=0;

    while(1){

	/* take line from first file */
	for(;;){
	    c=getc(fp1);
	    if(c==EOF){
		flag1++;
		break;
	    }
	    if(c=='\n'){
		putchar('\t');
		break;
	    }
	    putchar(c);
	}

	/* take line from second file */
	for(;;){
	    c=getc(fp2);
	    if(c==EOF){
		flag2++;
		break;
	    }
	    if(c=='\n'){
		break;
	    }
	    putchar(c);
	}

	/* quit silently if both files ended, complain if only one ended */
	if( flag1 && flag2 )break;
	/* put newline if you got anything on this line */
	putchar('\n');
	if( flag1 ){
		fprintf(stderr,"Warning: file 1 ended before file 2\n");
		exit(1);
	}
	if( flag2 ){
		fprintf(stderr,"Warning: file 2 ended before file 1\n");
		exit(1);
	}

    }
    fclose(fp1);
    fclose(fp2);
}
