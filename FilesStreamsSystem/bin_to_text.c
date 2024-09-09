/*  simple code for converting binary to text output  */
/*  bin_to_text [#fields] input.bin output.txt        */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main( int argc , char *argv[] ){
	FILE *fp_in,*fp_out;
	int nf,i;
	float *data;
	
	if( argc == 4 ) 
		fp_in = fopen(argv[2],"rb");
	else return(1);
	fp_out = fopen(argv[3],"w");
	
	nf = atoi(argv[1]);
	data = (float *)malloc(nf*sizeof(float));
	
	while ( fread(data,nf*sizeof(float),1,fp_in ) == 1 ) {
		for(i=0;i<nf;i++) fprintf(fp_out,"%f ",data[i]);
		fprintf(fp_out,"\n");
	}
	
	fclose(fp_in);
	fclose(fp_out);
	return(0);
}
