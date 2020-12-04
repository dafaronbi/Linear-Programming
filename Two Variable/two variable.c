#include <stdio.h>

void main(){
	//create file pointer
	FILE *fp;
	
	//create buffer for reading from MPS file
	char buff[1000];
	
	//read test MPS file
	fp = fopen("../MPS/test.mps", "r");
	
	fgets(buff,1000,fp);
	
	printf(buff);
	

}