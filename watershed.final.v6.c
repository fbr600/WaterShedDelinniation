
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "watershed.final.v6.h"
    int n = 1386; //columns	Catalunya
    int m = 1335; //rows	Catalunya

//    int n = 1008; //columns	Catalunya
//    int m = 970; //rows	Catalunya

//	int n = 3080; //columns	Catalunya
//	int m = 2966; //rows	Catalunya


//    int n = 7040; //columns	Andorra
//    int m = 5540; //rows	Andorra

int main()
{
	float *altitude;
	short *dir,*isboundary,*nprev;
	int *next, **ant,*ppin,*ppout,*index,*id,*flowacc,i,j,k;

	FILE *fin,*fout;

	//catalunya
	fin  = fopen("MDE200m_ICC_Aster.img","rb");
	fout = fopen("Catalunya_watershed.img","wb");

	//catalunya menys pixels
//	fin  = fopen("MDE200m_ICC_DENS_Aster.img","rb");
//	fout = fopen("Catalunya_watershed_DENS.img","wb");

	//catalunya més pixels resolució
//	fin  = fopen("MDE200m_gran_ICC_Aster.IMG","rb");
//	fout = fopen("Catalunya_watershed_DENS.img","wb");

    //andorra
//	fin  = fopen("MDE_Andorra.img","rb");
//	fout = fopen("Andorra_watershed.img","wb");

	altitude 	= ( float* )malloc((m*n)*sizeof( float ));
	dir 		= ( short* )malloc((m*n)*sizeof( short ));
	isboundary 	= ( short* )malloc((m*n)*sizeof( short ));
	next 		= ( int*   )malloc((m*n)*sizeof( int   ));
	ppin 		= ( int*   )malloc((m*n)*sizeof( int   ));
	ppout 		= ( int*   )malloc((m*n)*sizeof( int   ));
	nprev		= ( short* )malloc((m*n)*sizeof( short ));
	index 		= ( int *)malloc((m*n)*sizeof(int));
	id	 		= ( int *)malloc((m*n)*sizeof(int));
	flowacc		= ( int *)malloc((m*n)*sizeof(int));

	//Reading input data
	InitDataBinnary(altitude,fin,m,n);

	i = 45;

	//Cleaning input data
	for(i=0;i<m*n;i++) if(altitude[i]<0) altitude[i] = 0;

	//Flow Direction
	Flowdir(altitude,dir,m,n);/*OK dir en el mar o en nodata és 0*/

	//Set Next
	SetNext(dir,next,altitude,m,n);

	//Set Previous
	ant = SetPrev(ant,nprev,next,index,dir,m,n);

	//Set Basin Id
	j = SetID(id,next,ant,index,nprev,altitude,m,n);

	//Boundaries
	Boundaries(isboundary,id,m,n);

	//Pour points
	Ppoint(ppin,ppout,altitude,isboundary,next,dir,id,m,n,j);

	//Setting new directions (merging small basin)
	Newdirections(next,ant,dir,altitude,id,ppin,ppout,index,nprev,m,n,j);

	//Set previous
	free(ant);
	ant = NULL;
	ant = SetPrev(ant,nprev,next,index,dir,m,n);

	//Final Basin ID
	j = SetID(id,next,ant,index,nprev,altitude,m,n);

	//Flow Accumulation (Optional)
	//Flowaccum(flowacc,nprev,ant,index,m,n);

	//Write the final layer
	WriteDataBinnary(id,fout,m,n);
}
