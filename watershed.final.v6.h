
#include <stdio.h>
#include <stdlib.h>

#ifndef LIBRARYHEADERFILE_INCLUDED
#define LIBRARYHEADERFILE_INCLUDED

/******************************************************************************************/
/* HEADERS ********************************************************************************/
/******************************************************************************************/

/*modelling functions**********************************************************************/
void 	Flowdir		 ( float *altitude ,short *dir,int m, int n);
void	SetNext		 ( short *dir ,int *next,float *alt,int m,int n);
int**	SetPrev		 ( int   **ant ,short *nprev,int *next,int *index,short *dir,int m,int n);
void	Setnprev	 ( short *nprev ,int *next,int m,int n);
int		SetID		 ( int   *id ,int *next,int **ant,int *index,short *nprev,float *altitude,int m,int n);
void 	Completeup	 ( int   run ,int *id,int j,int **ant,int *index,short *nprev,int m,int n);
void 	Boundaries	 ( short *isboundary,int *id,int m,int n);
void	Ppoint		 ( int   *ppin ,int *ppout,float *altitude,short *isboundary,int *next,short *dir,int *id,int m,int n,int numreg);
void 	Newdirections( int   *next ,int **ant,short *dir,float *altitude,int *id,int *ppin,int *ppout,int *index,short *nprev,int m,int n,int numreg);
void	Flowaccum	 ( int *flowacc,short *nprev,int **ant,int *index,int m,int n);
int 	Accum		 ( int position,int **ant,short *nprev,int *index,int m,int n);

/*technical functions**********************************************************************/
void 	PointsNear		( float *altitude, int i, float *b,int m,int n);
void	PosofFloat		( float*,int*,int,int,int);
void 	Positions		( int *positions,int pos,int m,int n);
int		Dir				( int);
int 	WhichMin		( float *b);
int		Isintheborder	( int i, int m, int n);
void  	Setindex		( int *index,short *nprev,int m,int n);


/*reading and writing functions************************************************************/
void 	InitDataBinnary				( float *a, FILE *f,int m, int n);
void 	WriteMatrixDecimalCSVshort	( short *a,int m,int n,FILE *f);
void 	WriteMatrixDecimalCSVint	( int *a,int m,int n,FILE *f);
void	WriteDataBinnaryShort		( short *a,FILE *f,int m, int n);
void	WriteDataBinnary			( int *a,FILE *f,int m, int n);


/*global variables*************************************************************************/
float NODATA = -9999;
int pos[9];
int counter=0;

/******************************************************************************************/
/* FUNCTIONS ******************************************************************************/
/******************************************************************************************/

/*modelling functions**********************************************************************/

//matrix dir has at the end the SFD8 codification for directions
void 	Flowdir(float *altitude,short *dir,int m, int n)
{
	printf("Flow Direction\n");
	int i,k;
	float *b;
	b = (float*)malloc(9*sizeof(float));
	for(i=0;i<m*n;i++)
	{
		PointsNear(altitude,i,b,m,n);
		k = WhichMin(b);
		dir[i] = Dir(k);
	}
	free(b);
	return;
}

//next is the matrix in which positions has the next cell position
void	SetNext(short *dir,int *next,float *alt,int m,int n)
{
	printf("Next positions\n");
	int i,j,posx;
	float minalt = 10000;

	for(i=0;i<m*n;i++)
	{
		j=0;
		Positions(pos,i,m,n);
		if(dir[i]==0)
		{
			//if pos == -1 we can  not evaluate dir[pos] -> out of map or NOADATA
			if(pos[0]==-1) j++; else if( dir[pos[0]] == 2   ) j++;
			if(pos[1]==-1) j++; else if( dir[pos[1]] == 4   ) j++;
			if(pos[2]==-1) j++; else if( dir[pos[2]] == 8   ) j++;
			if(pos[3]==-1) j++; else if( dir[pos[3]] == 1   ) j++;
			if(pos[5]==-1) j++; else if( dir[pos[5]] == 16  ) j++;
			if(pos[6]==-1) j++; else if( dir[pos[6]] == 128 ) j++;
			if(pos[7]==-1) j++; else if( dir[pos[7]] == 64  ) j++;
			if(pos[8]==-1) j++; else if( dir[pos[8]] == 32  ) j++;

			if(j==8 || j==0) next[i]=-1;
			else {dir[i] = 3; next[i]=-1;}
		}
		else
		{
			if( dir[i]==1  ) next[i] = pos[5];
			if( dir[i]==2  ) next[i] = pos[8];
			if( dir[i]==4  ) next[i] = pos[7];
			if( dir[i]==8  ) next[i] = pos[6];
			if( dir[i]==16 ) next[i] = pos[3];
			if( dir[i]==32 ) next[i] = pos[0];
			if( dir[i]==64 ) next[i] = pos[1];
			if( dir[i]==128) next[i] = pos[2];

		}
		if(Isintheborder(i,m,n)==1) next[i]=-1;
	}
	for(i=0;i<m*n;i++)
	{
		if(dir[i]==3)
		{
			if(Isintheborder(i,m,n)==1) {dir[i] = 0; next[i]=-1; continue;}
			Positions(pos,i,m,n);
			minalt = 10000;
			j=0;
			posx=-1;
			if(pos[0]!=-1) if( dir[pos[0]] != 2   ) if(alt[next[pos[0]]] <minalt ) {minalt =alt[next[pos[0]]]; posx = next[pos[0]]; }
			if(pos[1]!=-1) if( dir[pos[1]] != 4   ) if(alt[next[pos[1]]] <minalt ) {minalt =alt[next[pos[1]]]; posx = next[pos[1]]; }
			if(pos[2]!=-1) if( dir[pos[2]] != 8   ) if(alt[next[pos[2]]] <minalt ) {minalt =alt[next[pos[2]]]; posx = next[pos[2]]; }
			if(pos[3]!=-1) if( dir[pos[3]] != 1   ) if(alt[next[pos[3]]] <minalt ) {minalt =alt[next[pos[3]]]; posx = next[pos[3]]; }
			if(pos[5]!=-1) if( dir[pos[5]] != 16  ) if(alt[next[pos[5]]] <minalt ) {minalt =alt[next[pos[5]]]; posx = next[pos[5]]; }
			if(pos[6]!=-1) if( dir[pos[6]] != 128 ) if(alt[next[pos[6]]] <minalt ) {minalt =alt[next[pos[6]]]; posx = next[pos[6]]; }
			if(pos[7]!=-1) if( dir[pos[7]] != 64  ) if(alt[next[pos[7]]] <minalt ) {minalt =alt[next[pos[7]]]; posx = next[pos[7]]; }
			if(pos[8]!=-1) if( dir[pos[8]] != 32  ) if(alt[next[pos[8]]] <minalt ) {minalt =alt[next[pos[8]]]; posx = next[pos[8]]; }
			next[i] = posx;
		}
	}
	return;
}

//each ant[i] is the i+1-st element that has previous,and will have exactly the number of
//anteriors.
//The translation between positions i and the ones that has previous is through the vector
//index index[i] is -1 if i has not previous, and i is the index[i]-th element with previous
//if i has previous
int**	SetPrev(int **ant,short *nprev,int *next,int *index,short *dir,int m,int n)
{
	int i,j,k,l,maxnprev=0,aux;
	printf("Previous positions\n");
	Setnprev(nprev,next,m,n);
	for(i=0;i<m*n;i++)
	{
		if(nprev[i]>maxnprev) {maxnprev = nprev[i];}
	}
	j=0;
	for(i=0;i<m*n;i++)
	{

		if(nprev[i]>0)
		{
			//if(counter==1) printf("i %d, nprev %d\n",i,nprev[i]);
			if(j==0)
			{
				ant = (int**)malloc(1*sizeof(int*));
				j++;
				ant[0] = (int*)malloc((1+nprev[i])*sizeof(int));
				ant[0][0] = i;
				for(k=0;k<nprev[i];k++) ant[0][k+1] = -1;
			}
			else
			{
				j++;
				ant = (int**)realloc(ant,j*sizeof(int*));
				ant[j-1] = (int*)malloc((1+nprev[i])*sizeof(int));
				ant[j-1][0] = i;
				for(k=0;k<nprev[i];k++) ant[j-1][k+1] = -1;
			}
			index[i] = j-1;
		}
		else
			index[i] = -1;
	}
	/*j is the number of cells with at least one previous	*/
	for(i=0;i<m*n;i++)
	{
		if(next[i]!=-1 && index[next[i]]!=-1)
		{
			k=0;
			while(ant[index[next[i]]][k]!=-1) {k++;};
			ant[index[next[i]]][k] = i;
 		}
	}
	return(ant);
}
//nprev[i] is the number of anterior elements that i has
void 	Setnprev(short *nprev,int *next,int m,int n)
{
	int i;
	for(i=0;i<m*n;i++)
	{
		if(Isintheborder(next[i],m,n)==1) continue;
		nprev[next[i]]++;
	}
	return;
}
//id is the layer with the labels of the basin
int		SetID(int *id,int *next,int **ant,int *index,short *nprev,float *altitude,int m,int n)
{
	printf("Set ID\n");
	int i,j,run,aux,*numid;
	j=0;
	for(i=0;i<m*n;i++) id[i] = 0;
	for(i=0;i<m*n;i++)
	{
		if(id[i]!=0) continue; /* if id has been set before */
		if(next[i]==-1 && nprev[i]==0 && (Isintheborder(i,m,n)==1 || altitude[i]==0)) {id[i]=-1;continue;} /* sea or map border */
 		if(nprev[i]==0) /*starting point*/
		{
			j++;
			id[i]=j;
			run = i;
			while(Isintheborder(run,m,n)==0 && next[run]!=-1 )
			{
				run = next[run];
 				Completeup(run,id,j,ant,index,nprev,m,n);
			};
		}
	}
	for(i=0;i<m*n;i++) if(Isintheborder(i,m,n)==1) id[i] = -1;
	numid = (int*)malloc((j+1)*sizeof(int));
	for(i=0;i<=j;i++) numid[i]=i;
	for(i=0;i<m*n;i++)
	{
		if(id[i]!=-1) numid[id[i]]=0;
	}
	aux = 0;
	for(i=0;i<=j;i++) if(numid[i]==0) aux++;
	free(numid);
	printf("%d regions\n",aux);
	return(j);
}
//this function assign an id up slope when two steams join
void 	Completeup(int run,int *id,int j,int **ant,int *index,short *nprev,int m,int n)
{
	int i;
	if(run==-1) return;
	id[run]=j;
	if(nprev[run]==0) return;
	for(i=0;i<nprev[run];i++)
	{
		if(id[ant[index[run]][i+1]] == j ) continue;
		else {Completeup(ant[index[run]][i+1],id,j,ant,index,nprev,m,n);}
	}
	return;
}
//Creates a new layer showing if each point is on the boundary of the basin or not
void 	Boundaries(short *isboundary,int *id,int m,int n)
{
	int i,k;
	printf("Boudaries \n");
	for(i=0;i<m*n;i++)
	{
		Positions(pos,i,m,n);
		for(k=0;k<9;k++)
			if(pos[k]!=-1)
				if(id[pos[k]]!=id[i])
					isboundary[i] = 1;
	}
	return;
}

//Setting pour points inside and out side for each basin
void	Ppoint(int *ppin,int *ppout,float *altitude,short *isboundary,int *next,short *dir, int *id,int m,int n,int numreg)
{
	printf("Pour Points\n");
	int i,j,k,*hashpp,*hashppout,posx;
	float *altitudein,*altitudeout;
	//struct cell *b,*c;
	hashpp = (int*)malloc((numreg+1)*sizeof(int));
	hashppout = (int*)malloc((numreg+1)*sizeof(int));
	altitudein = (float*)malloc((numreg+1)*sizeof(float));
	altitudeout = (float*)malloc((numreg+1)*sizeof(float));
	for(i=0;i<=numreg;i++) {altitudein[i] =10000; altitudeout[i]=10000;}
	for(i=0;i<m*n;i++)
	{
		if(Isintheborder(i,m,n)==1) continue;
		Positions(pos,i,m,n); /*split in two cases: the normal one and the links*/
		for(k=0;k<9;k++)
		{
			if(pos[k]!=-1 /**/&& Isintheborder(pos[k],m,n)==0)
			{
				if(isboundary[i]==1 && id[pos[k]]==id[i] && altitude[pos[k]]<altitudein[id[i]])
				{
					altitudein[id[i]] = altitude[pos[k]];
					hashpp[id[i]] = pos[k];
				}
				if(isboundary[i]==1 && id[pos[k]]!=id[i] && altitude[pos[k]]<altitudeout[id[i]] && next[pos[k]]!=-1) /**/
				{
					altitudeout[id[i]] = altitude[pos[k]];
					hashppout[id[i]] = pos[k];
				}
			}
		}
	}
	for(i=0;i<m*n;i++)
	{
		ppin[i] = hashpp[id[i]];
		ppout[i] = hashppout[id[i]];
	}
	/*Correction of the ppout in the cases of lakes and wide rivers the id of this regions are -1 but we
	don't want to connect this zones with the borders and sea points*/
	for(i=0;i<m*n;i++)
	{
		if(next[i]==-1 && altitude[i]>0.001 && Isintheborder(i,m,n)==0)
		{
			Positions(pos,i,m,n);
			/*****/
			for(k=0;k<9;k++)
			{
				if(next[pos[k]]==-1 && ppout[pos[k]]!=ppout[0] )
				{
					ppout[i] = ppout[pos[k]];
				}
			}
			if(ppout[i]!=ppout[0]) continue;
			for(k=0;k<9;k++)
			{
				if(ppout[pos[k]]!=ppout[0] )
					{ppout[i] = ppout[pos[k]]; printf("ASDF\n"); }
			}
			/*******/
			/*******/
		}
	}

	free(hashpp); free(hashppout);
	free(altitudein); free(altitudeout);
	return;
}
//this function assign new directions to the pits of the basins
void 	Newdirections(int *next,int **ant,short *dir,float *altitude,int *id,int *ppin,int *ppout,int *index,short *nprev,int m,int n,int numreg)
{
	counter ++;

	printf("New Directions\n");
	int i,j,k,l,a,b,c,d,aux,posx,posy;
	float minalt;

	for(i=0;i<m*n;i++)
	{
		aux = 0;
		if(next[i]==-1 && altitude[i]>0.0001 && Isintheborder(i,m,n)==0)
		{
			Positions(pos,i,m,n); /*isolated point with this characteristics*/
			for(k=0;k<9;k++)
				if(pos[k]!=-1 && next[pos[k]]==-1)
					aux++;

			if(aux==0) continue;

			if(aux==1)
			{
				a = ppin[i];
				b = -1;
				Positions(pos,a,m,n);
				minalt =10000;
				for(k=0;k<9;k++)
					if(pos[k]!=-1 && id[pos[k]]!=-1 && id[pos[k]]!=id[i] && altitude[pos[k]]<minalt)
					{	b = pos[k]; minalt = altitude[pos[k]];}

				c = ppout[i];
				d = -1;
				Positions(pos,a,m,n);
				minalt =10000;
				for(k=0;k<9;k++)
					if(pos[k]!=-1 && id[pos[k]]!=-1 && id[pos[k]]==id[a] && altitude[pos[k]]<minalt)
					{	d = pos[k]; minalt = altitude[pos[k]];}

				if(b!=-1 && d!=-1)
				{

					if(altitude[b]>altitude[d])
					{
						if(altitude[i]>=altitude[c])
						{
							next[i] = c;
							nprev[c]++;
							continue;
						}
					}
					else
					{
						if(altitude[i]>=altitude[b])
						{
							next[i]= b;
							nprev[b]++;
							continue;
						}
					}
				}
			}
			/*searching a lower point when the flat zones*/
			j=i;
			k=0;
			while(altitude[j]>=altitude[i])
			{
				k++;
				for(l=-k;l<=k;l++)
				{	posx = i-k*n+l; if(posx%n==0 ||posx%n==m-1 ||posx<0|| posx>m*n) continue;
					if(altitude[posx]==NODATA) continue;
					if(altitude[posx]<altitude[i]) j = posx;}
				for(l=-k+1;l<=k-1;l++)
				{	posx = i-l*n-k; if(posx%n==0 ||posx%n==m-1|| posx<0|| posx>m*n) continue;
					if(altitude[posx]==NODATA) continue;
					if(altitude[posx]<altitude[i]) j = posx;}
				for(l=-k+1;l<=k-1;l++)
				{	posx = i-l*n+k; if(posx%n==0 ||posx%n==m-1|| posx<0|| posx>m*n) continue;
					if(altitude[posx]==NODATA) continue;
					if(altitude[posx]<altitude[i]) j = posx;}
				for(l=-k;l<=k;l++)
				{	posx = i+k*n+l; if(posx%n==0 ||posx%n==m-1|| posx<0|| posx>m*n) continue;
					if(altitude[posx]==NODATA) continue;
					if(altitude[posx]<altitude[i]) j = posx;}
			};
			next[i] = j ;
			nprev[j]++;
		}
	}

	return;
}
//Flow Accumulation
void	Flowaccum	 ( int *flowacc,short *nprev,int **ant,int *index,int m,int n)
{
	int i,j,sum;
	for(i=0;i<m*n;i++) flowacc[i]=0;
	for(i=0;i<m*n;i++)
	{
		sum = 0;
		if(nprev[i]==0) {flowacc[i]=0; continue;}
		else
		{
			if(flowacc[i]!=0) continue;
			for(j=1;j<=nprev[i];j++)
			{
				sum = sum+Accum(ant[index[i]][j],ant,nprev,index,m,n);
			}
			flowacc[i] = sum+nprev[i];
		}
	}
	return;
}
//Recursive function that computes the flow Accumulation up slope
int 	Accum(int position,int **ant,short *nprev,int *index,int m,int n)
{
	int j,sum=0;
	if(position==-1) return(0);
	if(nprev[position]==0) return(0);
	else
	{
		for(j=1;j<=nprev[position];j++)
		{
			sum = sum+Accum(ant[index[position]][j],ant,nprev,index,m,n);
		}
		return(sum+nprev[position]);
	}
	return(0);
}

/*technical functions**********************************************************************/

//b is a vector with the altitude of the neighbours of i
//in the cases of out of limits or NODATA, it returns 10000, since a minimum is sought
//we want this points to not take importance
void 	PointsNear(float *altitude, int i, float *b,int m,int n)
{
	int *positions;
	short j;

	PosofFloat(altitude,pos,i,m,n);
	for(j=0;j<9;j++)
	{
		if(pos[j]==-1) b[j] = 10000;
		else b[j] = altitude[pos[j]];
	}
	return;
}

//positions returns -1 when the neighbours are outside the limits or when the altitude
//layer is NODATA in those points
void	PosofFloat(float *alt,int *positions,int pos,int m,int n)
{
	int i;
	Positions(positions,pos,m,n);
	for(i=0;i<9;i++)
	{
		if(positions[i]!=-1) if(alt[positions[i]]==NODATA) positions[i]=-1;
	}
	return;
}
//return the positions of the neighbours of pos
//As before, if the neighbours are outside the limits are set to -1
void 	Positions(int *positions,int pos,int m,int n)
{
	int x,y;
	x = pos%n;
	y = (pos -x)/n;
	positions[0] = pos-n-1;
	positions[1] = pos-n;
	positions[2] = pos-n+1;
	positions[3] = pos-1;
	positions[4] = pos;
	positions[5] = pos+1;
	positions[6] = pos+n-1;
	positions[7] = pos+n;
	positions[8] = pos+n+1;

	if(x==0 && y==0)
	{ 	positions[0] = -1;positions[1] = -1;positions[2] = -1;positions[3] = -1;positions[6] = -1; }
	if(x>0 && x<n-1 && y==0)
	{	positions[0] = -1;positions[1] = -1;positions[2] = -1; }
	if(x==n-1 && y==0)
	{	positions[0] = -1;positions[1] = -1;positions[2] = -1;positions[5] = -1;positions[8] = -1;	}
	if(x==0 && y>0 && y<m-1)
	{	positions[0] = -1;positions[3] = -1;positions[6] = -1;	}
	if(x==n-1 && y>0 && y<m-1)
	{	positions[2] = -1;positions[5] = -1;positions[8] = -1;	}
	if(x==0 && y==m-1)
	{ 	positions[0] = -1;positions[3] = -1;positions[6] = -1;positions[7] = -1;positions[8] = -1;	}
	if(x>0 && x<n-1 && y==m-1)
	{ 	positions[6] = -1;positions[7] = -1;positions[8] = -1;	}
	if(x==n-1 && y == m-1)
	{ 	positions[2] = -1;positions[5] = -1;positions[6] = -1;positions[7] = -1;positions[8] = -1;	}
	return;
}
//return which is the minimum of a vector of 9 real components
//if the 4th position does not contain the minimum, then the returned minimum is the first.
//otherwise, always 4 is returned
int 	WhichMin(float *b)
{
	int i,mini,j=0;
	float mn;
	mn = 10000;
	j=0;
	for(i=0;i<9;i++)
	{
		if(b[i]<mn && i!=4) {mn = b[i]; mini=i;}
	}
	if(mn>=b[4]) return(4);
	return(mini);
}
//given the position of the flow direction in terms of a 9-neighbours vector,
//is returned the code of SFD8
int 	Dir(int k)
{
	if(k==0) 	return(32);
	if(k==1) 	return(64);
	if(k==2)	return(128);
	if(k==3)	return(16);
	if(k==4)	return(0);
	if(k==5)	return(1);
	if(k==6)	return(8);
	if(k==7)	return(4);
	if(k==8)	return(2);
}
//given a position i of the matrix, it returns 0 if is not on the border of the matrix
//and return 1 if yes
int 	Isintheborder(int i,int m,int n)
{
	int x,y;
	x = i%n;
	y = (i -x)/n;
	if(x==n-1 || x==0 || y==m-1 || y==0) return(1);
	return(0);
}
//given the nprev information, it assign the index vector, containing:
//index[i] is -1 if i has not previous, and i is the index[i]-th element with previous
//if i has previous
void  	Setindex(int *index,short *nprev,int m,int n)
{
	int i,j;

	j = 0;
	for(i=0;i<m*n;i++)
	{
		if(nprev[i]>0) {j++; index[i] = j-1;}
		else index[i] = -1;
	}
}


/*reading and writing functions************************************************************/
void 	InitDataBinnary(float *a, FILE *f,int m, int n)
{
	int i;
	i=fread(a,4,m*n,f);
	return;
}
void 	WriteMatrixDecimalCSVshort(short *a,int m,int n,FILE *f)
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			fprintf(f,"%4d,",a[n*i+j]);
		fprintf(f,"\n");
	}
	return;
}
void 	WriteMatrixDecimalCSVint(int *a,int m,int n,FILE *f)
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			fprintf(f,"%4d,",a[n*i+j]);
		fprintf(f,"\n");
	}
	return;
}
void	WriteDataBinnaryShort(short *a,FILE *f,int m, int n)
{
	fwrite(a,sizeof(short),m*n,f);
	return;
}
void	WriteDataBinnary(int *a,FILE *f,int m, int n)
{
	fwrite(a,sizeof(int),m*n,f);
	return;
}
void WriteNeighbours(float *b,int i,int n)
{

	printf("%2.7f,%2.7f,%2.7f\n",b[i-n-1],b[i-n],b[i-n+1]);
	printf("%2.7f,%2.7f,%2.7f\n",b[i-1],b[i],b[i+1]);
	printf("%2.7f,%2.7f,%2.7f\n\n",b[i+n-1],b[i+n],b[i+n+1]);
	return;
}
#endif




