/* 
    Code edited from the original threshold calculation to include
    binding of SOX2 and OCT4 by Madeline Galbraith
    last modified March 2020
*/
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <cstring>
# include <iostream>

int main (int argc, char **argv);
double Hillshift (double x, double x0, double nx, double lamda);
void estimate_threshold(int num);
double randd2(double minN,double maxN);
double randu2(double minV, double maxV); 

//global parameters //!!!

double threshA = 0.0;
double threshB = 0.0;
double threshC = 0.0;
double threshD = 0.0;
double threshE = 0.0;
double threshF = 0.0;
double threshG = 0.0;
double threshH = 0.0;
double threshI = 0.0;

# define NEQN 9

/******************************************************************************/
int main (int argc, char **argv)
{  

  estimate_threshold(10000);
    
  return 0;
}

/******************************************************************************/
double Hillshift (double x, double x0, double nx, double lamda)
{
    double out;

    out = lamda + (1.0 - lamda) * (1.0/(1.0 + pow((x/x0), nx)));
    
    return out;
}

/******************************************************************************/
double randd2(double minN,double maxN)
{
	double u=0.0;
	double z=0.0;
	int cnt = 0;
	int i = 0;

	do {
	     u = randu2(minN-1,maxN+1);
	     
             for (i=0; i< (maxN-minN+1); i++){
			if (u>=minN-0.5+i && u < minN+0.5+i){
				z = minN +i;
				cnt = 1;
				break;
			}
		}
	} while( u==maxN+0.5 || cnt ==0);

	return z;
}/******************************************************************************/
double randu2(double minV, double maxV)  // uniform random number between 0 and 1
{
    double u;
    do {
	u = ((double)rand()/RAND_MAX);	
    } while (u==0);

    return (minV+(maxV-minV)*u);	
}

/******************************************************************************/
double median(double *x,int n)
{
	double tmp;
	int 	i=0;
	int 	j=0;
	
	//Sorting
	for( i=0; i<n-1; i++){
		for (j=i+1; j<n; j++){
			if( x[j] <x[i]){
				tmp = x[i];
				x[i] = x[j];
				x[j] = tmp;
			}
		}
	}

	if(n%2==0){
		//printf("%f to %f \n",  x[n/2],x[n/2-1]);
		return ((x[n/2]+x[n/2-1])/2.0);
	}
	
	return x[n/2];
}


void estimate_threshold(int num)
{
	int i    = 0;
	int j    = 0;

	double k,g,n,T,lambda,MA,MB,MC,MD,ME;
	MA = 0.0;
	MB = 0.0;
	MC = 0.0;
	MD = 0.0;
	ME = 0.0;
	double minP = 1.;
	double maxP = 100.;
	double minK = 0.1;
	double maxK = 1.;
	double minN = 1.;
	double maxN = 6.;
	double minD = 0.1;
	double maxD = 1.;
	double minB = 0.01;
	double maxB = 0.1;
	double minU = 0.01;
	double maxU = 0.1;
	double minF = 1.;
	double maxF = 100.;

	double *A; // The standalone dist
	double *B;
	double *C;
	double *D;

	double SF = 1.0;
	double f1 = (2.0 - SF*1.96)/2.;
	double f2 = (2.0 + SF*1.96)/2.;

	A = (double *)calloc(num,sizeof(double));
	B = (double *)calloc(num,sizeof(double));
	C = (double *)calloc(num,sizeof(double));
	D = (double *)calloc(num,sizeof(double));
	
	int numI[] =  {1,0,1,0,3,0,2,0,0};
	int numA[] =  {1,2,1,1,2,1,0,0,0};
	int numB[] =  {0,0,0,0,0,0,0,2,0};
	int numIc[] = {1,0,0,0,0,0,0,0,0};
	int numAc[] = {0,0,0,0,1,0,0,0,0};
	int numId[] = {0,0,1,0,0,0,0,0,0};
	int numAd[] = {1,0,0,2,0,0,0,0,0};

	// Unb: 0, 0, 0, 0, 0, 0, 1, 0, 1
	//A Gata6 //B Gcnf //C Cdx2 //D Klf4 //E Nanog //F Pbx1 //G Oct4 //H Oct4-Sox2 //I Sox2

	int maxSimNum = 10.;

	for( int simNum=0; simNum < maxSimNum; simNum+=1){
	  	srand(simNum*2.);
		//Standalone dist for all but complex
		for( i=0; i< num; i++){
			g = randu2(minP,maxP);
			k = randu2(minK,maxK);
			A[i] = g/k;
		}
		MA = median(A,num);
		

		//Find estimated threshold for a standalone gene with binding eg Sox2
		for(i=0; i<num; i++){
			double g1= randu2(minP,maxP);
			double k1= randu2(minK,maxK);
			double g2= randu2(minP,maxP);
			double k2= randu2(minK,maxK);
			double g3= randu2(minP,maxP);
			double k3= randu2(minK,maxK);
			double rb = randu2(minB,maxB);
			double ru = randu2(minU,maxU);
			double kk = randu2(minD,maxD);
			C[i] = (g1+ru*(g2/k2))/(k1+(rb+kk)*(g3/k3));
			
		}
		MD = median(C,num);
		//printf("C %f\n", MD);

		// Find the estimated threshold of the complex 
		for (i=0; i<num ;i++){
			B[i] = randu2(minB,maxB)/randu2(minU,maxU);

			for(j=0; j< numB[7]; j++){
				double g1= randu2(minP,maxP);
				double k1= randu2(minK,maxK);
				double g2= randu2(minP,maxP);
				double k2= randu2(minK,maxK);
				double g3= randu2(minP,maxP);
				double k3= randu2(minK,maxK);
				double rb = randu2(minB,maxB);
				double ru = randu2(minU,maxU);
			double kk = randu2(minD,maxD);

				B[i] = B[i]*(g1+ru*(g2/k2))/(k1+(rb+kk)*(g3/k3));
			}
		}
		MC = median(B,num);
		//printf("C %f\n", MC);

		//Find estimated threshold for OCt4
		for(i=0; i<num; i++){
			double g1= randu2(minP,maxP);
			double k1= randu2(minK,maxK);
			double g2= randu2(minP,maxP);
			double k2= randu2(minK,maxK);
			double g3= randu2(minP,maxP);
			double k3= randu2(minK,maxK);
			double rb = randu2(minB,maxB);
			double ru = randu2(minU,maxU);
			double kk = randu2(minD,maxD);
			C[i] = g1/(k1+(rb+kk)*(g3/k3));

			for(j =0; j<2; j++){
				g = randu2(minP,maxP);
				k = randu2(minK,maxK);
				n = randd2(minN,maxN);
				T = randu2(MA*f1,MA*f2);
				lambda = 1.0/randu2(minF,maxF);
				C[i] = C[i]*Hillshift(g/k,T,n,lambda);
			}
			C[i] = C[i]+(ru*(g2/k2))/(k1+(rb+kk)*(g3/k3));
		}//
		ME = median(C,num);
		//printf("D %f\n", ME);

		for( int ng=0; ng <9; ng++){
			if(ng !=7 && ng!=6 && ng!=8){
				for( i=0; i< num; i++){
					g = randu2(minP,maxP);
					k = randu2(minK,maxK);
					B[i] = g/k;

					for(j=0; j< numA[ng]; j++){
						g = randu2(minP,maxP);
						k = randu2(minK,maxK);
						n = randd2(minN,maxN);
						T = randu2(MA*f1,MA*f2);
						lambda = randu2(minF,maxF);
						
						B[i] = B[i]*Hillshift(g/k,T,n,lambda)/lambda;
					}

					for(j=0; j< numI[ng]; j++){
						g = randu2(minP,maxP);
						k = randu2(minK,maxK);
						n = randd2(minN,maxN);
						T = randu2(MA*f1,MA*f2);
						lambda = 1.0/randu2(minF,maxF);
					
						B[i] = B[i]*Hillshift(g/k,T,n,lambda);
					}
					for(j=0; j< numB[ng]; j++){
						g = randu2(minP,maxP);
						k = randu2(minK,maxK);
						B[i] = B[i]*g/k;
					}
			
					for(j=0; j< numAc[ng]; j++){//activated by complex
						double b = randu2(minB,maxB)/randu2(minU,maxU);
						for(int k=0; k<2; k++){
							double g1= randu2(minP,maxP);
							double k1= randu2(minK,maxK);
							double g2= randu2(minP,maxP);
							double k2= randu2(minK,maxK);
							double g3= randu2(minP,maxP);
							double k3= randu2(minK,maxK);
							double rb = randu2(minB,maxB);
							double ru = randu2(minU,maxU);
							double kk = randu2(minD,maxD);

							b = b*(g1+ru*(g2/k2))/(k1+(rb+kk)*(g3/k3));
						}

						n = randd2(minN,maxN);
						T = randu2(MC*f1,MC*f2);
						lambda = randu2(minF,maxF);
					
						B[i] = B[i]*Hillshift(b,T,n,lambda)/lambda;
					}

					for(j=0; j< numIc[ng]; j++){//inhibited by complex
						double b = randu2(minB,maxB)/randu2(minU,maxU);
						for(int k=0; k<2; k++){
							double g1= randu2(minP,maxP);
							double k1= randu2(minK,maxK);
							double g2= randu2(minP,maxP);
							double k2= randu2(minK,maxK);
							double g3= randu2(minP,maxP);
							double k3= randu2(minK,maxK);
							double rb = randu2(minB,maxB);
							double ru = randu2(minU,maxU);
							double kk = randu2(minD,maxD);

							b = b*(g1+ru*(g2/k2))/(k1+(rb+kk)*(g3/k3));
						}
						n = randd2(minN,maxN);
						T = randu2(MC*f1,MC*f2);
						lambda = 1.0/randu2(minF,maxF);
						
						B[i] = B[i]*Hillshift(b,T,n,lambda);
					}	

					for(j=0; j< numId[ng]; j++){//inhibited by oct4 or sox2
						double b = (randu2(minP,maxP)+randu2(minU,maxU)*randu2(minP,maxP)/randu2(minK,maxK))/(randu2(minK,maxK)+(randu2(minB,maxB)+randu2(minD,maxD))*randu2(minP,maxP)/randu2(minK,maxK));
						n = randd2(minN,maxN);
						T = randu2(MD*f1,MD*f2);
						lambda = 1.0/randu2(minF,maxF);
						
						B[i] = B[i]*Hillshift(b,T,n,lambda);
					}	
					for(j=0; j< numAd[ng]; j++){//inhibited by oct4 or sox2
						double b = (randu2(minP,maxP)+randu2(minU,maxU)*randu2(minP,maxP)/randu2(minK,maxK))/(randu2(minK,maxK)+(randu2(minB,maxB)+randu2(minD,maxD))*randu2(minP,maxP)/randu2(minK,maxK));
						n = randd2(minN,maxN);
						T = randu2(MD*f1,MD*f2);
						lambda = randu2(minF,maxF);
						
						B[i] = B[i]*Hillshift(b,T,n,lambda)/lambda;
					}	
				}
			}//end if ng!=7
			MB = median(B,num);
			//printf("B %d-%f\n", ng,MB);

			switch (ng){
				case 0: //A Gata6
					threshA = threshA + MB;
				break;
				case 1: //B Gcnf
					threshB = threshB + MB;
				break;
				case 2: //C Cdx2
					threshC = threshC + MB;
				break;
				case 3: //D Klf4
					threshD = threshD + MB;
				break;
				case 4: //E Nanog
					threshE = threshE + MB;
				break;
				case 5: //F Pbx1
					threshF = threshF + MB;
				break;
				case 6: //G Oct4
					threshG = threshG + ME;
				break;
				case 7: //H Oct4-Sox2
					threshH = threshH + MC;
					//printf("C %f\n", MC);
				break;
				case 8: //I Sox2
					threshI = threshI + MD;
				break;
			}

		}//end of genes
	}// end of averaging

	threshA = threshA/maxSimNum;
	threshB = threshB/maxSimNum;
	threshC = threshC/maxSimNum;
	threshD = threshD/maxSimNum;
	threshE = threshE/maxSimNum;
	threshF = threshF/maxSimNum;
	threshG = threshG/maxSimNum;
	threshH = threshH/maxSimNum;
	threshI = threshI/maxSimNum;

        printf("Gata6\t\t- %f\n", threshA);
        printf("Gcnf\t\t- %f\n", threshB);
        printf("Cdx2\t\t- %f\n", threshC);
        printf("Klf4\t\t- %f\n", threshD);
        printf("Nanog\t\t- %f\n", threshE);
        printf("Pbx1\t\t- %f\n", threshF);
     	printf("Oct4\t\t- %f\n", threshG);
     	printf("Oct4Sox2\t- %f\n", threshH);
     	printf("Sox2\t\t- %f\n", threshI);

	free(A);
	free(B);
}
