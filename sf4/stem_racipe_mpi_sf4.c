# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <cstring>
# include <mpi.h>
# include <iostream>

int main (int argc, char **argv);
int stem (double y_store[], int j, int num);
void model ( double t, double y[], double yp[] );
double Hillshift (double x, double x0, double nx, double lamda);
double randu(void);
double randud(void);
void set_parameters (char **argv);
double abs_test(double x);
double sumdelta(double y[], double ytmp[]);
double log2(double x);
int count_state (double y_store[], int num_ode, double thrd, double soln[]);

//global parameters //!!!
double ga = 0.0;
double gb = 0.0;
double gc = 0.0;
double gd = 0.0;
double ge = 0.0;
double gf = 0.0;
double gg = 0.0;
double gh = 0.0;
double gI = 0.0;

double rb = 0.0;
double ru = 0.0;

double ka = 0.0;
double kb = 0.0;
double kc = 0.0;
double kd = 0.0;
double ke = 0.0;
double kf = 0.0;
double kg = 0.0;
double kh = 0.0;
double kI = 0.0;

double lamdaaa = 0.0;
double lamdaha = 0.0;
double lamdaea = 0.0;
double lamdaga = 0.0;
double lamdaab = 0.0;
double lamdacb = 0.0;
double lamdacc = 0.0;
double lamdaec = 0.0;
double lamdagc = 0.0;
double lamdaId = 0.0;
double lamdaed = 0.0;
double lamdagd = 0.0;
double lamdaee = 0.0;
double lamdafe = 0.0;
double lamdahe = 0.0;
double lamdade = 0.0;
double lamdaae = 0.0;
double lamdace = 0.0;
double lamdaef = 0.0;
double lamdabg = 0.0;
double lamdacg = 0.0;
double lamdahg = 0.0;
double lamdahI = 0.0;
double lamdagh = 0.0;
double lamdaIh = 0.0;


double a0a = 0.0;
double h0a = 0.0;
double e0a = 0.0;
double g0a = 0.0;
double a0b = 0.0;
double c0b = 0.0;
double c0c = 0.0;
double e0c = 0.0;
double g0c = 0.0;
double I0d = 0.0;
double e0d = 0.0;
double g0d = 0.0;
double e0e = 0.0;
double f0e = 0.0;
double h0e = 0.0;
double d0e = 0.0;
double a0e = 0.0;
double c0e = 0.0;
double e0f = 0.0;
double b0g = 0.0;
double c0g = 0.0;
double h0I = 0.0;
double h0g = 0.0;
double g0h = 0.0;
double I0h = 0.0;

double naa = 0.0;
double nha = 0.0;
double nea = 0.0;
double nga = 0.0;
double nab = 0.0;
double ncb = 0.0;
double ncc = 0.0;
double nec = 0.0;
double ngc = 0.0;
double nId = 0.0;
double ned = 0.0;
double ngd = 0.0;
double nee = 0.0;
double nfe = 0.0;
double nhe = 0.0;
double nde = 0.0;
double nae = 0.0;
double nce = 0.0;
double nef = 0.0;
double nbg = 0.0;
double ncg = 0.0;
double nhg = 0.0;
double nhI = 0.0;
double ngh = 0.0;
double nIh = 0.0;

# define NEQN 9

/******************************************************************************/
int main (int argc, char **argv)
{  

  int i = 0;
  int j = 0;
  int h = 0;
  int h2= 0;
  int num_paras = 10000;
  int num_ode   = 1000;
  int cnt = 0;
  int cnt_store[10] = {0};
  double thrd = 1; 
  double *y_store, *soln;	
  y_store = (double *)calloc(num_ode*NEQN, sizeof(double));
  soln    = (double *)calloc(10*NEQN, sizeof(double));

  int cnt_loop = 0.0;

  // Initialize everything necessary for mpi
  int my_rank, comm_sz;
  MPI_Init(NULL,NULL);
  //get num processes
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  //Get rank of process
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  
  //output
  MPI_File f_p;
  MPI_File f_s_1;
  MPI_File f_s_2;
  MPI_File f_s_3;
  MPI_File f_s_4;
  MPI_File f_s_5;
  MPI_File f_s_6;
  MPI_File f_s_7;
  MPI_File f_s_8;
  MPI_File f_s_9;
  MPI_File f_s_10;
  
  char fpname [100]="";
  char fsname1[100]="";
  char fsname2[100]="";
  char fsname3[100]="";
  char fsname4[100]="";
  char fsname5[100]="";
  char fsname6[100]="";
  char fsname7[100]="";
  char fsname8[100]="";
  char fsname9[100]="";
  char fsname10[100]="";
  char temp_string[1000]="";
  std::string string_addto = "";
  
  strcpy(fpname,"stem_parameters_");
  strcat(fpname, argv[1]);
  strcat(fpname, ".dat");
  
  strcpy(fsname1,"stem_solution_");
  strcat(fsname1, argv[1]);
  strcat(fsname1, "_1.dat");
  
  strcpy(fsname2,"stem_solution_");
  strcat(fsname2, argv[1]);
  strcat(fsname2, "_2.dat");
  
  strcpy(fsname3,"stem_solution_");
  strcat(fsname3, argv[1]);
  strcat(fsname3, "_3.dat");
  
  strcpy(fsname4,"stem_solution_");
  strcat(fsname4, argv[1]);
  strcat(fsname4, "_4.dat");
  
  strcpy(fsname5,"stem_solution_");
  strcat(fsname5, argv[1]);
  strcat(fsname5, "_5.dat");
  
  strcpy(fsname6,"stem_solution_");
  strcat(fsname6, argv[1]);
  strcat(fsname6, "_6.dat");
  
  strcpy(fsname7,"stem_solution_");
  strcat(fsname7, argv[1]);
  strcat(fsname7, "_7.dat");
  
  strcpy(fsname8,"stem_solution_");
  strcat(fsname8, argv[1]);
  strcat(fsname8, "_8.dat");
  
  strcpy(fsname9,"stem_solution_");
  strcat(fsname9, argv[1]);
  strcat(fsname9, "_9.dat");
  
  strcpy(fsname10,"stem_solution_");
  strcat(fsname10, argv[1]);
  strcat(fsname10, "_10.dat");
  
  MPI_File_open(MPI_COMM_WORLD,fpname,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_p);
  MPI_File_open(MPI_COMM_WORLD,fsname1,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_1);
  MPI_File_open(MPI_COMM_WORLD,fsname2,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_2);
  MPI_File_open(MPI_COMM_WORLD,fsname3,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_3);
  MPI_File_open(MPI_COMM_WORLD,fsname4,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_4);
  MPI_File_open(MPI_COMM_WORLD,fsname5,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_5);
  MPI_File_open(MPI_COMM_WORLD,fsname6,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_6);
  MPI_File_open(MPI_COMM_WORLD,fsname7,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_7);
  MPI_File_open(MPI_COMM_WORLD,fsname8,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_8);
  MPI_File_open(MPI_COMM_WORLD,fsname9,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_9);
  MPI_File_open(MPI_COMM_WORLD,fsname10,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_s_10);
  
  srand(my_rank*2);

  int istart = my_rank*(num_paras/comm_sz);
  int istop  = (my_rank+1)*(num_paras/comm_sz);
  if(my_rank==comm_sz-1){
    istop=num_paras;
   }
    
  for (i = istart; i < istop; i++){
    
    set_parameters(argv);

    for (j = 0; j < num_ode; j++){
        cnt_loop = stem (y_store, j, atoi(argv[1]));
    }

    cnt  = count_state(y_store, num_ode, thrd, soln); 
    
    cnt_store[cnt-1] = cnt_store[cnt-1] + 1;  
    
    snprintf(temp_string, sizeof(temp_string), "%d\t%d\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\n", i+1, cnt, ga, gb, gc, gd, ge, gf, gg, gI, rb, ru, ka, kb, kc, kd, ke, kf, kg, kh, kI, lamdaaa, lamdaha, lamdaea, lamdaga, lamdaab, lamdacb, lamdacc, lamdaec, lamdagc, lamdaId, lamdaed, lamdagd, lamdaee, lamdafe, lamdahe, lamdade, lamdaae, lamdace, lamdaef, lamdabg, lamdacg, lamdahg, lamdahI, a0a, h0a, e0a, g0a, a0b, c0b, c0c, e0c, g0c, I0d, e0d, g0d, e0e, f0e, h0e, d0e, a0e, c0e, e0f, b0g, c0g, h0g,h0I, naa, nha, nea, nga, nab, ncb, ncc, nec, ngc, nId, ned, ngd, nee, nfe, nhe, nde, nae, nce, nef, nbg, ncg,nhg,nhI);
    
    MPI_File_write_shared(f_p, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
    fflush(stdout);
   
    switch(cnt){
    case 1  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_1, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 2  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_2, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 3  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_3, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 4  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_4, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 5  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_5, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 6  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_6, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 7  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_7, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 8 :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_8, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 9  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_9, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
    case 10  :
        string_addto = std::to_string(i+1)+"\t";
        for (h = 0; h < cnt; h++){
            h2 = 1;
            while (h2 <= NEQN) {
		string_addto+= std::to_string(log2(soln[NEQN*h + h2 - 1]))+"\t";
                h2++;
            }
       }
       string_addto+="\n";
       strcpy(temp_string, string_addto.c_str());
       MPI_File_write_shared(f_s_10, temp_string, strlen(temp_string), MPI_CHAR, MPI_STATUS_IGNORE);
       fflush(stdout);
       break;
	}
	
    
  }

  //fprintf(f_p, "---------------Conclusion---------------\n");
  
  for (i = 0; i < 10; i++){
     printf("%d\t%d\n", i+1, cnt_store[i]);
     //fprintf(f_p, "%d\t%d\n", i+1, cnt_store[i]);
  }

  MPI_File_close(&f_p);
  MPI_File_close(&f_s_1);
  MPI_File_close(&f_s_2);
  MPI_File_close(&f_s_3);
  MPI_File_close(&f_s_4);
  MPI_File_close(&f_s_5);
  MPI_File_close(&f_s_6);
  MPI_File_close(&f_s_7);
  MPI_File_close(&f_s_8);
  MPI_File_close(&f_s_9);
  MPI_File_close(&f_s_10);
  
  free(y_store);
  free(soln);
  
  //# undef NEQN
  
  MPI_Finalize();//clean up mpi

  return 0;
}

/******************************************************************************/
void set_parameters (char **argv) 
{ 
  ka = 0.1 + (1.0-0.1)*randu();
  kb = 0.1 + (1.0-0.1)*randu();
  kc = 0.1 + (1.0-0.1)*randu();
  kd = 0.1 + (1.0-0.1)*randu();
  ke = 0.1 + (1.0-0.1)*randu();
  kf = 0.1 + (1.0-0.1)*randu();
  kg = 0.1 + (1.0-0.1)*randu();
  kh = 0.1 + (1.0-0.1)*randu();
  kI = 0.1 + (1.0-0.1)*randu();
  
  lamdaaa = (1.0 + (100.0-1.0)*randu());
  lamdaha = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdaea = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdaga = (1.0 + (100.0-1.0)*randu());
  lamdaab = (1.0 + (100.0-1.0)*randu());
  lamdacb = (1.0 + (100.0-1.0)*randu());
  lamdacc = (1.0 + (100.0-1.0)*randu());
  lamdaec = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdagc = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdaId = (1.0 + (100.0-1.0)*randu());
  lamdaed = (1.0 + (100.0-1.0)*randu());
  lamdagd = (1.0 + (100.0-1.0)*randu());
  lamdaee = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdafe = (1.0 + (100.0-1.0)*randu());
  lamdahe = (1.0 + (100.0-1.0)*randu());
  lamdade = (1.0 + (100.0-1.0)*randu());
  lamdaae = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdace = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdaef = (1.0 + (100.0-1.0)*randu());
  lamdabg = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdacg = 1.0/(1.0 + (100.0-1.0)*randu());
  lamdahI = (1.0 + (100.0-1.0)*randu());
  lamdahg = (1.0 + (100.0-1.0)*randu());
  lamdagh = (1.0 + (100.0-1.0)*randu());
  lamdaIh = (1.0 + (100.0-1.0)*randu());


  a0a = 0.777*(0.02 + (1.98-0.02)*randu());
  a0b = 0.777*(0.02 + (1.98-0.02)*randu());
  a0e = 0.777*(0.02 + (1.98-0.02)*randu());
  b0g = 13.782*(0.02 + (1.98-0.02)*randu());
  c0b = 2.302*(0.02 + (1.98-0.02)*randu());
  c0c = 2.302*(0.02 + (1.98-0.02)*randu());
  c0g = 2.302*(0.02 + (1.98-0.02)*randu());
  c0e = 2.302*(0.02 + (1.98-0.02)*randu());
  d0e = 3.933*(0.02 + (1.98-0.02)*randu());
  e0a = 0.096*(0.02 + (1.98-0.02)*randu());
  e0c = 0.096*(0.02 + (1.98-0.02)*randu());
  e0f = 0.096*(0.02 + (1.98-0.02)*randu());
  e0d = 0.096*(0.02 + (1.98-0.02)*randu());
  e0e = 0.096*(0.02 + (1.98-0.02)*randu());
  f0e = 41.360*(0.02 + (1.98-0.02)*randu());
  g0a = 0.533*(0.02 + (1.98-0.02)*randu());
  g0c = 0.533*(0.02 + (1.98-0.02)*randu());
  g0d = 0.533*(0.02 + (1.98-0.02)*randu());
  g0h = 0.533*(0.02 + (1.98-0.02)*randu());
  h0a = 2.074*(0.02 + (1.98-0.02)*randu());
  h0e = 2.074*(0.02 + (1.98-0.02)*randu());
  h0g = 2.074*(0.02 + (1.98-0.02)*randu());
  h0I = 2.074*(0.02 + (1.98-0.02)*randu());
  I0d = 4.590*(0.02 + (1.98-0.02)*randu());
  I0h = 4.590*(0.02 + (1.98-0.02)*randu());

  naa = randud();
  nha = randud();
  nea = randud();
  nga = randud();
  nab = randud();
  ncb = randud();
  ncc = randud();
  nec = randud();
  ngc = randud();
  nId = randud();
  ned = randud();
  ngd = randud();
  nee = randud();
  nfe = randud();
  nhe = randud();
  nde = randud();
  nae = randud();
  nce = randud();
  nef = randud();
  nbg = randud();
  ncg = randud();
  nhg = randud();
  nhI = randud();
  nIh = randud();
  ngh = randud();

  ga = (1.0 + (100.0-1.0)*randu())/(lamdaaa*lamdaga);
  gb = (1.0 + (100.0-1.0)*randu())/(lamdaab*lamdacb);
  gc = (1.0 + (100.0-1.0)*randu())/(lamdacc);
  gd = (1.0 + (100.0-1.0)*randu())/(lamdaId*lamdaed*lamdagd);
  ge = (1.0 + (100.0-1.0)*randu())/(lamdafe*lamdahe*lamdade);
  gf = (1.0 + (100.0-1.0)*randu())/(lamdaef);
  gg = (1.0 + (100.0-1.0)*randu())/(lamdahg);
  rb = (0.001 + (1.0-0.001)*randu());
  ru = (0.1 + (1.0-0.1)*randu());
  gI = (1.0 + (100.0-1.0)*randu())/(lamdahI);
  
  
}

/******************************************************************************/
int count_state (double y_store[], int num_ode, double thrd, double soln[])  // distance and no iteration
{
    int i = 0;
    int j = 0;
    int h = 0;
    int count = 0;
    int cnt = 1;
    double delta = 0.0;
    double sumpow = 0.0;
    
    for (h = 1; h <= NEQN; h++){
        soln[h - 1] = y_store[h - 1]; 
    }
    
    for (i = 2; i <= num_ode; i++){
        count = 0;
        
        for (j = 1; j <= cnt; j++){
            h = 1;
            sumpow = 0.0;
            while (h <= NEQN) {
                sumpow = sumpow + pow((y_store[NEQN*(i-1) + h - 1] - soln[NEQN*(j-1) + h - 1]), 2);
                h++;
            }
        
            delta = sqrt(sumpow);
            
            if (delta > thrd){
                count = count + 1;
            } 
        }    
        
        if (count == cnt){
            cnt = cnt + 1;
            
            if (cnt < 10) {
                for (h = 1; h <= NEQN; h++){
                    soln[(cnt-1)*NEQN + h - 1] = y_store[NEQN*(i-1) + h - 1]; 
                }
            }
            else{
                cnt = 10;
                break;
            }
        }
    }
    
    return cnt;
}

/******************************************************************************/
int stem (double y_store[], int j, int num)
{
  int    i_step = 1;
  int    n_step = 1000; 
  int    i = 0;
  double testdelta = 0.0;
  double t = 0.0;
  double dt = (double)100/n_step;
  double t_start = 0.0;
  double t_stop = 0.0;
  double y[NEQN]  = {0.0};;
  double yp[NEQN] = {0.0};;
  double ytmp[NEQN] = {2000.0};
 
  int cnt_loop = 0;
  
  y[0] = exp2(log2(ga*(lamdaha*lamdaea)/ka) 		   + (log2(ga*(lamdaaa*lamdaga)/ka) 		- log2(ga*(lamdaha*lamdaea)/ka))*randu());
  y[1] = exp2(log2(gb/kb) 				   		   + (log2(gb*(lamdaab*lamdacb)/kb) 			- log2(gb/kb))*randu());
  y[2] = exp2(log2(gc*(lamdaec*lamdagc)/kc) 		   + (log2(gc*(lamdacc)/kc) 				- log2(gc*(lamdaec*lamdagc)/kc))*randu());
  y[3] = exp2(log2(gd/kd) 				   		   + (log2(gd*(lamdaId*lamdaed*lamdagd)/kd) 	- log2(gd/kd))*randu());
  y[4] = exp2(log2(ge*(lamdaee*lamdaae*lamdace)/ke) + (log2(ge*(lamdafe*lamdahe*lamdade)/ke) 	- log2(ge*(lamdaee*lamdaae*lamdace)/ke))*randu());
  y[5] = exp2(log2(gf/kf) 						   + (log2(gf*(lamdaef)/kf) 					- log2(gf/kf))*randu());
  y[6] = exp2(log2(gg*(lamdabg*lamdacg)/kg) 		   + (log2(gg*lamdahg/kg) 				- log2(gg*(lamdabg*lamdacg)/kg))*randu());
  //y[7] = exp2(log2(rb/ru)*randu());
  y[7] = exp2(log2(gh/kh)*randu());
  y[8] = exp2(log2(gI/kI) 						   + (log2(gI*lamdahI/kI) 					- log2(gI/kI))*randu());
 
  testdelta = sumdelta(y, ytmp);
  
  while (testdelta != 0 && cnt_loop < 20) {
    t_start = t_stop;
    t_stop  = t_stop + 100;
    
    cnt_loop = cnt_loop + 1;
  
    for ( i_step = 1; i_step <= n_step; i_step++ )
    { 
        for (i = 0; i < NEQN; i++){
            ytmp[i] = y[i];
        }
        
		model ( t, ytmp, yp );
		
		for (i = 0; i < NEQN; i++){
            y[i] = ytmp[i] + yp[i]*dt;
        }
        
        t = t + dt;
        
    }
    
    testdelta = sumdelta(y, ytmp);
  }
  
  for (i = 0; i < NEQN; i++){
    y_store[NEQN*j + i] = y[i];
  }
  
  return cnt_loop;
}

/******************************************************************************/
double log2 (double x)
{
    return (log10(x)/log10(2));   
}

/******************************************************************************/
double abs_test (double x)
{
    if (x < 0)
        return -1.0*x;
    else
        return x;
}

/******************************************************************************/
double sumdelta (double y[], double ytmp[])
{
    int i = 0;
    double out = 0.0;
    
    for (i = 0; i < NEQN; i++){
        out = out + (y[i] - ytmp[i])*(y[i] - ytmp[i]);
    }

    return sqrt(out);
}

/******************************************************************************/
double Hillshift (double x, double x0, double nx, double lamda)
{
    double out;

    out = lamda + (1.0 - lamda) * (1.0/(1.0 + pow((x/x0), nx)));
    
    return out;
}

/******************************************************************************/
void model ( double t, double y[], double yp[] )
{
  yp[0] = ga*Hillshift(y[0], a0a, naa, lamdaaa)*Hillshift(y[7], h0a, nha, lamdaha)*Hillshift(y[4], e0a, nea, lamdaea)*Hillshift(y[6], g0a, nga, lamdaga)  - ka * y[0]; //A Gata6
  yp[1] = gb*Hillshift(y[0], a0b, nab, lamdaab)*Hillshift(y[2], c0b, ncb, lamdacb)  - kb * y[1]; //B Gcnf
  yp[2] = gc*Hillshift(y[2], c0c, ncc, lamdacc)*Hillshift(y[4], e0c, nec, lamdaec)*Hillshift(y[6], g0c, ngc, lamdagc) - kc * y[2];  //C Cdx2
  yp[3] = gd*Hillshift(y[8], I0d, nId, lamdaId)*Hillshift(y[4], e0d, ned, lamdaed)*Hillshift(y[6], g0d, ngd, lamdagd) - kd * y[3];  //D Klf4
  yp[4] = ge*Hillshift(y[4], e0e, nee, lamdaee)*Hillshift(y[5], f0e, nfe, lamdafe)*Hillshift(y[7], h0e, nhe, lamdahe)*Hillshift(y[3], d0e, nde, lamdade)*Hillshift(y[0], a0e, nae, lamdaae)*Hillshift(y[2], c0e, nce, lamdace) - ke * y[4];  //E Nanog
  yp[5] = gf*Hillshift(y[4], e0f, nef, lamdaef) - kf * y[5];  //F Pbx1
  yp[6] = gg*Hillshift(y[7], h0g, nhg, lamdahg)*Hillshift(y[1], b0g, nbg, lamdabg)*Hillshift(y[2], c0g, ncg, lamdacg) - kg * y[6]- rb*y[6]*y[8]+ ru * y[7]; //G Oct4
  yp[7] = rb*y[6]*y[8]-ru*y[7];  //H Oct4-Sox2
  yp[8] = gI*Hillshift(y[7], h0I, nhI, lamdahI) - kI*y[8]-rb*y[6]*y[8]+ru*y[7];  //I Sox2

  return;
}

/******************************************************************************/
double randu(void)  // uniform random number between 0 and 1
{
    double u;
	u = ((double)rand()/RAND_MAX);	

	return u;	
}

/******************************************************************************/
double randud(void)  // discrete integer uniform random number between 0 and 6
{
    double u = 0.0;
    double z = 0.0;		
	
	do {
		u = 6.0*((double)rand()/RAND_MAX);
	
		if (u < 1){
	    	z = 1.0;
		}
		else if ( u >= 1 && u < 2){
		    z = 2.0;
		}
		else if ( u >= 2 && u < 3){
		    z = 3.0;
		}
		else if ( u >= 3 && u < 4){
		    z = 4.0;
		} 
		else if ( u >= 4 && u < 5){
		    z = 5.0;
		}
		else if ( u >= 5 && u < 6){
		    z = 6.0;
		}
	} while (u == 6.0);

	return z;	
}
