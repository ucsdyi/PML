#include <time.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>
using namespace std;


#define MaxLPhi 5000000
#define MaxLP 2000000

#define SzRandAr 3000000//10000000
#define SPcond 100000000

#define NaNlog -100
#define SAVETEMP "20"
#define PML_FILENAME "PMLFile"
#define PROFILE_FILENAME "proFile"
#define INITIAL_DIST_FILENAME ""//initial_dist"
#define MinOp 500//50
//#define MaxEM 20
//#define MAXSZ 10000
//#define SRPLNaNLog -1e6

#define SHRINK_EM 1

double estCondlog(int, int, int, int, int, int, int,    double *, int *,  double *, double *);
int maxinitialize(int, int, int, int, int, double *, int *, int, double *);
int randominitialize(int, int, int, int, int, double *, int *, int,  double *);
void mergesort(int, double *);
void mergesortdouble(int, int *, int *);
double logsum(double, double);
double  mpd(int *, int, double *, int, char *, int);


int PML(int MAXSZ=10000, int maximum_EM=20, int EM_n=100)
{
    /*
    FILE *mpf;
    mpf=fopen("otherOutputs","w");
    if(mpf==NULL){
       //printf("FILE OPEN ERROR on otherOutputs");
        exit(1);
    }
    */
    
    int *phi=(int *)malloc((MaxLPhi+1)*sizeof(int));
    if(phi==NULL){
       //printf("Unable to malloc for phi[]\n");
        exit(1);
    }
    for(int i=0;i<=MaxLPhi;i++){
        phi[i]=0;
    }
    
    char matlabFile[100] =PML_FILENAME;
    char proFile[100] = PROFILE_FILENAME;

    
    // parameters for saving temporary PML
    char *tempFilenamePrefix = (char*) malloc(50*sizeof(char));
    sprintf(tempFilenamePrefix,"%s",matlabFile);
    
    char str_SAVETEMP[100] = SAVETEMP;
    
    // -----------------------------------
    
    FILE *toMatlab;
    
    //proFile contains phi[1],phi[2],...
    FILE *readprofile;
    readprofile = fopen(proFile,"r");
    if(readprofile==NULL){
       //printf("Can't open file aab!");
        exit(1);
    }
    int profile_length=0;
    int next_prevalence;
    for(profile_length=0;;){
        int ret= fscanf(readprofile,"%d", &next_prevalence);
        if(ret==EOF)  break;
        profile_length++;
        
        phi[profile_length]=next_prevalence;
    }
    //printf("profile = %d\n",phi[0]);
    phi[0]=profile_length;
    fclose(readprofile);

    
    
    //Writing profile to Output_File
    /*fprintf(mpf,"profile length: %d\nprofile:",phi[0]);
    for(int i=1;i<=phi[0];i++)
        fprintf(mpf," %d",phi[i]);
    fprintf(mpf,"\n");
    */
    
    //Display profile on screen
   //printf("------------------------------\n");
   //printf("profile length: %d\nprofile:",phi[0]);
    int headPart = 10;
    int tailPart = phi[0]-10;
    if (headPart>phi[0]) headPart = phi[0];
    if(tailPart <= headPart) tailPart = headPart+1;
    
    //for(int i=1;i<= headPart;i++)
       //printf(" %d",phi[i]);
    //if(headPart < phi[0])//printf(" ...");
    //for(int i=tailPart;i<=phi[0];i++)
       //printf(" %d",phi[i]);
    
   //printf("\n");
   //printf("------------------------------\n");
    //================================================
    
    
    //the program begins =============================
    int k=0,n=0;
    for (int i=1;i<=phi[0];i++){
        k+=phi[i];
        n+=i*phi[i];
    }
   //printf("n=%d,k=%d\n",n,k);
    
    //int MAXSZ = floor((k * log(n)) + 0.5); //***
    
    //upper bound on MPPlog
    double uMPPlog=0.0f;
    for(int i=1;i<=n;i++){
        uMPPlog+=log((double)i);
    }
    for( int i=1;i<=phi[0];i++){
        double ilog=0.0f;
        
        for(int j=1;j<=i;j++){
            ilog+=log((double)j);
        }
        double philog=0.0f;
        if(phi[i]>0){
            for(int j=1;j<=phi[i];j++){
                philog+=log((double)j);
            }
        }
        uMPPlog-=philog+phi[i]*ilog;
    }
   //printf("uMMlog=%.4e;\n", -uMPPlog);
    
    //fprintf(mpf, "uMMPlog=%.4e;\n",-uMPPlog);
    int pl[20];
    double errphi[20];
    double MPPlog[20];
    //    int i;
    for(int i=1; i<2; i++)
    {
        pl[i]=MAXSZ;
        double *plog=(double *)malloc((MaxLP+1)*sizeof(double));
        
        for(int j=0;j<=MaxLP;j++){
            plog[j]=NaNlog;
        }
        if(pl[i]> MaxLP)
        {
            pl[i] = MaxLP; //printf("Initial support size too large, reset to MaxLP\n");
        }
        
        //------------------------------------------------------------
        char initiald[100]= INITIAL_DIST_FILENAME;
        if(strlen(initiald) > 0) {
            FILE *initd;
            initd  = fopen(initiald, "r");
            if(initd==NULL){
                printf("Can't open file!");
                exit(1);}
            int prob_length = 0;
            float next_prob;
            for (prob_length = 0;;)
            {
                int ret  = fscanf(initd,"%f", &next_prob);
                if(ret==EOF) break;
                prob_length++;
                plog[prob_length]=log(next_prob);
            }
            fclose(initd);
            printf("%d probabilities read\n", prob_length);
            
            if(pl[i] > prob_length){
                for (int j = prob_length+1;j<=pl[i];j++){
                    plog[j]=NaNlog;
                }
            }
            else{
                double tailsum = 0;
                for(int j=pl[i]+1;j<=prob_length; j++){
                    tailsum += exp(plog[j]);
                }
                tailsum = log(1 - tailsum);
                for(int j=1; j<= pl[i];j++){
                    plog[j] -= tailsum;
                }
            printf("Warning: extra probabilities truncated (the rest normalized).\n");
                //pl[i] = prob_length; // reset initial support size
            }
        }
        else{
            //printf("the initial distribution is uniform with small random disturbations");
            for(int j=1; j<=pl[i]; j++){
                plog[j]=log(1.0/pl[i]*(1+0.01*rand()/RAND_MAX));
            }
            
        }
        //------------------------------------------------------------
       /*printf("Initial support size: %d\n", pl[i]);
        //print the intial distribution
       printf("Initial distribution:\n");
        for(int j=0;j<=5;j++)
        {
        printf("p[%d]=%.5f\n",j,exp(plog[j]));
        }
       printf("...\n");*/
        
        
        MPPlog[i]=mpd(phi,pl[i],plog,maximum_EM,tempFilenamePrefix,EM_n);
       //printf("mpd(..) done!\n");
       // printf(MPPlog[i]);
        // errphi[i]=monitorphi(phi,n,pl[i],plog);
        errphi[i]=0;
        
       //printf("pl(%d)=%d;\n", i+1, pl[i]);
       //printf("mpplog(%d)=%f;\n",i+1, MPPlog[i]);
       //printf("errphi(%d)=%f;\n",i+1, errphi[i]);
        
        
        //fprintf(mpf,"pl(%d)=%d;\n", i+1, pl[i]);
        //fprintf(mpf,"mpplog(%d)=%f;\n",i+1, MPPlog[i]);
        //fprintf(mpf,"errphi(%d)=%f;\n",i+1, errphi[i]);
        
        
        //Open file for writing distribution for use in Matlab
        //toMatlab = 'aaaa';
        toMatlab=fopen(matlabFile,"w");
        if(toMatlab==NULL){
           //printf("FILE OPEN ERROR on toMatlab");
            exit(1);
        }
        
        for( int j=1;j<=pl[i];j++){
            double next_prob;
            if(plog[j]>NaNlog)
            {
                next_prob = exp(plog[j]);
                fprintf(toMatlab,"%.4e\n",next_prob);
            }
            else
                break;
        }
        
        free(plog);
    }
    
    free(tempFilenamePrefix);
    //fclose(mpf);
    fclose(toMatlab);
    //====================================
    free(phi);
    
    
   //printf("\nFiles:\n");
   //printf("(1) Profile output to \"%s\"\n",proFile);
   //printf("(2) PML distribution output to \"%s\"\n",
    //       matlabFile);
   //printf("(3) Some other information output to \"otherOutputs\"\n");

    return 0;
}

extern "C"{
    int PMLplus(int argv0, int argv1, int argv2){
        return PML(argv0, argv1, argv2);}
}


int  randominitialize(int BWantPc, int phio, int k, int pl, int apl, double *plog, int *f, int ipcinphio, double *epliolog)
{
    if(pl<k){
       //printf("\n pl=%d<k=%d not implemented", pl, k);
        exit(1);
    }
    int *lf=(int *)malloc((apl+1)*sizeof(int));
    int *rf=(int *)malloc((apl+1)*sizeof(int));
    int *randAr=(int *)malloc((apl+1)*sizeof(int));
    if(lf==NULL||randAr==NULL){
       //printf("Unable to malloc for randAr (randominitialize)\n");
        exit(1);
    }
    for(int i=1; i<=apl; i++){
        lf[i]=i;
        rf[i]=0;
        randAr[i]=rand();
    }
    int pcinphio=ipcinphio;
    if(BWantPc==0){
        mergesortdouble(pl, randAr+1, lf+1);
        for(int i=1;i<=pl;i++){
            f[i]=lf[i];
        }
        pcinphio=0;
    }
    else{
        mergesortdouble(apl, randAr+1, lf+1);
        for(int i=1;i<=phio;i++){
            f[i]=lf[i];
            rf[lf[i]]=1;
            if(f[i]>pl){
                pcinphio++;
            }
        }
        int prf=1;
        int restpl=pl-phio;
        for(int i=1; i<=restpl; i++){
            while(rf[prf]!=0){ prf++;}
            lf[i]=prf;
            prf++;
        }
        for(int i=1; i<=restpl; i++){
            randAr[i]=rand();
        }
        if(restpl>0){
            mergesortdouble(restpl, randAr+1, lf+1);
        }
        for(int i=1; i<=restpl; i++){
            f[phio+i]=lf[i];
            rf[lf[i]]=1;
        }
        prf=1;
        for(int i=pl+1; i<=apl; i++){
            while(rf[prf]!=0){
                prf++;
            }
            f[i]=prf;
            prf++;
        }
    }//else
    
    free(lf);
    free(rf);
    free(randAr);
    
    return(pcinphio);
}

double  mpd(int * phi, int pl, double *plog, int maximum_EM, char *tempFilenamePrefix, int EM_n)
//Subroutine for searching maximum pattern dist(MPD).
//The return value is the maximum pattern prob(MPP).
//phi is the input profile and phi[0] is the profile length.
//pl is the prob. length or to say the assumed support size
//p is the intial prob. in descending order carries the returned MPD.
//p[0] is the possible continous part
//p[1] is largest discrete prob. and so on.
{
    
    int num_saved = 0;
    double rMPPlog=0;
    
    //Phil, phi length
    int phil=phi[0];
    
    plog[0]=NaNlog;
    //k is the number of appeared elements and n is the pattern length.
    int k=0,n=0;
    for (int i=1;i<=phil;i++){
        k+=phi[i];
        n+=i*phi[i];
    }
    
    //Determine whether use continuous part
    int BWantPc=0;

    //prob. length dl may be less than initial profle length ipl
    //extra prob. length(>phi[1]) is for injection for p[0]
    //apl(All prob. length) is the dynamic prob. length plus possible extra length.
    int epl=phi[1];
    int apl=pl+BWantPc*epl;
    
    if (phi[0]>MaxLPhi){
       //printf("The number of appeared elements %d excees MaxPhi %d\n", k, MaxLPhi);
        return(0);
    }
    if (apl>MaxLP){
       //printf("All profile length %d excees MaxLP %d\n", apl, MaxLP);
        return(0);
    }
    
    //the array for accumulator cp[]
    double* cp=(double *)malloc((apl+1)*sizeof(double));
    if(cp==NULL) {printf("Unable to malloc for cp[].\n"); exit(1);}
    
    //the array for last visted (time) v[]
    double* lv=(double *)malloc((apl+1)*sizeof(double));
    if(lv==NULL) {printf("Unable to malloc for lv[].\n"); exit(1);}
    
    //the array for injection mapping f[]
    int* f=(int *)malloc((apl+1)*sizeof(int));
    if(f==NULL) {printf("Unable to malloc for f[].\n"); exit(1);}
    
    //the array for multificity mu[]
    int* mu=(int *)malloc((apl+1)*sizeof(int));
    if(mu==NULL) {printf("Unable to malloc for mu[].\n"); exit(1);}
    
    //and the array for logeplio[], log(epl-i+1)
    double* epliolog=(double *)malloc((epl+1)*sizeof(double));
    if(epliolog==NULL) {printf("Unable to malloc for epliolog[].\n"); exit(1);}
    
    //malloc random array 1, 2 and log
    long* randAr1=(long *)malloc((SzRandAr+1)*sizeof(long));
    if(randAr1==NULL) {printf("Unable to malloc for randAr1[].\n"); exit(1);}
    
    long* randAr2=(long *)malloc((SzRandAr+1)*sizeof(long));
    if(randAr2==NULL) {printf("Unable to malloc for randAr2[].\n"); exit(1);}
    
    double* randArlog=(double *)malloc((SzRandAr+1)*sizeof(double));
    if(randArlog==NULL) {printf("Unable to malloc for randArlog[].\n"); exit(1);}
    
    
    if (lv==NULL||cp==NULL||f==NULL||mu==NULL||epliolog==NULL
        ||randAr1==NULL||randAr2==NULL||randArlog==NULL){
       //printf("\n malloc error");
        exit(1);
    }
    
    
    //set mu[]
    int s=0;
    
    for (int i=1;i<=phi[0];i++){
        for(int j=0;j<phi[i];j++){
            s++;
            mu[s]=i;
        }
    }
    
    //set extra mu to zero
    for (int i=k+1;i<=apl;i++){
        mu[i]=0;
    }
    
    //Set the random array
    double minrand=exp(NaNlog)*RAND_MAX;
    
    // generate accept probabilities
    for(int i=1;i<=SzRandAr;i++){
        int j;
        do{ j=rand();
        }while((double)j<minrand);
        randArlog[i]=log((double)j/(double)RAND_MAX);
    }
    
    
    //Running EM algorithm untill MaxEM iterations
    printf("EM algorithm for n=%d begins ... ", EM_n);
    int EM;
    int sz = pl;
    
    /*
    FILE *num_em;
    num_em=fopen("EM","w");
    fprintf(num_em,"0\n");
    fclose(num_em);
    */
    
    for(EM=1;EM<=maximum_EM;EM++){
        // set probability length to support size
        // i.e., removing zero probabilities
        pl=sz;
        apl=sz + BWantPc*epl;
        
        if (EM % SHRINK_EM == 0) {
            // Generate a random neighbor state
            for(int i=1;i<=SzRandAr;i++){
                randAr1[i] = rand() % k + 1;
                randAr2[i] = rand() % apl + 1;
            }
        }
        
        //clean cp[] and lv[]
        for (int i=0;i<=apl;i++){
            cp[i]=0;
            lv[i]=0;
        }
        
        int vMinOp=MinOp;
        //accelerate in the first two iterations
        /*if(EM<=2){   ???
            vMinOp/=10;
        }*/
        
        int pcinphio;
        //initlize f[] and possible nonzero pcinphio
        if(pl>=k){
            pcinphio=randominitialize(BWantPc,phi[1], k, pl, apl, plog, f, 0, epliolog);
        }
        else{
            pcinphio=maxinitialize(BWantPc, phi[1], k, pl, apl, plog, f,0, epliolog);
        }
        
        //Time in sampling
        double t=0.0;
        
        //Run Markov chain sampling between Min
        do{
            //set pointers for random arrayes.
            long *pAr1=randAr1;
            long *pAr2=randAr2;
            double *pArlog=randArlog;
            
            //Run MCsampling to collect pf
            for(int i=1;i<=SzRandAr;i++){
                
                t+=1.0;
                pAr1++;
                pAr2++;
                pArlog++;
                
                double rlog=*pArlog;
                
                //Jump to next state
                long j1 = randAr1[i];  //changed
                long j2 = randAr2[i];  //changed
                while( j2 == j1)
                    j2 = rand () % apl + 1;
                
                //Calculate reject prob and test
                //if fail, then stay.
                double dplog=plog[f[j2]]-plog[f[j1]];
                double dmu=mu[j1]-mu[j2];
                double dprodlog=dplog*dmu;

                if(rlog>dprodlog){
                    continue;
                }
                
                //accumulate cp
                cp[f[j1]]+=(t-lv[f[j1]])*mu[j1];
                if(mu[j2]>0){
                    cp[f[j2]]+=(t-lv[f[j2]])*mu[j2];
                }
                
                //swapping f[j1] and f[j2];
                int swp=f[j1];
                f[j1]=f[j2];
                f[j2]=swp;
                
                //setting lv[f[j1]] and lv[f[j2]]
                lv[f[j1]]=t;
                lv[f[j2]]=t;
            }
            
        }while(t<vMinOp);
        
        //sum the rest cp
        for(int i=1;i<=apl;i++){
            cp[f[i]]+=(double)(t-lv[f[i]])*(double)mu[i];
        }
        
        //sum cp[pl+1]+...+cp[apl] to get cp[0];
        cp[0]=0;
        
        //calculate plog
        double ntlog=log(n)+log(t);
        
        for(int i=0;i<=pl;i++){
            if(cp[i]==0){
                plog[i]=NaNlog;
            }
            else{
                plog[i]=log(cp[i]);
            }
            plog[i]-=ntlog;
        }
        
        //merge sort plog[1:1:pl]
       //printf("\nBegin mergesorting...");
        mergesort(pl,plog+1);
       //printf("Mergesorting done!\n");
        
        double minp = NaNlog;
        int last_sz = sz;
        for( int i=pl; i>=k; i--){
            if(plog[i] > NaNlog){
                sz=i;
                minp = exp(plog[i]);
                break;
            }
        }
        
        if (EM==1)
            std::cout << "[=";
        else if (EM==maximum_EM)
            std::cout << "=]\n";
        else
            std::cout << "=";
        std::cout.flush();
     /*
        num_em=fopen("EM","w");
        fprintf(num_em,"%d\n",EM);
        fflush(num_em);
        fclose(num_em);
     */
    }//for EM
    
    /**/
    //fclose(MPC);
    free(cp);
    free(lv);
    free(f);
    free(mu);
    free(epliolog);
    free(randAr1);
    free(randAr2);
    free(randArlog);
    return(rMPPlog);
}

void mergesortdouble(int m, int * randAr, int *lf)
{
    int n=1<<((int)(log(m)/log(2.0))+1);
    void submerge(int, int *brandAr, int *blf);
    int * brandAr=(int *)malloc(n*sizeof(int));
    int * blf=(int *)malloc(n*sizeof(int));
    if (brandAr==NULL||blf==NULL){
       //printf("Unable to malloc for blf (mergesortdouble)\n");
        exit(1);
    }
    for(int i=0; i<m; i++){
        brandAr[i]=randAr[i];
        blf[i]=lf[i];
    }
    for( int i=m; i<n; i++){
        brandAr[i]=-1;
        blf[i]=-1;
    }
    
    submerge(n, brandAr, blf);
    for(int i=0; i<m; i++){
        randAr[i]=brandAr[i];
        lf[i]=blf[i];
    }
    free(brandAr);
    free(blf);
    return;
}

void submerge(int n, int *randAr, int *lf)
{
    int * brandAr=(int *)malloc(n*sizeof(int));
    int * blf=(int *)malloc(n*sizeof(int));
    if (brandAr==NULL||blf==NULL){
       //printf("Unable to malloc for blf (submerge)\n");
        exit(1);
    }
    for(int i=0; i<n; i++){
        brandAr[i]=randAr[i];
        blf[i]=lf[i];
    }
    int hn=n/2;
    int *brandArh=brandAr+hn;
    int *blfh=blf+hn;
    if (hn>=2){
        submerge(hn, brandAr, blf);
        submerge(hn, brandArh, blfh);
    }
    int i1=0;
    int i2=0;
    int i=0;
    while(i1<hn&&i2<hn){
        if(brandAr[i1]>brandArh[i2]){
            randAr[i]=brandAr[i1];
            lf[i]=blf[i1];
            i1++;
        }
        else{
            randAr[i]=brandArh[i2];
            lf[i]=blfh[i2];
            i2++;
        }
        i++;
    }
    if(i1<hn){
        while(i1<hn){
            randAr[i]=brandAr[i1];
            lf[i]=blf[i1];
            i++;
            i1++;
        }
    }
    else{
        while(i2<hn){
            randAr[i]=brandArh[i2];
            lf[i]=blfh[i2];
            i++;
            i2++;
        }
    }//else
    free(brandAr);
    free(blf);
    return;
}

void mergesort(int pl, double *plog){
    
    void submerge(int, double *);
    int n=1<<((int)(log(pl)/log(2.0))+1);
    double * bplog=(double *)malloc(n*sizeof(double));
    if (bplog==NULL){
       //printf("Unable to malloc for bplog (mergesort).\n");
        exit(1);
    }
    for(int i=0; i<pl; i++){
        bplog[i]=plog[i];
    }
    for(int i=pl; i<n; i++){
        bplog[i]=NaNlog;
    }
    
    submerge(n, bplog);
    for(int i=0; i<pl; i++){
        plog[i]=bplog[i];
    }
    free(bplog);
    return;
}
void submerge(int n, double *plog)
{
    double * bplog=(double *)malloc(n*sizeof(double));
    if (bplog==NULL){
       //printf("Unable to malloc for bplog (submerge)");
        exit(1);
    }
    for(int i=0; i<n; i++){
        bplog[i]=plog[i];
    }
    int hn=n/2;
    double *bplogh=bplog+hn;
    if (hn>=2){
        submerge(hn, bplog);
        submerge(hn, bplogh);
    }
    int i1=0;
    int i2=0;
    int i=0;
    while(i1<hn&&i2<hn){
        if(bplog[i1]>bplogh[i2]){
            plog[i]=bplog[i1];
            i1++;
        }
        else{
            plog[i]=bplogh[i2];
            i2++;
        }
        i++;
    }
    if(i1<hn){
        while(i1<hn){
            plog[i]=bplog[i1];
            i++;
            i1++;
        }
    }
    else{
        while(i2<hn){
            plog[i]=bplogh[i2];
            i++;
            i2++;
        }
    }//else
    free(bplog);
    return;
}

int maxinitialize(int BWantPc, int phio, int k, int pl, int apl, double *plog, int *f, int ipcinphio,double *epliolog)
{
    int pcinphio=0;
    if(pl>=k){
        for(int i=k;i>=1;i--){
            f[i]=k+1-i;
        }
        for(int i=k+1;i<=apl;i++){
            f[i]=i;
        }
        pcinphio=ipcinphio;
        if(BWantPc==1){
            for(int i=1;i<=phio;i++){
                if(plog[f[i]]<plog[0]-epliolog[pcinphio+1])
                {
                    int tmp=f[i];
                    f[i]=f[pl+i];
                    f[pl+i]=tmp;
                    pcinphio++;
                }
                else{
                    break;
                }
            }//for
        }//if BWantPc
    }//if
    if(pl<k){
        for(int i=k;i>=k-pl+1;i--){
            f[i]=(k-i)+1;
        }
        for(int i=k-pl;i>=1;i--){
            f[i]=k-i+1;
        }
        for(int i=k+1;i<=apl;i++){
            f[i]=i;
        }
        pcinphio=k-pl;
        for(int i=k-pl+1;i<=phio;i++){
            if(plog[f[i]]<plog[0]-epliolog[i]){
                int tmp=f[i];
                f[i]=f[pl+i];
                f[pl+i]=tmp;
                pcinphio++;
            }
        }
    }//for
    return(pcinphio);
}

double logsum(double a, double b)
{
    double c;
    double d=a-b;
    if(d>0){
        c=a;
        d=-d;
    }
    else{
        c=b;
    }
    if(d<NaNlog){
        return(c);
    }
    else{
        c+=log(1+exp(d));
        return(c);
    }
}

double estCondlog(int phio, int k, int fk, int pl, int apl, int BWantPc, int ipcinphio,
                  double * plog, int * mu,  double * epliolog, double *randArlog)
{
    if(k==1&&pl==1&&BWantPc==0){
        return(0.0);
    }
    if(pl==0&&BWantPc==1){
        return(log(1/(double)apl));
    }
    
    //the array for injection mapping f[]
    int* f=(int *)malloc((apl+1)*sizeof(int));
    if(f==NULL){
       //printf("\n Malloc Error");
    }
    //initlize f[]
    int pcinphio;
    if(pl>=k ){
        pcinphio=randominitialize(BWantPc, phio, k, pl, apl, plog, f, ipcinphio, epliolog);
    }
    else{
        pcinphio=maxinitialize(BWantPc, phio, k, pl, apl, plog, f, ipcinphio, epliolog);
    }
    
    //Time in sampling
    double t=0.0;
    
    double ft=0.0;
    
    //Run Markov chain sampling between Min
    do{
        //set pointers for random arrayes.
        //Run MCsampling to collect pf
        double *pArlog=randArlog;
        for(int i=1; i<=SzRandAr; i++){
            t+=1;
            if(fk==1&&f[k]==1){
                ft++;
            }
            if(fk>pl&&f[k]>pl){
                ft++;
            }
            
            //Jump to next state
            int j1;
            int j2;
            
            j1 = rand() % k + 1;
            do{
                j2 = rand() % apl + 1;
            } while(j2 == j1);
            
            //------------------------------
            
            pArlog++;
            double rlog=*pArlog;
            
            
            int dpcinphio=0;
            if(BWantPc==1){
                dpcinphio=(f[j2]>pl)-(f[j1]>pl)
                +(j2<=k)*((f[j1]>pl)-(f[j2]>pl));
            }
            
            //If going to assign u>=2 to pc, then stay.
            if (BWantPc==1){
                if((j1>phio)&&(f[j2]>pl)){
                    continue;
                }
                if((j2>phio)&&(j2<=k)&&(f[j1]>pl)){
                    continue;
                }
            }
            
            //Calculate reject prob and test
            //if fail, then stay.
            double dplog=plog[f[j2]]-plog[f[j1]];
            double dmu=mu[j1]-mu[j2];
            double dprodlog=dplog*dmu;
            if(BWantPc==1){
                if(dpcinphio==1){
                    dprodlog-=epliolog[pcinphio+1];
                }
                else{
                    if(dpcinphio==-1){
                        dprodlog+=epliolog[pcinphio];
                    }
                }//else
            }//if==1
            if(rlog>dprodlog){
                continue;
            }
            
            //swapping f[j1] and f[j2];
            int swp=f[j1];
            f[j1]=f[j2];
            f[j2]=swp;
            
            
            //setting new pcinphio
            pcinphio+=dpcinphio;
            
        }//for i=1:1:SzRandAr
        
    }while(t<SPcond);
    
    
    
    double Pcond=(double)ft/(double)t;
    
   //printf(" %d",k);
    free(f);
    return(log(Pcond));
}
