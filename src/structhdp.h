#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<set>
#include<ctime>
#include<cstdlib>
#include<map>
#include<climits>
#include<vector>
#include<algorithm>
#include<sstream>
#include<iomanip>

using namespace std;

extern const int INIT_ALLOCATED_POPS, MAX_CENTROIDS, MAX_ITER, INIT_NUM_TOPICS, MAXALLELES, MAXK;
extern const double THRESH;
extern int PRINTCT,SATPRINTCT, MAX_GIBBS_ITER, BURNIN, INTERVAL,STIRLING_SIZE, NUMINDS, NUMLOCI, PLOIDY;
extern float GAMMA_A,GAMMA_B,ALPHA_A,ALPHA_B,H;
#define DEBUG 1
//#define DEBUG_0 1
extern long SEED;
extern string DATAFILE,OTHERPARAMS,OUTPUT_PREFIX,STIRLINGFILE;
extern map<string, int> gibbs_sampling_params;
extern map<string, float> gibbs_probability_params;

struct datamatrix{
	int ***data,***datacopy; //data will be destructively overwritten using the uniqalleles map for ease of access
	int ninds,nloci,ploidy;
};

struct model{
	int nloci,ploidy;
	int npops;
	double alpha_a,alpha_b,gamma_a,gamma_b;
	double gamma,alpha,h;
	int* numcentroids; //dim = nloci
	map<int,int>* uniqalleles; //stores unique alleles at each locus to reduce computation
	double** stirlingmatrix;
};

struct satellite{
	int ninds,nloci,npops,ploidy;
	int allocated_pops;
	//int nddotdot; // equal to nloci*ploidy - could be made to count only non-missing entries
	vector<int> njdotdot; // count only non-missing entries
	int* numcentroids;
	int** zmatrix;
	vector< vector< vector<int> > > count_lwk;
	vector< vector<int> > njdotk;
	vector< vector<int> > mjk;
	vector<int> mdotk;
	vector<int> mjdot;
	int mdotdot;
	vector<double> betas;
};


//data is ninds*nloci*ploidy
//double log1p(double);

/*model* initialise_model(int*** data,int nloci,int npops,int ploidy,int ninds);
void putinitialvalues(model* m,int*** data,int ninds);
void infermodel(model* m,int*** data,int ininds,double**** varc,double**** varz,double** vartheta);
void updatemodel(model* m,int*** data,int ninds,double**** varc,double**** varz,double** vartheta);
double computelogf(int x, int mu, double delta);
void updateparams(model* m,satellite* info);
void learnmodel(string filename);*/

/*satellite* init_from_model(model* m,int numinds,int*** data);
satellite* makecopy(satellite* current);*/
/*double sumlogs(double a, double b);
void printmodel(model* m,string s);
void dumpmodel(model* m,string s);
void printsatellite(satellite* info,string s);
void dumpsatellite(satellite* info,string s);*/

double logfcomputeold(int allele,int locus,int pop,satellite* info);
/*double zeroonerand();
double rangedrand(double min,double max);
double logit(double x);
double sigmoid(double x);
double psi(double x);*/
datamatrix* readdata(string filename);
