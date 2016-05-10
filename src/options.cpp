#include <getopt.h>
#include "structhdp.h"
#include "options.h"

void print_usage(){
	cout<<"Usage: structhdp uses the following input parameters in any order:\n";
	cout<<"-d <datafile> \n";
	cout<<"-n <number of individuals>\n";
	cout<<"-m <number of loci>\n";
	cout<<"-p <ploidy>\n";
	cout<<"-g <other-settings-file>\n";
	cout<<"-o <output-prefix> (optional, default='structhdp')\n";
	cout<<"-r <random seed> (optional, default=1)\n";
}

void valuegiven(int count,int len,char* arr[]){
        if(count+1==len){
                cout<<"No value provided for switch "<<arr[count]<<endl;
                exit(1);
        }
        else if(arr[count+1][0]=='-'){
                cout<<"No value provided for switch "<<arr[count]<<endl;
                exit(1);
        }
}

void verifyoption(int option,string message,string optswitch){
        if(option==-1){
                cout << "Must provide "<<message<<" with switch "<<optswitch<<endl;
                cout << "For help, use the switch \" -h\"\n";
                exit(1);
        }
}

void readopts(int len,char* arr[]){
        int opt=0,long_index=0;
        string temp;
        int dgiven=-1,ngiven=-1,ogiven=-1,pgiven=-1,rgiven=-1,ggiven=-1,mgiven=-1;
	SEED=1;
	OUTPUT_PREFIX="structhdp";
	static struct option long_options[]=
		{
		{"data_input",required_argument,0,'d'},
		{"output_prefix",required_argument,0,'o'},
		{"numinds",required_argument,0,'n'},
		{"numloci",required_argument,0,'m'},
		{"ploidy",required_argument,0,'p'},
		{"help",no_argument,0,'h'},
		{"random_seed",required_argument,0,'r'},
		{"settings_file",required_argument,0,'g'},
		{0,0,0,0}
		};
	while ((opt = getopt_long(len, arr,"d:o:n:m:p:hr:g:", 
                   long_options, &long_index )) != -1) {
        	switch (opt) {
                        case 'h':
                                //help option
				print_usage();
                                exit(0);
                        case 'd':
                                // Data file
                                DATAFILE=string(optarg); dgiven=1;
                                break;
                        case 'o':
                                //Output prefix
                                OUTPUT_PREFIX=string(optarg); ogiven=1;
                                break;
                        case 'n':
                                // No. of individuals
                                NUMINDS=atoi(optarg); ngiven=1;
                                break;
                        case 'm':
                                // No. of loci
                                NUMLOCI=atoi(optarg); mgiven=1;
                                break;
                        case 'r':
                                // Random seed
                                SEED=atoi(optarg); rgiven=1;
                                break;
                        case 'p':
                                // Ploidy
                                PLOIDY=atoi(optarg); pgiven=1;
                                break;
                                //break;
                        case 'g':
                                // Global params
                                OTHERPARAMS=string(optarg); ggiven=1;
                                break;
			case '?':
                        default:
                                cout<<"Invalid input"<<endl;
                                exit(1);
        	}
    	}
	if (optind < len) {
		cout<<"Non-option elements in input\n";
		exit(EXIT_FAILURE);
	}
        verifyoption(dgiven,"Datafile","-d");
        verifyoption(ngiven,"Number of individuals","-n");
        verifyoption(mgiven,"Number of loci","-m");
        verifyoption(pgiven,"Ploidy","-p");
        verifyoption(ggiven,"Other parameter file","-g");
        if(rgiven==-1)
                cout<<"No seed provided, using default value: "<<SEED<<endl;
        if(ogiven==-1)
                cout<<"No output prefix provided, using default value: "<<OUTPUT_PREFIX<<endl;

}
