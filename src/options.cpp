#include <getopt.h>
#include "structhdp.h"
#include "options.h"

void set_param_values(){
	MAX_GIBBS_ITER=gibbs_sampling_params["gibbs_max_iter"];
	BURNIN=gibbs_sampling_params["gibbs_burnin"];
	if(BURNIN>=MAX_GIBBS_ITER){
		cout<<"Error: BURNIN ("<<BURNIN<<") must be smaller than MAX_GIBBS_ITER ("<<MAX_GIBBS_ITER<<")\n";
		exit(1);
	}
	INTERVAL=gibbs_sampling_params["gibbs_interval"];
	ALPHA_A=gibbs_probability_params["gibbs_alpha_a"];
	ALPHA_B=gibbs_probability_params["gibbs_alpha_b"];
	GAMMA_A=gibbs_probability_params["gibbs_gamma_a"];
	GAMMA_B=gibbs_probability_params["gibbs_gamma_b"];
	H=gibbs_probability_params["gibbs_h"];
}

void initialize_param_maps(){
	gibbs_sampling_params["gibbs_max_iter"]=10000;
	gibbs_sampling_params["gibbs_burnin"]=5000;
	gibbs_sampling_params["gibbs_interval"]=50;
	gibbs_probability_params["gibbs_alpha_a"]=1.0;
	gibbs_probability_params["gibbs_alpha_b"]=1.0;
	gibbs_probability_params["gibbs_gamma_a"]=1.0;
	gibbs_probability_params["gibbs_gamma_b"]=1.0;
	gibbs_probability_params["gibbs_h"]=0.5;
}

void print_usage(){
	cout<<"Usage: structhdp uses the following input parameters in any order:\n";
	cout<<"-d <datafile> \n";
	cout<<"-n <number of individuals>\n";
	cout<<"-m <number of loci>\n";
	cout<<"-p <ploidy>\n";
	cout<<"-g <other-settings-file>\n";
	cout<<"--stirling_file <filename>\n";
	cout<<"--stirling_size <size of stirling matrix in file>\n";
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
        int dgiven=-1,ngiven=-1,ogiven=-1,pgiven=-1,rgiven=-1,ggiven=-1,mgiven=-1,sgiven=-1,fgiven=-1;
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
		{"gibbs_max_iter",required_argument,0,'M'},
		{"gibbs_burnin",required_argument,0,'B'},
		{"gibbs_interval",required_argument,0,'I'},
		{"gibbs_alpha_a",required_argument,0,'W'},
		{"gibbs_alpha_b",required_argument,0,'X'},
		{"gibbs_gamma_a",required_argument,0,'Y'},
		{"gibbs_gamma_b",required_argument,0,'Z'},
		{"gibbs_h",required_argument,0,'H'},
		{"stirling_file",required_argument,0,'F'},
		{"stirling_size",required_argument,0,'S'},
		{0,0,0,0}
		};
	while ((opt = getopt_long(len, arr,"d:o:n:m:p:hr:g:M:B:I:W:X:Y:Z:H:F:S:", 
                   long_options, &long_index )) != -1) {
        	switch (opt) {
			case 'F':
				STIRLINGFILE=string(optarg); fgiven=1;
				break;
			case 'S':
				STIRLING_SIZE=atoi(optarg); sgiven=1;
				break;
			case 'M':
			case 'B':
			case 'I':
				gibbs_sampling_params[long_options[long_index].name]=atoi(optarg);
				cout<<"Using custom value for "<<long_options[long_index].name<<"="<<optarg<<endl;
				break;
			case 'W':
			case 'X':
			case 'Y':
			case 'Z':
			case 'H':
				gibbs_probability_params[long_options[long_index].name]=atof(optarg);
				cout<<"Using custom value for "<<long_options[long_index].name<<"="<<optarg<<endl;
				break;
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
        verifyoption(sgiven,"Size of stirling matrix","--stirling_size");
        verifyoption(fgiven,"Stirling file","--stirling_file");
        verifyoption(dgiven,"Datafile","-d");
        verifyoption(ngiven,"Number of individuals","-n");
        verifyoption(mgiven,"Number of loci","-m");
        verifyoption(pgiven,"Ploidy","-p");
        //verifyoption(ggiven,"Other parameter file","-g");
        if(rgiven==-1)
                cout<<"No seed provided, using default value: "<<SEED<<endl;
        if(ogiven==-1)
                cout<<"No output prefix provided, using default value: "<<OUTPUT_PREFIX<<endl;

}

void readsettings(string filename){
	ifstream settingsfile;
	cout<<"Settings filename is "<<filename<<endl;
	settingsfile.open(filename.c_str());
	string settingname,settingvalues;
	float settingvaluef;
	int settingvaluei;
	while(!settingsfile.eof()){
		settingsfile>>settingname;
		if(settingname.compare("MAX_GIBBS_ITER")==0){
			settingsfile>>settingvaluei;
			MAX_GIBBS_ITER=settingvaluei;
		}
		else if(settingname.compare("BURNIN")==0){
			settingsfile>>settingvaluei;
			BURNIN=settingvaluei;
		}
		else if(settingname.compare("INTERVAL")==0){
			settingsfile>>settingvaluei;
			INTERVAL=settingvaluei;
		}
		else if(settingname.compare("ALPHA_A")==0){
			settingsfile>>settingvaluef;
			ALPHA_A=settingvaluef;
		}
		else if(settingname.compare("ALPHA_B")==0){
			settingsfile>>settingvaluef;
			ALPHA_B=settingvaluef;
		}
		else if(settingname.compare("GAMMA_A")==0){
			settingsfile>>settingvaluef;
			GAMMA_A=settingvaluef;
		}
		else if(settingname.compare("GAMMA_B")==0){
			settingsfile>>settingvaluef;
			GAMMA_B=settingvaluef;
		}
		else if(settingname.compare("H")==0){
			settingsfile>>settingvaluef;
			H=settingvaluef;
		}
		else if(settingname.compare("STIRLINGFILE")==0){
			settingsfile>>settingvalues;
			STIRLINGFILE=settingvalues;
		}
		else if(settingname.compare("STIRLING_SIZE")==0){
			settingsfile>>settingvaluei;
			STIRLING_SIZE=settingvaluei;
		}
		else{
			cout<<"Unknown setting "<<settingname;
			exit(1);
		}
			
	}
	settingsfile.close();
}

