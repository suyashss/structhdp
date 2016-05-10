#include "structhdp.h"
#include "options.h"

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
        int count=1,len1;
        string temp;
        int dgiven=-1,ngiven=-1,ogiven=-1,pgiven=-1,rgiven=-1,ggiven=-1,kgiven=-1,mgiven=-1;
	SEED=1;
	OUTPUT_PREFIX="structhdp";
        while(count!=len){
                temp=string(arr[count]);
                len1=temp.size();
                if(len1!=2){
                        cout<<"Invalid switch "<<arr[count]<<endl;
                        exit(1);
                }
                if(arr[count][1]!='h')
                        valuegiven(count,len,arr);
                switch (arr[count][1]) {
                        case 'h':
                                //help option
                                cout<<"Usage: structhdp uses the following input parameters in any order:\n";
                                cout<<"-d <datafile> \n";
                                cout<<"-n <number of individuals>\n";
                                cout<<"-m <number of loci>\n";
                                cout<<"-p <ploidy>\n";
                                cout<<"-g <other-settings-file>\n";
                                cout<<"-o <output-prefix> (optional, default='structhdp')\n";
				cout<<"-r <random seed> (optional, default=1)\n";
                                exit(0);
                        case 'd':
                                // Data file
                                DATAFILE=string(arr[count+1]); dgiven=1;
                                break;
                        case 'o':
                                //Output prefix
                                OUTPUT_PREFIX=string(arr[count+1]); ogiven=1;
                                break;
                        case 'n':
                                // No. of individuals
                                NUMINDS=atoi(arr[count+1]); ngiven=1;
                                break;
                        case 'm':
                                // No. of loci
                                NUMLOCI=atoi(arr[count+1]); mgiven=1;
                                break;
                        case 'r':
                                // Random seed
                                SEED=atoi(arr[count+1]); rgiven=1;
                                break;
                        case 'p':
                                // Ploidy
                                PLOIDY=atoi(arr[count+1]); pgiven=1;
                                break;
                                //break;
                        case 'g':
                                // Global params
                                OTHERPARAMS=string(arr[count+1]); ggiven=1;
                                break;
                        default:
                                cout<<"Unknown switch "<<arr[count]<<endl;
                                exit(1);
                }
                count+=2;
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


