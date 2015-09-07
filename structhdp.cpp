#include "utils.h"
#include "structhdp.h"
//#define DEBUG_2_2 1

gsl_rng * RANDOM_NUMBER;

model* initialise_model(int*** data,int nloci,int npops,int ploidy,int ninds){
	model* m=new model();
	//	INIT  no of pops, no of loci, ploidy
	m->npops=npops;
	m->nloci=nloci;
	m->ploidy=ploidy;
	cout<<"Model description is: npops="<<m->npops<<" nloci="<<m->nloci<<" ploidy="<<m->ploidy<<endl;
	cout<<"Number of inds is "<<ninds<<endl;
	//	INIT centroid arrays
	m->numcentroids=new int[nloci];
	m->uniqalleles=new map<int,int>[nloci];
	int allelecount,thislocus;
	for(int i=0;i<m->nloci;i++){
		m->uniqalleles[i] =map<int,int>(); // inited allele unique map for this locus
#ifdef DEBUG_0
		cout<<"Doing locus "<<i<<endl;
#endif
		allelecount=0;
		for(int j=0;j<ninds*m->ploidy;j++){
			thislocus=data[j/m->ploidy][i][j%m->ploidy];
			if(thislocus>=0){
				if(m->uniqalleles[i].find(thislocus)==m->uniqalleles[i].end()){
					m->uniqalleles[i][thislocus]=allelecount; //added new unique allele to list for locus i
					data[j/m->ploidy][i][j%m->ploidy]=allelecount++;
				}
				else{
					data[j/m->ploidy][i][j%m->ploidy]=(*(m->uniqalleles[i].find(thislocus))).second;
				}
				// destructively modifies data array, datacopy array left untouched.
			}
		}
		m->numcentroids[i]=m->uniqalleles[i].size();
	}
	cout<<"Printing numcentroids\n";
	for(int i=0;i<m->nloci;i++)
		cout<<m->numcentroids[i]<<" ";
	cout<<endl;
	m->alpha_a=ALPHA_A;//1.0;
	m->alpha_b=ALPHA_B;//1.0;
	m->alpha=(m->alpha_a*m->alpha_b);

	m->gamma_a=GAMMA_A;//1.0;
	m->gamma_b=GAMMA_B;//1.0;
	m->gamma=(m->gamma_a*m->gamma_b);

	m->h=H;//0.5;
	stringstream ss;
	cout<<"Stirlingfile is "<<STIRLINGFILE<<endl;
	ifstream ifile;
	double temp;
	ifile.open(STIRLINGFILE.c_str());
	for(int i=0;i<STIRLING_SIZE;i++)
		for(int j=0;j<STIRLING_SIZE;j++){
			ifile>>temp;
			m->stirlingmatrix[i][j]=temp;
		}
	cout<<"Read stirling matrix\n";
	ifile.close();
#ifdef DEBUG
	cout<<"Returning brand new model\n";
#endif
	return m;
}

void init_z_fixed_topics(satellite* info,model* m,int*** data,int kinit){
	int chosenind;
	vector< vector< vector<int> > > freqs;
	freqs.resize(kinit);
	for(int i=0;i<kinit;i++){
		freqs[i].resize(info->nloci);
		for(int j=0;j<info->nloci;j++){
			freqs[i][j].resize(info->numcentroids[j],0);
		}
	}
	int assigned[info->ninds];
	for(int i=0;i<info->ninds;i++)
		assigned[i]=-1;
	for(int k=0;k<kinit;k++){
		do{
			chosenind=runiform_int(info->ninds);
		}
		while(assigned[chosenind]>-1);
		assigned[chosenind]=k;
		info->njdotk[chosenind][k]=info->njdotdot[chosenind];
		for(int j=0;j<info->nloci;j++){
			for(int p=0;p<info->ploidy;p++){
				if(data[chosenind][j][p]>=0){
					freqs[k][j][data[chosenind][j][p]]++;
					info->zmatrix[info->ploidy*chosenind+p][j]=k;
					info->count_lwk[j][data[chosenind][j][p]][k]++;
				}
			}
		}
	}
	vector< vector<int> > ind_freq;
	ind_freq.resize(info->nloci);
	for(int j=0;j<info->nloci;j++){
		ind_freq[j].resize(info->numcentroids[j],0);
	}
	int chosenk;
	double simvals[kinit],maxval;
	for(int i=0;i<info->ninds;i++){
		if(assigned[i]>-1)
			continue;
		for(int j=0;j<info->nloci;j++){
			ind_freq[j].assign(info->numcentroids[j],0);
		}
		for(int j=0;j<info->nloci;j++){
			for(int p=0;p<info->ploidy;p++){
				if(data[i][j][p]>=0)
					ind_freq[j][data[i][j][p]]++;
			}
		}
		for(int k=0;k<kinit;k++){
			simvals[k]=0;
			for(int j=0;j<info->nloci;j++){
				simvals[k]+=similarity(&ind_freq[j][0],&freqs[k][j][0],info->numcentroids[j]);
			}
		}
		maxval=max(simvals,kinit,&chosenk);
		for(int j=0;j<info->nloci;j++){
			for(int p=0;p<info->ploidy;p++){
				if(data[i][j][p]>=0){
					freqs[chosenk][j][data[i][j][p]]++;
					info->zmatrix[info->ploidy*i+p][j]=chosenk;
					info->count_lwk[j][data[i][j][p]][chosenk]++;
				}
			}
		}
		info->njdotk[i][chosenk]+=info->njdotdot[i];
	}

}

satellite* init_from_model(model* m,int numinds,int*** data){
	//Populate everything 
	satellite* info=new satellite();
	info->nloci=m->nloci;
	info->ninds=numinds;
	info->ploidy=m->ploidy;
	info->npops=m->npops;
	info->allocated_pops=INIT_ALLOCATED_POPS;
	info->numcentroids=new int[info->nloci];
	//Initialize number of centroids
	for(int i=0;i<m->nloci;i++){
		info->numcentroids[i]=m->numcentroids[i];
	}
	info->njdotdot.resize(info->ninds,0);
	for(int i=0;i<info->ninds;i++)
		for(int j=0;j<info->nloci;j++)
			for(int p=0;p<info->ploidy;p++)
				if(data[i][j][p]>=0)
					info->njdotdot[i]++;

	info->njdotk.resize(info->ninds);
	info->mjk.resize(info->ninds);
	info->mjdot.resize(info->ninds,0);
	for(int i=0;i<info->ninds;i++){
		info->njdotk[i].resize(info->allocated_pops,0);
		info->mjk[i].resize(info->allocated_pops,0);
	}
	info->count_lwk.resize(info->nloci);
	for(int i=0;i<info->nloci;i++){
		info->count_lwk[i].resize(info->numcentroids[i]);
		for(int j=0;j<info->numcentroids[i];j++){
			info->count_lwk[i][j].resize(info->allocated_pops,0);
		}
	}
	info->zmatrix=new int*[info->ninds*info->ploidy];
	int allele;
	cout<<"Updating counts\n";
	for(int i=0;i<info->ninds*info->ploidy;i++){
		info->zmatrix[i]=new int[info->nloci];
	}
	info->npops=INIT_NUM_TOPICS;
	init_z_fixed_topics(info,m,data,INIT_NUM_TOPICS);
	info->mdotk.resize(info->allocated_pops,0);
	info->betas.resize(info->allocated_pops+1,0);
	info->mdotdot=0;
	for(int i=0;i<INIT_NUM_TOPICS;i++){
		info->mdotk[i]=5;
		info->mdotdot+=info->mdotk[i];
		info->betas[i]=(1.0)/(1+INIT_NUM_TOPICS);	
	}
	info->betas[INIT_NUM_TOPICS]=(1.0)/(1+INIT_NUM_TOPICS);	
	cout<<"Created satellite info\n";
	return info;
}

double logfcomputeold(int allele,int locus,int pop,satellite* info,model* m){
	double sum=0;
	for(int i=0;i<info->numcentroids[locus];i++){
		sum+=(info->count_lwk[locus][i][pop]+m->h);
	}
	double ansnum= (info->count_lwk[locus][allele][pop]+m->h);
	double ans=log(ansnum/sum);
	if(isnan(ans)){
		cout<<"Logfcomputeold diagnostics\n";
		cout<<"Allele "<<allele<<endl;
		cout<<"Locus "<<locus<<endl;
		cout<<"Pop no "<<pop<<endl;
		cout<<"Countval "<<info->count_lwk[locus][allele][pop]<<endl;
		cout<<"Numerator "<<ansnum<<endl;
		cout<<"Sum "<<sum<<endl;
		cout<<"Ans "<<ans<<endl;
		while(1){}
	}
	return ans;
}

double log_stirling_num(model* curr,int n,int m,bool flag){
	if(n>=STIRLING_SIZE || m>=STIRLING_SIZE){
		cout<<"Stirling called on a number larger than capacity: "<<n<<" "<<m<<"\n";
	}
	if(n==0 || m==0 || m> n)
		return -100.0;
	double ans=curr->stirlingmatrix[n-1][m-1];  //due to zero-indexing
	if(flag){
		cout<<"Stiring ("<<n<<" "<<m<<") = "<<ans<<endl;
	}
	return ans;
}


void samplehyperparams(satellite* info,model* m,int*** data){
	//Update alpha
	double  shape = m->alpha_a;
	double  scale = m->alpha_b;

	int n = info->mdotdot;
	double rate, sum_log_w, sum_s;

	for (int step = 0; step < 30; step++)
	{
		sum_log_w = 0.0;
		sum_s = 0.0;
		for (int d = 0; d < info->ninds; d++)
		{
			sum_log_w += log(rbeta(m->alpha + 1, info->njdotdot[d]));
			sum_s += (double)rbernoulli(info->njdotdot[d] / (info->njdotdot[d] + m->alpha));
		}
		rate = 1.0 / scale - sum_log_w;
		m->alpha = rgamma(shape + n - sum_s, 1.0 / rate);
		//cout<<"Step "<<step<<", alpha = "<<m->alpha<<endl;
	}

	if(m->alpha<1e-10)
		m->alpha=1e-10;
	//Update gamma
	unsigned int cc;
	int k=info->npops;
	shape = m->gamma_a;
	scale = m->gamma_b;

	double eta = rbeta(m->gamma + 1, n);
	double pi = shape + k - 1;
	rate = 1.0 / scale - log(eta);
	pi = pi / (pi + rate * n);

	cc = rbernoulli(pi);
	if (cc == 1)
		m->gamma = rgamma(shape + k, 1.0 / rate);
	else
		m->gamma = rgamma(shape + k - 1, 1.0 / rate);
}

void samplebeta(satellite* info,model* m,int*** data){
	double d_alpha[info->npops+1];
	for(int i=0;i<info->npops;i++)
		d_alpha[i]=info->mdotk[i];
	d_alpha[info->npops]=m->gamma;
	rdirichlet(info->npops+1,d_alpha,&(info->betas[0]));
}

void samplem(satellite* info,model* m,int*** data){
	int njdotkcurr,chosenm;
	double* parr;
	info->mdotdot=0;
	for(int k=0;k<info->npops;k++)
		info->mdotk[k]=0;
	for(int i=0;i<info->ninds;i++){
		info->mjdot[i]=0;
		for(int k=0;k<info->npops;k++){
			info->mjk[i][k]=0;
			njdotkcurr=info->njdotk[i][k];
			if(njdotkcurr==0){
				continue;
			}
			vector<double> pvec(njdotkcurr,-1000);
			for(int m1=0;m1<njdotkcurr;m1++){
				pvec[m1]=lgamma(m->alpha*info->betas[k])-lgamma(m->alpha*info->betas[k]+njdotkcurr)+log_stirling_num(m,njdotkcurr,m1+1,false)+(m1+1)*log(m->alpha)+(m1+1)*log(info->betas[k]);
			}		
			log_normalize(pvec,pvec.size());
			parr=new double[njdotkcurr];
			for(int m1=0;m1<njdotkcurr;m1++){
				parr[m1]=exp(pvec[m1]);
			}
			chosenm=rmultinomial(parr,njdotkcurr);
			delete parr;
			pvec.clear();
			info->mjk[i][k]=(chosenm+1);
			info->mdotk[k]+=(chosenm+1);
			info->mdotdot+=(chosenm+1);
			info->mjdot[i]+=(chosenm+1);
		}
	}
}

void samplez(satellite* info,model* m,int*** data){
	int chosenk,allele,oldk;
	double b,temp;
	bool addedpop=false;
	for(int i=0;i<info->ninds*info->ploidy;i++){
		for(int j=0;j<info->nloci;j++){
			allele=data[i/info->ploidy][j][i%info->ploidy];
			if(allele<0)
				continue;
			oldk=info->zmatrix[i][j];
			info->njdotk[i/info->ploidy][oldk]--;
			info->count_lwk[j][allele][oldk]--;
			if(info->count_lwk[j][allele][oldk]<0){
				cout<<"This should not happen\n";
				cout<<j<<" "<<allele<<" "<<oldk<<" "<<info->count_lwk[j][allele][oldk]<<endl;
				while(1){}
			}
			vector<double> pvec(info->npops+1,-1000);
			for(int k=0;k<info->npops;k++){
				temp=logfcomputeold(allele,j,k,info,m);
				pvec[k]=log(info->njdotk[i/info->ploidy][k]+m->alpha*info->betas[k])+temp;
				if(isnan(pvec[k])){
					cout<<"Pvec diagnostics\n";
					cout<<"info->njdotk[i/info->ploidy][k] "<<info->njdotk[i/info->ploidy][k]<<endl;
					cout<<"m->alpha "<<m->alpha<<endl;
					cout<<"info->betas[k] "<<info->betas[k]<<endl;
					cout<<"logfcomputeold(allele,j,k,info,m) "<<temp<<endl;
					while(1){}
				}
			}
			pvec[info->npops]=log(m->alpha*info->betas[info->npops])-log(info->numcentroids[j]);
			log_normalize(pvec,pvec.size());
			double* parr=new double[pvec.size()];
			for(int k=0;k<info->npops+1;k++)
				parr[k]=exp(pvec[k]);
			chosenk=rmultinomial(parr,info->npops+1);
			info->zmatrix[i][j]=chosenk;
			if(chosenk<info->npops){
				info->njdotk[i/info->ploidy][chosenk]++;
				info->count_lwk[j][allele][chosenk]++;
			}
			else{ // new topic
				addedpop=true;
				b=rbeta(1,m->gamma);
				info->betas[info->npops+1]=(1-b)*info->betas[info->npops];
				info->betas[info->npops]*=b;
				info->npops++;
				if(chosenk<info->allocated_pops){
					info->njdotk[i/info->ploidy][chosenk]++;
					info->count_lwk[j][allele][chosenk]++;
				}
			}
			delete parr;
			pvec.clear();
		}
	}
}

double computelikelihood(model* m,satellite* info,int*** data){
	double ll=0;
	int nldotk;
	ll+=lgamma(m->gamma)-lgamma(m->gamma+info->mdotdot); 
	ll+=info->npops*log(m->gamma);
	for(int k=0;k<info->npops;k++){
		ll+=lgamma(info->mdotk[k]);
		for(int l=0;l<info->nloci;l++){
			nldotk=0;
			for(int w=0;w<info->numcentroids[l];w++){
				nldotk+=info->count_lwk[l][w][k];
				ll+=(lgamma(m->h+info->count_lwk[l][w][k])-lgamma(m->h));
			}
			ll+=(lgamma(info->numcentroids[l]*m->h)-lgamma(info->numcentroids[l]*m->h+nldotk));
		}
	}
	for(int i=0;i<info->ninds;i++){
		ll+=(lgamma(m->alpha)-lgamma(m->alpha+info->njdotdot[i]));
		ll+=info->mjdot[i]*log(m->alpha);
		for(int k=0;k<info->npops;k++)
			if(info->njdotk[i][k]>0)
				ll+=log_stirling_num(m,info->njdotk[i][k],info->mjk[i][k],false);
	}
	double shape = m->gamma_a;
        double scale = m->gamma_b;
        ll += (shape-1)*log(m->gamma) - m->gamma/scale - shape*log(scale) - lgamma(shape);

        shape = m->alpha_a;
        scale = m->alpha_b;
        ll += (shape-1)*log(m->alpha) - m->alpha/scale - shape*log(scale) - lgamma(shape);

	return ll;
}


datamatrix* readdata(string filename){
	ifstream ifile;
	int temp;
	cout<<"Filename is "<<filename<<endl;
	ifile.open(filename.c_str());
	datamatrix* d=new datamatrix();
	d->ninds=NUMINDS;
	d->nloci=NUMLOCI;
	d->ploidy=PLOIDY;
	d->data=new int**[d->ninds];
	d->datacopy=new int**[d->ninds];
	for(int i=0;i<d->ninds;i++){
		d->data[i]=new int*[d->nloci];
		d->datacopy[i]=new int*[d->nloci];
		for(int j=0;j<d->nloci;j++){
			d->data[i][j]=new int[d->ploidy];
			d->datacopy[i][j]=new int[d->ploidy];
		}
	}
	for(int i=0;i<d->ninds;i++){
		for(int p=0;p<d->ploidy;p++){
			for(int j=0;j<d->nloci;j++){
				ifile>>temp;
				d->data[i][j][p]=temp;
				d->datacopy[i][j][p]=temp;
			}
		}
	}
	ifile.close();
	return d;
}

void clearmatrix(datamatrix* d){
	for(int i=0;i<d->ninds;i++){
		for(int p=0;p<d->ploidy;p++){
			delete[] d->data[i][p];
		}
		delete[] d->data[i];
	}
	delete[] d->data;
	delete d;
}

void compact_state(satellite* info,model* m){
	int num_topics_old=info->npops;
	int* k_to_new_k = new int[num_topics_old];
	int k, new_k;
	for (k = 0, new_k = 0; k < num_topics_old; k++)
	{
		if (info->mdotk[k] > 0)
		{
			k_to_new_k[k] = new_k;
			swap_vec_element(info->mdotk,  new_k, k);
			swap_vec_element(info->betas,  new_k, k);
			for(int i=0;i<info->ninds;i++){
				swap_vec_element(info->njdotk[i],new_k,k);
				swap_vec_element(info->mjk[i],new_k,k);
			}
			for(int j=0;j<info->nloci;j++)
				for(int w=0;w<info->numcentroids[j];w++)
					swap_vec_element(info->count_lwk[j][w],new_k,k);	
			new_k ++;
		}
	}
	for(int i=0;i<info->ninds*info->ploidy;i++)
		for(int j=0;j<info->nloci;j++){
			info->zmatrix[i][j]=k_to_new_k[info->zmatrix[i][j]];
		}
	info->npops = new_k;
	delete k_to_new_k;
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
		else{
			cout<<"Unknown setting "<<settingname;
			exit(1);
		}
			
	}
	settingsfile.close();
}

void learnmodel(string filename){
#ifdef DEBUG
	cout<<"Starting file read\n";
#endif
	datamatrix* d=readdata(filename);
#ifdef DEBUG
	cout<<"Read file\n";
#endif

	model* currmodel=initialise_model(d->data,d->nloci,INIT_NUM_TOPICS,d->ploidy,d->ninds);
#ifdef DEBUG
	cout<<"Inited model\n";
#endif
	double currlkhd,oldlkhd;
	satellite* info=init_from_model(currmodel,d->ninds,d->data);
	cout.precision(4);
	cout<<fixed;

	ofstream stateout;
	stateout.open((OUTDIR+"/state.log").c_str());
	stateout<<"iter num.topics likelihood alpha gamma\n";

	double ksum=0.0;
	double knum=0;
	ofstream kout;
	kout.open((OUTDIR+"/klog.txt").c_str());
	for(int iter=1;iter<=MAX_GIBBS_ITER;iter++){
		samplez(info,currmodel,d->data);
		samplem(info,currmodel,d->data);
		compact_state(info,currmodel);
		samplebeta(info,currmodel,d->data);
		samplehyperparams(info,currmodel,d->data);
		currlkhd=computelikelihood(currmodel,info,d->data);
		cout<<"iter= "<<setw(5)<<iter;
		cout<<", K= "<<setw(3)<<info->npops;
		cout<<" ,#Pops= "<<setw(5)<<info->mdotdot;
		cout<<", alpha= "<<currmodel->alpha<<", gamma= "<<currmodel->gamma;
		cout<<", ll= "<<currlkhd;
		cout<<endl;
		stateout<<iter<<" "<<info->npops<<" "<<currlkhd<<" "<<currmodel->alpha<<" "<<currmodel->gamma<<"\n";
		if(iter%INTERVAL==0 && iter>BURNIN){
			kout<<info->npops<<endl;
		}
	}
	kout.close();

	samplem(info,currmodel,d->data);
	compact_state(info,currmodel);
	cout<<"Printing final ancestry proportions\n";

	ofstream propout;
	propout.open((OUTDIR+"/proportions.txt").c_str());

	double propsum=0;
	for(int i=0;i<info->ninds;i++){
		propsum=0;
		for(int j=0;j<info->npops;j++){
			propsum+=info->njdotk[i][j];
		}
		for(int j=0;j<info->npops;j++){
			cout<<setw(4)<<info->njdotk[i][j]/propsum<<" ";
			propout<<info->njdotk[i][j]/propsum<<" ";
		}
		cout<<endl;
		propout<<endl;
	}
	stateout.close();
	propout.close();
	clearmatrix(d);
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
        int count=1,len1;
        string temp;
        int dgiven=-1,ngiven=-1,ogiven=-1,pgiven=-1,rgiven=-1,ggiven=-1,kgiven=-1,mgiven=-1;
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
                                cout<<"Usage: structhdp requires the following parameters in any order:\n";
                                cout<<"-d <datafile>\n";
                                cout<<"-o <output-directory>\n";
                                cout<<"-n <number of individuals>\n";
                                cout<<"-m <number of loci>\n";
                                cout<<"-p <ploidy>\n";
                                cout<<"-g <other-settings-file>\n";
                                cout<<"You can provide a random seed (optional) with \"-r <integer seed>\" \n";
                                exit(0);
                        case 'd':
                                // Data file
                                DATAFILE=string(arr[count+1]); dgiven=1;
                                break;
                        case 'o':
                                //Output directory
                                OUTDIR=string(arr[count+1]); ogiven=1;
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
        verifyoption(ogiven,"Output directory","-o");
        if(rgiven==-1)
                cout<<"No seed provided, using default value: "<<SEED<<endl;
}



int main(int argc,char** argv){
	time_t t;
	time(&t);
	SEED=(long)t;
	readopts(argc,argv);
	RANDOM_NUMBER = gsl_rng_alloc(gsl_rng_taus);
    	gsl_rng_set(RANDOM_NUMBER, (long) SEED); // init the seed
	readsettings(OTHERPARAMS);
	learnmodel(DATAFILE);
	return 0;
}
