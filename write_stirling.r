
logsum<-function(log_a,log_b){
if (log_a < log_b){
        v = log_b+log1p( exp(log_a-log_b));
	}
    else{
        v = log_a+log1p(exp(log_b-log_a));
	}
    v
}

logstirling<-function(n,m){
	if(m>n || (m==0 && n > 0)){
		value=-1000;
	}
	else if(m==n){
		value=0;
	}
	else{
		value=logsum(logstirling(n-1,m-1),log(n-1)+logstirling(n-1,m));
	}		
	value;
}

logstirlingmatrix<-function(n){
#Compute logstirling numbers: stirling(n,1)... stirling(n,n)
value=matrix(0,nrow=n,ncol=n);
value[1,1]=0;
if(n>1){
for(i in 2:n){
	value[1,i]=-1000
	value[i,1]=lfactorial(i-1)
}
for(i in 2:n){
	cat("I is",i,"\n")
	for(j in 2:n){
		if(j>i){
			value[i,j]=-1000;
		}
		else{
			value[i,j]=logsum(value[i-1,j-1],log(i-1)+value[i-1,j]);
		}		
	}
}
}
value
}

writelogstirling<-function(n,filename){
	value=logstirlingmatrix(n)
	write(t(value),file=filename,ncolumns=n)
}

args<-commandArgs(T)

if(length(args)!=1){
	stop("Usage: Rscript write_stirling.r <n>\n")
}

outfile=paste("logstirling_",args[[1]],".txt",sep="")
cat("Output file is",outfile,"\n")
writelogstirling(as.numeric(args[[1]]),outfile)
