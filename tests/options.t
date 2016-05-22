
These are some simple tests using the example dataset.

  $ STRUCTHDP=$TESTDIR/../structhdp
  $ DATA=$TESTDIR/../Data

Expected command output:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 -r 1 --stirling_file $TESTDIR/../logstirling_500.txt --stirling_size 500 --gibbs_max_iter 1000 --gibbs_burnin 500 --gibbs_alpha_a 0.95 | head -n 5
  Using custom value for gibbs_max_iter=1000
  Using custom value for gibbs_burnin=500
  Using custom value for gibbs_alpha_a=0.95
  Starting file read
  Filename is /Users/suyashshringarpure/Work/structhdp/tests/../Data/tiny.stru
