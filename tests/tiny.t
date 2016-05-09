
These are some simple tests using the example dataset.

  $ STRUCTHDP=$TESTDIR/../structhdp
  $ DATA=$TESTDIR/../Data

Basic invocation of help:
  $ $STRUCTHDP -h
  Usage: structhdp uses the following input parameters in any order:
  -d <datafile> 
  -n <number of individuals>
  -m <number of loci>
  -p <ploidy>
  -g <other-settings-file>
  -o <output-prefix> (optional, default='structhdp')
  -r <random seed> (optional, default=1)

Normal command termination:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 -g $TESTDIR/../settings.txt -r 1 > /dev/null
  $ echo $?
  0

Expected command output:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 -g $TESTDIR/../settings.txt -r 1 | tail -n 20
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  1.0000 0.0000 0.0000 
  0.0000 1.0000 0.0000 
  0.0000 1.0000 0.0000 
  0.0000 0.0000 1.0000 
