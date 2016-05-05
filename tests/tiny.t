
These are some simple tests using the example dataset.

  $ STRUCTHDP=$TESTDIR/../structhdp
  $ DATA=$TESTDIR/../Data

Basic invocation of help:
  $ $STRUCTHDP -h
  Usage: structhdp requires the following parameters in any order:
  -d <datafile>
  -o <output-directory>
  -n <number of individuals>
  -m <number of loci>
  -p <ploidy>
  -g <other-settings-file>
  You can provide a random seed (optional) with "-r <integer seed>" 

Normal command termination:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 -g $TESTDIR/../settings.txt -r 1 > /dev/null
  $ echo $?
  0

Expected command output:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 -g $TESTDIR/../settings.txt -r 1 | tail -n 20
  0.9231 0.0769 
  1.0000 0.0000 
  0.8462 0.1538 
  0.9286 0.0714 
  0.9231 0.0769 
  0.9286 0.0714 
  0.8571 0.1429 
  0.8571 0.1429 
  0.9231 0.0769 
  0.9286 0.0714 
  1.0000 0.0000 
  0.6923 0.3077 
  1.0000 0.0000 
  1.0000 0.0000 
  0.9286 0.0714 
  0.9286 0.0714 
  0.9286 0.0714 
  0.9286 0.0714 
  1.0000 0.0000 
  0.7692 0.2308 
