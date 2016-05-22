
These are some simple tests using the example dataset.

  $ STRUCTHDP=$TESTDIR/../structhdp
  $ DATA=$TESTDIR/../Data

Stirling size check:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 --stirling_file /home/suyash/structHDP/logstirling_500.txt --stirling_size 10 -r 1 
  Starting file read
  Filename is /Users/suyashshringarpure/Work/structhdp/tests/../Data/tiny.stru
  Read file
  Error: STIRLING_SIZE (10) is less than Nloci*Ploidy (7*2=14)- the program will fail
  Please use the write_stirling.r to generate a larger STIRLINGFILE. See README for more detailed instructions.
  [1]

Stirling small matrix:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 --stirling_file /home/suyash/structHDP/logstirling_20.txt --stirling_size 20 -r 1  > /dev/null
  $ echo $?
  0
