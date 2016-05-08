
These are some simple tests using the example dataset.

  $ STRUCTHDP=$TESTDIR/../structhdp
  $ DATA=$TESTDIR/../Data

Sterling size check:
  $ $STRUCTHDP -d $DATA/tiny.stru -o output -n 20 -m 7 -p 2 -g $DATA/bad_settings.txt -r 1 
  Settings filename is /Users/suyashshringarpure/Work/structhdp/tests/../Data/bad_settings.txt
  Starting file read
  Filename is /Users/suyashshringarpure/Work/structhdp/tests/../Data/tiny.stru
  Read file
  Error: STIRLING_SIZE (10) is less than Nloci*Ploidy (7*2=14)- the program will fail
  Please use the write_stirling.r to generate a larger STIRLINGFILE. See README for more detailed instructions.
  [1]
