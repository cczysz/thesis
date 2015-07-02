def snapresults(test_rsids,pop,rsq):
  """
  Use curl to submit a POST form to SNAP in order to retrieve pairwise linkage data for input rsids
  """

  # rsidString = "%0D%0A".join([item[0] for item in test_rsids])
  rsidString = "%0D%0A".join(test_rsids)

  hapMapPanel = pop

  if pop is "JPT" or pop is "CHB":
    hapMapPanel = "CHBJPT" #SNAP combines these two groups

  searchString = """
    SearchPairwise=
    &snpList=%s
    &hapMapRelease=onekgpilot
    &hapMapPanel=%s
    &RSquaredLimit=0.8
    &distanceLimit=500000
    &downloadType=File
    &includeQuerySnp=on
    &arrayFilter=query
    &columnList[]=DP
    &columnList[]=GP
    &columnList[]=AM
    &submit=search
    """ % (rsidString, hapMapPanel)#, rsq

  subprocess.call(["curl","-s","-d",searchString,"-o","SNAPResults_"+rsq+".txt","http://www.broadinstitute.org/mpg/snap/ldsearch.php"])

  time.sleep(0.5) #wait to prevent server overload

import subprocess
import time

if __name__ == '__main__':
  snapresults()
