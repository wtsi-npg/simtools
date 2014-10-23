
// temporary program (for now) to test intensity normalization


#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "commands.h"
#include "Manifest.h"
#include "Egt.h"
#include "Fcr.h"

using namespace std;

int main(int argc, char *argv[])
{
   // hack to test normalization method
  Commander *commander = new Commander();
  double xraw = 10399; // 1KG_1_100177980 285293_A01_SC_SEPI5488306
  double yraw = 672;
  double xnorm;
  double ynorm;
  string gtcfile = "/nfs/users/nfs_i/ib5/data/fcr_devel/9298751015_R01C01.gtc";
  string manfile0 = "/nfs/users/nfs_i/ib5/data/fcr_devel/HumanCoreExome-12v1-0_A.normalized.bpm.csv";
  string snpName = "1KG_1_100177980";
  Gtc *gtc = new Gtc();
  Manifest *manifest = new Manifest();
  gtc->open(gtcfile, Gtc::XFORM | Gtc::INTENSITY);
  manifest->open(manfile0);
  unsigned int norm_id = 1e9;
  for (unsigned int i=0;i<manifest->snps.size();i++) {
    if (manifest->snps[i].name == snpName) {
      norm_id = manifest->normIdMap[manifest->snps[i].normId];
      break;
    }
  }
  if (norm_id == 1e9) {
    throw;
  }
  commander->normalizeIntensity(xraw, yraw, xnorm, ynorm, norm_id, gtc);
  cout << "X = " << xnorm << ", Y = " << ynorm << endl;

  double theta;
  double r;
  Fcr *fcr0 = new Fcr();
  fcr0->cartesianToPolar(xnorm, ynorm, theta, r);
  cout << "theta = " << theta << ", R = " << r << endl;    

  return 0;

}
