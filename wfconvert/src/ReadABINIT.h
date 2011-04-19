#ifndef READ_ABINIT_H
#define READ_ABINIT_H

#include <Common/Blitz.h>

class ABINIT_PSPinfo
{
public:
  string title;     // 132 characters
  double znuclpsp;
  double zionpsp;
  int pspso;
  int pspdat;
  int pspcod;
  int pspxc;
  int lmax;
  int lloc;
  int lmn_size;
};

class ABINIT_header
{
public:
  int recHeaderLength;

  string codvsn;  // 6 characters
  int headform;
  int fform;
  int bantot;
  int date;
  int ixc;
  int intxc;
  int natom;
  int ngfft[3];
  int nkpt;
  int nspden;
  int nspinor;
  int nsppol;
  int nsym;
  int npsp;
  int ntypat;
  int occopt;
  int pertcase;
  int usepaw;
  double ecut, ecutdg, ecutsm, ecut_eff;
  Vec3 qptn;
  Vec3 rprimd[3];
  double stmbias;
  double tphysel;
  double tsmear;
  int usewvl;
  vector<int> istwfk;
  vector<int> nband;
  vector<int> npwarr;
  vector<int> so_psp;
  vector<int> symafm;
  vector<Mat3> symrel;       // size=nsym
  vector<int> typat;         // size=natom
  vector<int> nrhoijsel;     // size=nspden
  Array<int,2> rhoijselct;   // size=*,nspden
  vector<Vec3> kpt;          // size=nkpt
  vector<double> occ;        // size=bantot;
  vector<Vec3> tnons;        // size=nsym
  vector<double> znucltypat; // size=ntypat
  vector<double> wtk;        // size=nkpt
  string title;              // size=132
  vector<ABINIT_PSPinfo> PSPinfo;
  double residm;
  vector<Vec3> xred;         // size=natom
  double etotal;
  double fermie;
  void Read (FILE *fin, bool print=true);
  void SkipRecHeader(FILE *fin);
};


#endif
