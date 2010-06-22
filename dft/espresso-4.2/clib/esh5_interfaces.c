/*
 * Copyright (C) 2004 PWSCF group 
 * Copyright (C) 2007 QMCPACK developers
 *
 * @author Jeongnim Kim http://www.mcc.uiuc.edu/qmcpack/
 * @brief Implements generic hdf5 interfaces for plane wave codes and qmcpack
 *
 * - esh5_open_file: open hdf5 file
 * - esh5_close_file : close hdf5 file
 * - esh5_open_eigg : open eigenstates
 * - esh5_close_eigg : close eigenstates
 * - esh5_open_eigr : open eigenstates_nx_ny_nz
 * - esh5_close_eigr : close eigenstates_nx_ny_nz
 * 
 */
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <assert.h>
#include "c_defs.h"
#include "hdf5.h"
#include "hdf5_hl.h"

/* file handler */
static hid_t h_file=-1;
/* handler for electrons or atoms*/
static hid_t h_ptcls=-1;
/* kpoint handler */
static hid_t h_kpoint=-1;
/* spin handler */
static hid_t h_spin=-1;
/* density handler */
static hid_t h_density=-1;
/* number of fft grid */
static int num_grid[3];
/* number of real-space grids */
static int h_ngridtot=0;
/* check for gamma */
static int is_gamma=0;
/* number of atoms */
static int num_atoms=0;
/* number of atom species */
static int num_species=0;
/* number of electrons */
static int num_els[2];
/* number of spin channels */
static int num_spins=1;
/* number of bands */
static int num_bands=0;
/* number of gvectors */
static int num_gvectors=0;
/* current k-point */
static int kpoint_now=-1;
/* is complex orbital */
static int psi_r_is_complex=1;
/* append data */
static int append_h5=0;
static int iteration=0;
static H5E_auto_t err_func;
static void *client_data=0;


/** create a file and write version & application
 * @param fname name of the output file
 * @param length size of the file name
 *
 * h_file is initialized.
 */
void F77_FUNC_(esh5_open_file,ESH5_OPEN_FILE)(const char* fname, const int* length, int* old)
{
  H5Eget_auto (&err_func, &client_data);
  H5Eset_auto (NULL, NULL);

  append_h5=*old;

  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ; 

  if(h_file>=0) H5Fclose(h_file); 
  h_file = H5Fopen(hfname,H5F_ACC_RDWR,H5P_DEFAULT);
  //if((append_h5)||(iteration))
  //{
  //  printf("esh5 open existing %s\n",hfname);
  //  h_file = H5Fopen(hfname,H5F_ACC_RDWR,H5P_DEFAULT);
  //}
  if(h_file<0)
  {
    printf("esh5 create %s\n",hfname);
    h_file = H5Fcreate(hfname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    /* impelements version 1.00 hdf5 format */
    int version[]={2,0,0};
    hsize_t dim=3;
    herr_t ret=H5LTmake_dataset(h_file,"version",1,&dim,H5T_NATIVE_INT,version);
    hsize_t ns=1;
    {
      hid_t strtype = H5Tcopy (H5T_C_S1);
      ret = H5Tset_size (strtype, 7); /* create string of length 5 */
      ret=H5LTmake_dataset(h_file,"format",1,&ns,strtype,"ES-HDF");
    }

    hid_t h_app = H5Gcreate(h_file,"application",0);
    {
      hid_t strtype = H5Tcopy (H5T_C_S1);
      ret = H5Tset_size (strtype, 8); /* create string of length 5 */
      ret=H5LTmake_dataset(h_app,"code",1,&ns,strtype,"espresso");
    }
    version[0]=4;
    version[2]=4;
    ret=H5LTmake_dataset(h_app,"version",1,&dim,H5T_NATIVE_INT,version);
    H5Gclose(h_app);
  }
  free(hfname);
  iteration = iteration+1;
}

void F77_FUNC_(esh5_close_file,ESH5_CLOSE_FILE)()
{
  if(h_file>=0) H5Fclose(h_file);
  h_file=-1;
  H5Eset_auto (err_func, client_data);
}

/** create electrons and create sub groups
 * @param nels_up 
 * @param nels_down
 * @param nspins number_of_spins
 * @param nkpts number_of_kpoints
 * @param nband number of electron states
 * @param ngr 3D mesh
 */
void F77_FUNC_(esh5_open_electrons,ESH5_OPEN_ELECTRONS)
  ( const int* nels_up, const int* nels_down , const int* nspins
    , const int* nkpts , const int *nband , const int* ngr
  )
{
  //save the values
  num_els[0]=*nels_up;
  num_els[1]=*nels_down;
  num_spins=*nspins;
  num_grid[0]=ngr[0];
  num_grid[1]=ngr[1];
  num_grid[2]=ngr[2];
  num_bands=*nband;

  h_ptcls = H5Gopen(h_file,"electrons");
  if(h_ptcls<0)
  {
    printf("Creating electrons\n");
    h_ptcls = H5Gcreate(h_file,"electrons",0);

    //save the number of up and down electrons
    const hsize_t dim1=1;
    const hsize_t dim2=2;
    const hsize_t dim3=3;
    herr_t ret=H5LTmake_dataset(h_ptcls,"number_of_electrons",1,&dim2,H5T_NATIVE_INT,num_els);
    ret=H5LTmake_dataset(h_ptcls,"number_of_spins",1,&dim1,H5T_NATIVE_INT,nspins);
    ret=H5LTmake_dataset(h_ptcls,"number_of_kpoints",1,&dim1,H5T_NATIVE_INT,nkpts);
    ret=H5LTmake_dataset(h_ptcls,"psi_r_mesh",1,&dim3,H5T_NATIVE_INT,ngr);
    //ret=H5LTmake_dataset(h_ptcls,"psi_r_is_complex",1,&dim1,H5T_NATIVE_INT,is_complex);

    //create kpoint/spin/state groups
    for(int ik=0; ik<*nkpts; ++ik)
    {
      char twistname[16];
      sprintf(twistname,"kpoint_%i",ik);
      hid_t h1 = H5Gcreate(h_ptcls,twistname,0);
      for(int ispin=0; ispin<num_spins; ispin++)
      {
        char spinname[16];
        sprintf(spinname,"spin_%i",ispin);
        hid_t h2 = H5Gcreate(h1,spinname,0);
        ret=H5LTmake_dataset(h2,"number_of_states",1,&dim1,H5T_NATIVE_INT,nband);
        for(int ib=0; ib<*nband; ib++) 
        {
          char bandname[16];
          sprintf(bandname,"state_%i",ib);
          hid_t h3 = H5Gcreate(h2,bandname,0);
          H5Gclose(h3);
        }
        H5Gclose(h2);
      }
      H5Gclose(h1);
    }
  }
}


void F77_FUNC_(esh5_close_electrons,ESH5_CLOSE_ELECTRONS) ()
{
  //if(append_h5)
  //{
    //write if psi_r is complex
  hsize_t dim1=1;
  printf("psi_r_is_complex");
  herr_t ret=H5LTmake_dataset(h_ptcls,"psi_r_is_complex",1,&dim1,H5T_NATIVE_INT,&psi_r_is_complex);
  //}

  H5Gclose(h_ptcls);
  h_ptcls=-1;
}

/** open kpoint 
 * @param ik the kpoint index
 */
void F77_FUNC_(esh5_open_kpoint,ESH5_OPEN_KPOINT)(const int* ik)
{
  kpoint_now=(*ik)-1;
  char kname[32];
  sprintf(kname,"kpoint_%i",kpoint_now);
  h_kpoint = H5Gopen(h_ptcls,kname);
  if (h_kpoint < 0) {
    fprintf (stderr, "Creating %s\n", kname);
    h_kpoint = H5Gcreate(h_ptcls,kname,0);
  }
 // assert (h_kpoint >= 0);
}
///* close kpoint */
void F77_FUNC_(esh5_close_kpoint,ESH5_CLOSE_KPOINT)()
{
  H5Gclose(h_kpoint);
}


/* write kpoint data */
void F77_FUNC_(esh5_write_kpoint_data,ESH5_WRITE_KPOINT_DATA)
(const double* xk, const double* wgt, const int* ngk_g)
// (const double* xk, const double* wgt, const int* ngk_g, const hsize_t* gints)
{
  hsize_t dim3=3;
  hsize_t dim1=1;
  hsize_t dim_g[2];
  dim_g[0] = *ngk_g;
  dim_g[1] = 3;
   
  if (iteration<2)
  {
    herr_t ret=H5LTmake_dataset(h_kpoint,"reduced_k",1,&dim3,H5T_NATIVE_DOUBLE,xk);
    ret=H5LTmake_dataset(h_kpoint,"weight",1,&dim1,H5T_NATIVE_DOUBLE,wgt);
    ret=H5LTmake_dataset(h_kpoint,"number_of_gvectors",1,&dim1,H5T_NATIVE_INT,ngk_g);
//     ret=H5LTmake_dataset(h_kpoint,"gvectors",2,dim_g,H5T_NATIVE_INT, gints);
  }
}

/** open spin
 * @param ispin the sin index
 */
void F77_FUNC_(esh5_open_spin,ESH5_OPEN_SPIN)(const int* ispin)
{
  char sname[32];
  sprintf(sname,"spin_%i",(*ispin)-1);
  h_spin=H5Gopen(h_kpoint,sname);
  if (h_spin < 0) {
    fprintf (stderr, "Creating %s\n", sname);
    h_spin=H5Gcreate(h_kpoint,sname,0);
  }
  assert (h_spin >= 0);
}

/* close kpoint */
void F77_FUNC_(esh5_close_spin,ESH5_CLOSE_SPIN)()
{
  H5Gclose(h_spin);
}


/* write eigen values 
 * @param ispin spin index
 * @param eigval eigen values 
 * @param nband number of bans
 */
void F77_FUNC_(esh5_write_eigvalues,ESH5_WRITE_EIGVALUES)(const double* eigval)
{
  hsize_t dim3=(hsize_t)num_bands;
  herr_t ret=H5LTmake_dataset(h_spin,"eigenvalues",1,&dim3,H5T_NATIVE_DOUBLE,eigval);
  H5Fflush(h_spin,H5F_SCOPE_GLOBAL);
  assert (ret >= 0);
}



/* write eigen value and eigen vector for (ibnd, ispin) */
void F77_FUNC_(esh5_write_psi_g,ESH5_WRITE_PSI_G)(const int* ibnd
    , const double* eigv, const int* ngtot
    )
{
  char aname[64];
  sprintf(aname,"state_%i/psi_g",(*ibnd)-1);
  hsize_t dims[2];
  dims[0] = (hsize_t)*ngtot;
  dims[1] = 2;
  //  fprintf(stderr, "aname = %s  ", aname);
  //  fprintf (stderr, "  ngtot = %d\n", *ngtot);
  herr_t ret=H5LTmake_dataset(h_spin,aname,2,dims,H5T_NATIVE_DOUBLE,eigv);
  assert (ret >= 0);

  if(h_ptcls>-1)
  {
    hid_t pid=H5Dopen(h_ptcls,"psi_r_is_complex");
    if(pid<0)
    {
      const hsize_t dim1=1;
      ret=H5LTmake_dataset(h_ptcls,"psi_r_is_complex",1,&dim1,H5T_NATIVE_INT,&psi_r_is_complex);
    }
    else
      H5Dclose(pid);
  }
}

/* write eigen value and eigen vector for (ibnd, ispin) */
void F77_FUNC_(esh5_write_psi_r,ESH5_WRITE_PSI_R)(const int* ibnd
    , const double* eigr, const int* use_complex
    )
{
  //need to flag this if they are not the same
  psi_r_is_complex=*use_complex;
  char aname[64];
  sprintf(aname,"state_%i/psi_r",(*ibnd)-1);
  hsize_t dims_out=(hsize_t)(psi_r_is_complex)?4:3;
  hsize_t dims[4];
  dims[0] = num_grid[0];
  dims[1] = num_grid[1];
  dims[2] = num_grid[2];
  dims[3] = 2;
  herr_t ret=H5LTmake_dataset(h_spin,aname,dims_out,dims,H5T_NATIVE_DOUBLE,eigr);
  assert (ret >= 0);
}


/** open density group and write its grid properties
 * @param g_qmc G in reduced coordinates
 * @param igtog index map
 * @param ngtot number of g vectors
 * @param nr1s grid of the first direction
 * @param nr2s grid of the second direction
 * @param nr3s grid of the third direction
 */
void F77_FUNC_(esh5_open_density,ESH5_OPEN_DENSITY)(const int* gint
    , const int* ngm, int *nr1s, int *nr2s, int *nr3s)
{
  num_grid[0]=*nr1s;
  num_grid[1]=*nr2s;
  num_grid[2]=*nr3s;
  num_gvectors=*ngm;

  h_density = H5Gcreate(h_ptcls,"density",0);
  hsize_t dim3=3;
  herr_t ret=H5LTmake_dataset(h_density,"mesh",1,&dim3,H5T_NATIVE_INT,num_grid);

  // {
  //   int *ig=(int*)malloc(3*num_gvectors*sizeof(int));
  //   for(int i=0,i3=0; i<num_gvectors; ++i)
  //   {
  //     int cur=3*(igtog[i]-1);
  //     ig[i3++]=(int)g[cur++];
  //     ig[i3++]=(int)g[cur++];
  //     ig[i3++]=(int)g[cur++];
  //   }

  //   hsize_t gdims[2];
  //   gdims[0] = (hsize_t)num_gvectors; 
  //   gdims[1] = (hsize_t)3;
  //   ret=H5LTmake_dataset(h_density,"gvectors",2,gdims,H5T_NATIVE_INT,ig);
  //   assert (ret >= 0);
  //   free(ig);
  // }
  hsize_t gdims[2];
  gdims[0] = (hsize_t)num_gvectors; 
  gdims[1] = (hsize_t)3;
  ret=H5LTmake_dataset(h_density,"gvectors",2,gdims,H5T_NATIVE_INT,gint);

  hsize_t dim1=1;
  ret=H5LTmake_dataset(h_density,"number_of_gvectors",1,
		       &dim1,H5T_NATIVE_INT,ngm);
}

/** open density group and write its grid properties
 * @param nr1s grid of the first direction
 * @param nr2s grid of the second direction
 * @param nr3s grid of the third direction
 */
void F77_FUNC_(esh5_open_density_r,ESH5_OPEN_DENSITY_R)(int *nr1s, int *nr2s, int *nr3s
    )
{
  printf("ARE YOU GONE MAD \n");
  num_grid[0]=*nr1s;
  num_grid[1]=*nr2s;
  num_grid[2]=*nr3s;

  h_density = H5Gcreate(h_ptcls,"density",0);
  hsize_t dim3=3;
  herr_t ret=H5LTmake_dataset(h_density,"mesh",1,&dim3,H5T_NATIVE_INT,num_grid);
}

void F77_FUNC_(esh5_close_density,ESH5_CLOSE_DENSITY)()
{
  H5Gclose(h_density);
}

/* write eigen value and eigen vector for (ibnd, ispin) */
void F77_FUNC_(esh5_write_density_r,ESH5_WRITE_DENSITY_R)(const int* ispin,const double* rho)
{
  char aname[8];
  sprintf(aname,"spin_%i",(*ispin)-1);
  /*hid_t h2 = H5Gcreate(h_density,aname,0);*/
  hid_t h2 = H5Gopen(h_density,aname);
  /* write eigenvector */
  hsize_t dims[3];
  for(int i=0; i<3; ++i) dims[i] = num_grid[i];
  herr_t ret=H5LTmake_dataset(h2,"density_r",3,dims,H5T_NATIVE_DOUBLE,rho);
  H5Gclose(h2);
}

void F77_FUNC_(esh5_write_density_g,ESH5_WRITE_DENSITY_G)
     (const int* ispin , const double* rhog)
{
  char aname[8];
  sprintf(aname,"spin_%i",(*ispin)-1);
  /*hid_t h2 = H5Gopen(h_density,aname);*/
  hid_t h2 = H5Gcreate(h_density,aname,0);
  hsize_t dims_g[2];
  dims_g[0]=num_gvectors;
  dims_g[1]=2;
  herr_t ret=H5LTmake_dataset(h2,"density_g",2,dims_g,H5T_NATIVE_DOUBLE,rhog);
  H5Gclose(h2);
}

/** write basisset: number of plane waves, plane wave coefficients
 */
void F77_FUNC_(esh5_write_gvectors,ESH5_WRITE_GVECTORS)
(const int* itmp, const int* igwk, const int* ngk_g)
{
  if (iteration<2)
  {
  int ngtot=*ngk_g;
  //int ng=*ngtot;
  int *igmapped=(int*)malloc(3*ngtot*sizeof(int));

  for(int ig=0,i3=0; ig<ngtot; ++ig)
  {
    int j3=(igwk[ig]-1)*3;
    igmapped[i3++]=itmp[j3++];
    igmapped[i3++]=itmp[j3++];
    igmapped[i3++]=itmp[j3++];
  }
  //hid_t h1 = H5Gcreate(h_file,"basis",0);
  //hsize_t dim=1;
  //herr_t ret=H5LTmake_dataset(h1,"num_planewaves",1,&dim,H5T_NATIVE_INT,ngtot);
  hsize_t dims[2];
  dims[0] = ngtot;
  dims[1] = 3;
  //herr_t ret=H5LTmake_dataset(h_kpoint,"gvectors",2,dims,H5T_NATIVE_INT,itmp);
  herr_t ret=H5LTmake_dataset(h_kpoint,"gvectors",2,dims,H5T_NATIVE_INT,igmapped);
  //ret=H5LTmake_dataset(h1,"planewaves",2,dims,H5T_NATIVE_DOUBLE,gcart);

  free(igmapped);
  //H5Gclose(h1);
  }
}

void F77_FUNC_(esh5_write_supercell,ESH5_WRITE_SUPERCELL)(const double* lattice)
{
  hid_t h1 = H5Gcreate(h_file,"supercell",0);
  hsize_t dims[]={3,3};
  herr_t ret=H5LTmake_dataset(h1,"primitive_vectors",2,dims,H5T_NATIVE_DOUBLE,lattice);
  H5Gclose(h1);
}

void F77_FUNC_(esh5_open_atoms,ESH5_OPEN_ATOMS)(const int* nat, const int *nspecies)
{
  h_ptcls = H5Gcreate(h_file,"atoms",0);
  hsize_t dim1=1;
  herr_t ret=H5LTmake_dataset(h_ptcls,"number_of_atoms",1,&dim1,H5T_NATIVE_INT,nat);
  ret=H5LTmake_dataset(h_ptcls,"number_of_species",1,&dim1,H5T_NATIVE_INT,nspecies);
  num_atoms=*nat;
  num_species=*nspecies;
}
void F77_FUNC_(esh5_close_atoms,ESH5_CLOSE_ATOMS)()
{
  H5Gclose(h_ptcls);
}

void F77_FUNC_(esh5_write_species,ESH5_WRITE_SPECIES)(const int* itype
    , const char* sname, const int* length
    , const double* atomic_number, const double* valcharge)
{
  char aname[16];
  sprintf(aname,"species_%i",(*itype)-1);
  hid_t h1 = H5Gcreate(h_ptcls,aname,0);
  hsize_t dim1=1;
  int int_charge = (int) round(*valcharge);
  herr_t ret=H5LTmake_dataset(h1,"valence_charge",1,&dim1,H5T_NATIVE_INT,&int_charge);
  ret=H5LTmake_dataset(h1,"atomic_number",1,&dim1,H5T_NATIVE_INT,atomic_number);

  char species_name[8];
  memcpy(species_name,sname,*length);
  species_name[*length] = '\0' ; 

  hid_t strtype = H5Tcopy (H5T_C_S1);
  ret = H5Tset_size (strtype, (*length)+1); /* create string of length 5 */
  ret=H5LTmake_dataset(h1,"name",1,&dim1,strtype,species_name);

  H5Gclose(h1);
}

void F77_FUNC_(esh5_write_species_ids,ESH5_WRITE_SPEICES_IDS)(const int* ids_in)
{
  int *ids=(int*)malloc(num_atoms*sizeof(int));
  for(int i=0; i<num_atoms; ++i) ids[i]=ids_in[i]-1;
  hsize_t dim1=num_atoms;
  herr_t ret=H5LTmake_dataset(h_ptcls,"species_ids",1,&dim1,H5T_NATIVE_INT,ids);
  free(ids);
}

void F77_FUNC_(esh5_write_positions,ESH5_WRITE_POSITIONS)(const double* r)
{
  hsize_t dims[2];
  dims[0]=num_atoms;
  dims[1]=3;
  herr_t ret=H5LTmake_dataset(h_ptcls,"positions",2,dims,H5T_NATIVE_DOUBLE,r);
}


/** write basisset: number of plane waves, plane wave coefficients
void F77_FUNC_(esh5_write_basis,ESH5_WRITE_BASIS)(const double* g, const int* igtog, const int* ngtot)
{
  int ng=*ngtot;
  int *ig=(int*)malloc(3*ng*sizeof(int));
  for(int i=0,i3=0; i<ng; i++)
  {
    int cur=3*(igtog[i]-1);
    ig[i3++]=(int)g[cur++];
    ig[i3++]=(int)g[cur++];
    ig[i3++]=(int)g[cur++];
  }

  hid_t h_basis = H5Gcreate(h_file,"basis",0);
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_basis, "num_planewaves", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ngtot);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  hsize_t dims[2];
  dims[0] = ng;
  dims[1] = 3;
  dataspace  = H5Screate_simple(2, dims, NULL);
  dataset =  H5Dcreate(h_basis, "planewaves", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ig);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Gclose(h_basis);

  free(ig);
}
 */
void F77_FUNC_(esh5_write_parameters,ESH5_WRITE_PARAMETERS)(
    const int* nelec, const int* nspin, const int* nband, const int* nk,
    const double* ecut, const double* alat, const double* at)
{
  hid_t h_param = H5Gcreate(h_file,"parameters",0);
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_param, "num_spins", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nspin);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "num_electrons", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nelec);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "num_bands", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nband);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "num_twists", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nk);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  int iscomplex=1;
  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "complex_coefficients", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&iscomplex);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "maximum_ecut", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ecut);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  double lattice[9];
  for(int i=0; i<9; i++) lattice[i]=(*alat)*at[i];
  hsize_t dims[]={3,3};
  dataspace  = H5Screate_simple(2, dims, NULL);
  dataset =  H5Dcreate(h_param, "lattice", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,lattice);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Gclose(h_param);
}

///* open mainbody:eigenstates */
//void F77_FUNC_(esh5_open_eigg,ESH5_OPEN_EIGG)()
//{
//  if(h_main>=0) H5Gclose(h_main);
//  h_main = H5Gcreate(h_file,"electrons",0);
//  //h_main = H5Gcreate(h_file,"eigenstates",0);
//}
//
///* close eigenstates */
//void F77_FUNC_(esh5_close_eigg,ESH5_CLOSE_EIGG)()
//{
//  if(h_main>=0) H5Gclose(h_main);
//  h_main=-1;
//}

void F77_FUNC_(esh5_write_rho,ESH5_WRITE_RHO)(const double* rho, const double* rhog, const int* ngm)
{
  hid_t h1 = H5Gcreate(h_ptcls,"density",0);

  hsize_t dim3=3;
  herr_t ret=H5LTmake_dataset(h1,"mesh",1,&dim3,H5T_NATIVE_INT,num_grid);

  hid_t h2 = H5Gcreate(h1,"spin_0",0);
  /* write eigenvector */
  hsize_t dims[3];
  dims[0] = num_grid[0];
  dims[1] = num_grid[1];
  dims[2] = num_grid[2];

  ret=H5LTmake_dataset(h2,"density_r",3,dims,H5T_NATIVE_DOUBLE,rho);
  hsize_t dims_g[2];
  dims_g[0]=*ngm;
  dims_g[1]=2;
  ret=H5LTmake_dataset(h2,"density_g",1,dims_g,H5T_NATIVE_DOUBLE,rhog);
  H5Gclose(h2);
  H5Gclose(h1);
  /*
  hsize_t gdims[2];
  gdims[0]=ngm;
  gdims[1]=2;
  dataspace  = H5Screate_simple(2, gdims, NULL);
  dataset =  H5Dcreate(h_file, "chargedensity_g", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,rhog);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  */

  /* testing with paraview/vtk
  if(is_gamma)
  {
    char vtkname[32];
    sprintf(vtkname,"band%i.vtk",(*ibnd)-1);
    FILE *vtk=fopen(vtkname,"w");

    fprintf(vtk,"# vtk DataFile Version 3.0\n");
    fprintf(vtk,"vtk output\n");
    fprintf(vtk,"ASCII\n");
    fprintf(vtk,"DATASET STRUCTURED_POINTS\n");
    fprintf(vtk,"DIMENSIONS %i %i %i\n",h_ngrid[0],h_ngrid[1],h_ngrid[2]);
    fprintf(vtk,"ORIGIN 0 0 0\n");
    fprintf(vtk,"SPACING 1 1 1\n");
    fprintf(vtk,"\nPOINT_DATA %i\n",h_ngridtot);
    fprintf(vtk,"SCALARS scalars float\n");
    fprintf(vtk,"LOOKUP_TABLE default\n");

    for(int i=0,i2=0; i<h_ngridtot;i+=10)
    { 
      for(int j=0; j<10; j++,i2+=2) fprintf(vtk,"%12.6e ",eigr[i2]*eigr[i2]);
      fprintf(vtk,"\n");
    }
    fprintf(vtk,"\n");
    fclose(vtk);
  }
  */
}

void F77_FUNC_(esh5_write_rhog,ESH5_WRITE_RHOG)(const double* rhog, const int* ngm)
{
}
