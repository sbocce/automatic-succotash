// BSD 3-Clause License
// 
// Copyright (c) 2018, Stefano Boccelli
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>      // std::stringstream
#include <math.h>       /* round, floor, ceil, trunc */
#include <algorithm>    // std::min

using namespace std;

// Prototypes for functions called in MAIN
void   initVDF(const string s_VDFtype, vector< vector<double> > &f);

void   writeSol(const char* filename,vector<double> &x,vector<double> &v,vector< vector<double> > *p_ff);

double limiterFun(double r);

double minmodFun(double v1, double v2, double v3);

void   loadValuesFromFile(const char* filename, const vector<double> &v_xcen, vector<double> &val);

// ------------------------------------------------

int main()
{
  // ##################################
  // ## SIMULATION PARAMETERS
  //
  // ---  Space and time discretizations
  double dt = 1.0e-9;
  size_t Nt = 10000;

  double xmin = -0.0;
  double xmax = 0.04;
  size_t Nx = 100;

  double vmin = -1000;
  double vmax = 20000;
  size_t Nv = 100;
  
  // --- Numerics
  string reconstruction = "mieussens"; // "none" - "linear" - "mieussens"

  // --- Physical quantities and thruster parameters
  double m_part = 21.8e-26;  // [kg] - Xe+ mass of the particle 
  double q_part = 1.602e-19; // [C]  - Xe+ charge 
  double vIoniz = 300.0;     // [m/s] - velocity of particles created by ionization

  string s_VDF_type_init = "zero"; // "zero" or "block"

  // --- Solution writing 
  size_t writeeach = 500;
  // ##################################

  // ###################################
  // BEGIN. 
  // Create vectors for cells interfaces in space and velocity
  vector<double> v_xint(Nx+1, 0.);
  vector<double> v_xcen(Nx, 0.);
  double dx = (xmax - xmin)/Nx;

  for(size_t i = 0; i < Nx+1; ++i)
    v_xint[i] = xmin + dx*i;

  for(size_t i = 0; i < Nx; ++i)
    v_xcen[i] = (xmin+dx/2) + dx*i;

  vector<double> v_vint(Nv+1, 0.);
  vector<double> v_vcen(Nv, 0.);
  double dv = (vmax - vmin)/Nv;

  for(size_t i = 0; i < Nv+1; ++i)
    v_vint[i] = vmin + dv*i;

  for(size_t i = 0; i < Nv; ++i)
    v_vcen[i] = (vmin+dv/2) + dv*i;

  // Initialize solution and fluxes to zero
  // Two containers are used: ff_1, initialized to the initial value of
  // the VDF and ff_2.
  // During the computation, these objects are updated alternately by
  // using two pointers and swapping them. This allows to update the 
  // solution immediately as we compute the flux: no need to save fluxes
  // into (big) flux matrices and access them at a later time.
  // Size is: ff[Nspace cells][Nvel cells]
  
  vector< vector<double> > ff_1(Nx, vector<double>(Nv, 0.)); // Solution 
  vector< vector<double> > ff_2(Nx, vector<double>(Nv, 0.)); // 

  initVDF(s_VDF_type_init, ff_1); // init according to model

  vector< vector<double> > *p_ff    = &ff_1; // points to current solution
  vector< vector<double> > *p_ffNew = &ff_2; // points to updated solution
 
  // --- Initialize additional fields

  vector<double> E_ext(Nx, 0.); // external electric field
  loadValuesFromFile("data/Edata.dat",v_xcen,E_ext);
  vector<double> v_a(Nx, 0.); // acceleration due to external electric field
  for(size_t i = 0; i < Nx; ++i) 
    v_a[i] = q_part*E_ext[i]/m_part;

  // --- Ionization source, happens at the velocity of neutrals
  vector<double> SS(Nx, 0.); // source terms
  loadValuesFromFile("data/Sdata.dat",v_xcen,SS);
  double DeltaVioniz = vIoniz - vmin;
  int ionizPosVel = floor(DeltaVioniz/dv);
  if(ionizPosVel < 0) 
    cerr << "Attention! Ionization source out of velocity domain!" << endl;

  // ============================
  // ===== TIME INTEGRATION =====
  // ============================
  size_t count         = 0;
  size_t filenamecount = 0;

  for(size_t t_id = 0; t_id < Nt; ++t_id) {

    cout << "Solving step: " << t_id << " of " << Nt << endl;

    #pragma omp parallel for
    for(size_t i = 1; i < Nx-1; ++i) { // BOUNDARIEEEESSSSSSS EXCLUDEEEDDD

      double a = v_a[i]; // local acceleration (external forces)

      for(size_t j = 1; j < Nv-1; ++j) {

        //double v = v_vcen[j]; // velocity in current cell
        double v = vmin+dv/2 + dv*j; // Computing is quicker than accessing memory

        double fR_xplus, fL_xplus, fR_xminus, fL_xminus;
        double fR_vplus, fL_vplus, fR_vminus, fL_vminus;

        // Compute space fluxes
        double Fx_plus;
        double Fx_minus;

        // Compute velocity fluxes
        double Fv_plus;
        double Fv_minus; 

        // Build left and right states at the interfaces:
        // Space interfaces are xplus and xminus. L and R are left and right
        // state values. Interfaces along the velocity axis are vplus and vminus.
        // ...and compute numerical fluxes

        if(reconstruction.compare("none") == 0) { // No reconstruction at interface

          fL_xplus  = p_ff->operator[](i)[j];
          fR_xplus  = p_ff->operator[](i+1)[j];
          fL_xminus = p_ff->operator[](i-1)[j];
          fR_xminus = fL_xplus;

          fL_vplus  = p_ff->operator[](i)[j];
          fR_vplus  = p_ff->operator[](i)[j+1];
          fL_vminus = p_ff->operator[](i)[j-1];
          fR_vminus = fL_vplus;

          // Compute space fluxes
          Fx_plus  = 0.5*(v*fL_xplus  + v*fR_xplus  - abs(v)*(fR_xplus  - fL_xplus));
          Fx_minus = 0.5*(v*fL_xminus + v*fR_xminus - abs(v)*(fR_xminus - fL_xminus));
  
          // Compute velocity fluxes
          Fv_plus  = 0.5*(a*fL_vplus  + a*fR_vplus  - abs(a)*(fR_vplus  - fL_vplus));
          Fv_minus = 0.5*(a*fL_vminus + a*fR_vminus - abs(a)*(fR_vminus - fL_vminus));

        // =====  Linear reconstruction  =====
        } else if(reconstruction.compare("linear") == 0) {

          // Near the borders go with a first order method (no reconstruction)
          if( (i <= 1) || (i >= Nx-2) || (j <= 1) || (j >= Nv-2) ) {

            fL_xplus  = p_ff->operator[](i)[j];
            fR_xplus  = p_ff->operator[](i+1)[j];
            fL_xminus = p_ff->operator[](i-1)[j];
            fR_xminus = fL_xplus;
  
            fL_vplus  = p_ff->operator[](i)[j];
            fR_vplus  = p_ff->operator[](i)[j+1];
            fL_vminus = p_ff->operator[](i)[j-1];
            fR_vminus = fL_vplus;

          // Far from boundaries, we can go with second order reconstruction
          } else {  
 
            // Extract useful quantities (not to get crazy)..
            double f_i_j   = p_ff->operator[](i)[j];
            double f_i_jp1 = p_ff->operator[](i)[j+1];
            double f_i_jm1 = p_ff->operator[](i)[j-1];
            double f_ip1_j = p_ff->operator[](i+1)[j];
            double f_im1_j = p_ff->operator[](i-1)[j];
            double f_i_jp2 = p_ff->operator[](i)[j+2];
            double f_i_jm2 = p_ff->operator[](i)[j-2];
            double f_ip2_j = p_ff->operator[](i+2)[j];
            double f_im2_j = p_ff->operator[](i-2)[j];

            // Ratios of consecutive gradients
            double rL_xplus = (f_i_j - f_im1_j)/(f_ip1_j - f_i_j   + 1.0e-25);
            double rR_xplus = (f_ip1_j - f_i_j)/(f_ip2_j - f_ip1_j + 1.0e-25);
  
            double rL_xminus = (f_im1_j - f_im2_j)/(f_i_j - f_im1_j + 1.0e-25);
            double rR_xminus = rL_xplus;
  
            double rL_vplus = (f_i_j - f_i_jm1)/(f_i_jp1 - f_i_j   + 1.0e-25);
            double rR_vplus = (f_i_jp1 - f_i_j)/(f_i_jp2 - f_i_jp1 + 1.0e-25);
  
            double rL_vminus = (f_i_jm1 - f_i_jm2)/(f_i_j - f_i_jm1 + 1.0e-25);
            double rR_vminus = rL_vplus;
  
            // Then reconstruct left and right fluxes applying the slope limiter
            fL_xplus  = f_i_j   + 0.5*limiterFun(rL_xplus)*(f_i_j - f_im1_j); 
            fR_xplus  = f_ip1_j - 0.5*limiterFun(rR_xplus)*(f_ip2_j - f_ip1_j);
  
            fL_xminus = f_im1_j + 0.5*limiterFun(rL_xminus)*(f_im1_j - f_im2_j);
            fR_xminus = f_i_j   - 0.5*limiterFun(rR_xminus)*(f_ip1_j - f_i_j);
  
            fL_vplus  = f_i_j   + 0.5*limiterFun(rL_vplus)*(f_i_j - f_i_jm1);
            fR_vplus  = f_i_jp1 - 0.5*limiterFun(rR_vplus)*(f_i_jp2 - f_i_jp1);
  
            fL_vminus = f_i_jm1 + 0.5*limiterFun(rL_vminus)*(f_i_jm1 - f_i_jm2);
            fR_vminus = f_i_j   - 0.5*limiterFun(rR_vminus)*(f_i_jp1 - f_i_j);

          }

          // Compute space fluxes
          Fx_plus  = 0.5*(v*fL_xplus  + v*fR_xplus  - abs(v)*(fR_xplus  - fL_xplus));
          Fx_minus = 0.5*(v*fL_xminus + v*fR_xminus - abs(v)*(fR_xminus - fL_xminus));
  
          // Compute velocity fluxes
          Fv_plus  = 0.5*(a*fL_vplus  + a*fR_vplus  - abs(a)*(fR_vplus  - fL_vplus));
          Fv_minus = 0.5*(a*fL_vminus + a*fR_vminus - abs(a)*(fR_vminus - fL_vminus));

        // =====  Mieussens-style minmod  =====
        } else if(reconstruction.compare("mieussens") == 0) { // Linear reconstruction

          fL_xplus  = p_ff->operator[](i)[j];
          fR_xplus  = p_ff->operator[](i+1)[j];
          fL_xminus = p_ff->operator[](i-1)[j];
          fR_xminus = fL_xplus;

          fL_vplus  = p_ff->operator[](i)[j];
          fR_vplus  = p_ff->operator[](i)[j+1];
          fL_vminus = p_ff->operator[](i)[j-1];
          fR_vminus = fL_vplus;

          double Phix_plus, Phix_minus, Phiv_plus, Phiv_minus;

          if( (i <= 1) || (i >= Nx-2) || (j <= 1) || (j >= Nv-2) ) { // near the boundaries

            Phix_plus  = 0.0;
            Phix_minus = 0.0;
            Phiv_plus  = 0.0;
            Phiv_minus = 0.0;

          } else { // far from boundaries

            double df_im3half = p_ff->operator[](i-1)[j] - p_ff->operator[](i-2)[j];
            double df_imhalf  = p_ff->operator[](i)[j]   - p_ff->operator[](i-1)[j];
            double df_iphalf  = p_ff->operator[](i+1)[j] - p_ff->operator[](i)[j];
            double df_ip3half = p_ff->operator[](i+2)[j] - p_ff->operator[](i+1)[j];

            double df_jm3half = p_ff->operator[](i)[j-1]   - p_ff->operator[](i)[j-2];
            double df_jmhalf  = p_ff->operator[](i)[j]   - p_ff->operator[](i)[j-1];
            double df_jphalf  = p_ff->operator[](i)[j+1] - p_ff->operator[](i)[j];
            double df_jp3half = p_ff->operator[](i)[j+2] - p_ff->operator[](i)[j+1];

            Phix_plus  = minmodFun(df_imhalf, df_iphalf, df_ip3half);
            Phix_minus = minmodFun(df_im3half, df_imhalf, df_iphalf);

            Phiv_plus  = minmodFun(df_jmhalf, df_jphalf, df_jp3half);
            Phiv_minus = minmodFun(df_jm3half, df_jmhalf, df_jphalf);

          }

          // Compute space fluxes
          Fx_plus  = 0.5*(v*fL_xplus  + v*fR_xplus  - abs(v)*(fR_xplus  - fL_xplus  - Phix_plus));
          Fx_minus = 0.5*(v*fL_xminus + v*fR_xminus - abs(v)*(fR_xminus - fL_xminus - Phix_minus));
  
          // Compute velocity fluxes
          Fv_plus  = 0.5*(a*fL_vplus  + a*fR_vplus  - abs(a)*(fR_vplus  - fL_vplus  - Phiv_plus));
          Fv_minus = 0.5*(a*fL_vminus + a*fR_vminus - abs(a)*(fR_vminus - fL_vminus - Phiv_minus));

        // Reconstruction type unknown
        } else {
          cerr << "ATTENTION! Reconstruction type not recognized. Aborting!" << endl;
          exit(1);
        }

        // Compute ionization source
        double Snow;
        if(j == ionizPosVel) 
        { Snow = SS[i]; }
        else
        { Snow = 0.0; }

        // Write solution in new pointer)
        p_ffNew->operator[](i)[j] = p_ff->operator[](i)[j] - dt/dx*(Fx_plus - Fx_minus)
                                          - dt/dv*(Fv_plus - Fv_minus)
                                          + dt/dv*Snow; 
      }
    }

    // ------ Write solution to file
    count++;
    if(count >= writeeach) { 
      // Create filename
      char filenamestr[512];
      int n = sprintf(filenamestr, "./output/file_%08d.dat", filenamecount);

      // Write file
      writeSol(filenamestr, v_xcen, v_vcen, p_ffNew);

      // Update "writing variables"
      filenamecount++;
      count = 0;
    }

    // ------ SWAP POINTERS
    // Now p_ff will point to the updated solution, while p_ffNew will be used to write the new sol
    swap(p_ff, p_ffNew);

  }

  return 0;
}

// ------------------------------------------------------

void writeSol(const char* filename, vector<double> &x, vector<double> &v, 
                              vector< vector<double> > *p_ff)
{
// Exports the solution into file

  std::ofstream sol_file;
  sol_file.open(filename);

  sol_file << "# Format: x, v, f(x,v)" << endl; // write header
  sol_file << "# Nx = " << x.size() << " - Nv = " << v.size() << endl;

  for(size_t i = 0; i < p_ff->size(); ++i) {
    for(size_t j = 0; j < p_ff->operator[](0).size(); ++j) {
      sol_file << x[i] << "   " << v[j] << "   " << p_ff->operator[](i)[j] << endl;
    }
  }

  sol_file << std::endl;
  sol_file.close(); 
}

// ------------------------------------------------------

double limiterFun(double r)
{
  // return (r + abs(r))/(1 + abs(r)); // van Leer 
  return max(0.0, min(1.0, r)); // minmod
  // return (r*r + r)/(r*r + 1.0); // van Albada
  // return 1.0;
  
}

// ------------------------------------------------------

double minmodFun(double v1, double v2, double v3)
{
  if((v1 >=0.0) && (v2 >= 0.0) && (v3 >= 0.0)) {

    double min_1_2 = min(v1, v2);
    return min(min_1_2, v3);

  } else if ((v1 < 0.0) && (v2 < 0.0) && (v3 < 0.0)) {

    double min_1_2 = min(abs(v1), abs(v2));
    return -min(min_1_2, abs(v3));
    
  } else {
    return 0.0;
  }
}

// ------------------------------------------------------

void initVDF(const string s_VDFtype, vector< vector<double> > &f)
{
// Initializes the VDF "f" according to constant string "s_VDFtype".

  // ======= Initialize all elements to zero
  if(s_VDFtype.compare("zero")==0) { // Do nothing
  
  // ======= Square block of elements to 1, all others to 0
  } else if (s_VDFtype.compare("block")==0) {

    // Create a block equal to one, with all the rest equal to zero
    size_t Nx = f.size();
    size_t Nv = f.at(0).size();

    for(size_t i = round(Nx/3); i < round(Nx/2); ++i) 
      for(size_t j = round(Nv/4); j < round(Nv/2); ++j)
        f[i][j] = 1.0;
 
  // ======= Type not recognized 
  } else {
    cerr << "Attention! Distribution type " << s_VDFtype 
         << " for initialization not recognized. ABORTING" << endl;
    exit(0);
  }
 
}

// ------------------------------------------------------

void loadValuesFromFile(const char* filename, const vector<double> &v_xcen, vector<double> &v)
{
// This function reads data from the file called "filename" and fills the vectors 
// "x" and "v".
// Every line of the provided file is read. If there are empty spaces or comments
// this will most likely generate an error (if lucky) or introduce some garbage
// (if unlucky).

  // 0) Resize array to proper size
  v.resize(v_xcen.size());

  // 1) Read a file
  string lineRead;
  ifstream filein(filename);
  if(filein.is_open() == 0)
  {
    std::cerr << "I could not open the file '" << filename
              << "' Does it exist?\nAborting!\n";
    exit(EXIT_FAILURE);
  }
 
  // 2) Fill vectors
  // (Keep reading until the EOF)

  stringstream nowstream;
  vector<double> x_tmp, v_tmp;
  while(getline(filein,lineRead))
  {
    nowstream.clear();          // clear stream variable nowstream
    nowstream.str(lineRead);    // put newly read line into nowstream

    // ...and unpack it into my variables
    double xNow, vNow;
    nowstream >> xNow >> vNow; 

    x_tmp.push_back(xNow);
    v_tmp.push_back(vNow);
  }

  // 3) Interpolate on the vector v_xcen and save solution
  for(size_t i = 0; i < v_xcen.size(); ++i) {

    double xNow = v_xcen[i]; // current position
    
    // find the locations in the file tabulated positions
    size_t pos = 0;
    double xTMP = x_tmp[0]; // initial tabulated point
    while(xTMP < xNow) { // way too small
      pos += 1;
      xTMP = x_tmp[pos];
    }
    
    size_t pos_1 = pos-1;
    size_t pos_2 = pos;

    // Interpolate value
    v[i] = v_tmp[pos_1] + (v_tmp[pos_2]-v_tmp[pos_1])/(x_tmp[pos_2] - x_tmp[pos_1] + 1.0e-20)*(xNow - x_tmp[pos_1]);

  }
}