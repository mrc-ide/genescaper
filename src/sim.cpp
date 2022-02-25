
#include "sim.h"
#include "misc_v13.h"
#include "probability_v17.h"

using namespace std;

//------------------------------------------------
// Simulate from simple Wright-Fisher model
Rcpp::List sim_wrightfisher_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // get inputs from Rcpp format to base C++ format
  int N = rcpp_to_int(args("N"));
  int K = rcpp_to_int(args("K"));
  int L = rcpp_to_int(args("L"));
  vector<int> alleles = rcpp_to_vector_int(args("alleles"));
  double mu = rcpp_to_double(args("mu"));
  vector<vector<double>> mig_mat = rcpp_to_matrix_double(args("mig_mat"));
  vector<int> t_out = rcpp_to_vector_int(args("t_out"));
  int n_t_out = t_out.size();
  int initial_method = rcpp_to_int(args("initial_method"));
  vector<vector<double>> initial_params = rcpp_to_matrix_double(args("initial_params"));
  bool silent = rcpp_to_bool(args["silent"]);
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // objects for storing results
  vector<vector<vector<vector<int>>>> allele_counts_store(n_t_out);
  
  // initialise allele counts in each deme based on input method
  vector<vector<vector<int>>> allele_counts(K);
  if (initial_method == 1) {
    double theta = 2*N*mu;
    for (int k = 0; k < K; ++k) {
      allele_counts[k] = vector<vector<int>>(L);
      for (int l = 0; l < L; ++l) {
        vector<double> p = rdirichlet1(theta / (double)alleles[l], alleles[l]);
        allele_counts[k][l] = rmultinom1(N, p);
      }
    }
  } else if (initial_method == 2) {
    for (int k = 0; k < K; ++k) {
      allele_counts[k] = vector<vector<int>>(L);
      if (k == 0) {
        for (int l = 0; l < L; ++l) {
          vector<double> p = rdirichlet2(initial_params[l]);
          allele_counts[k][l] = rmultinom1(N, p);
        }
      } else {
        allele_counts[k] = allele_counts[0];
      }
    }
  }
  
  
  // create list of all pairwise demes that have positive migration rates
  vector<pair<pair<int, int>, double>> mig_list;
  for (int k1 = 0; k1 < (K - 1); ++k1) {
    for (int k2 = (k1 + 1); k2 < K; ++k2) {
      if (mig_mat[k1][k2] > 0) {
        pair<int, int> deme_pair = {k1, k2};
        mig_list.push_back(make_pair(deme_pair, mig_mat[k1][k2]));
      }
    }
  }
  
  // option to store initial values
  int t_out_next = 0;
  if (t_out[t_out_next] == 0) {
    allele_counts_store[t_out_next] = allele_counts;
    t_out_next++;
  }
  
  // loop through generations
  for (int t = 1; t < (max(t_out) + 1); ++t) {
    
    // report progress
    if (!silent) {
      update_progress(args_progress, "pb", t, max(t_out));
    }
    
    // apply migration
    for (int j = 0; j < mig_list.size(); ++j) {
      
      // draw number of migrants
      int n_mig = rbinom1(N, mig_list[j].second);
      
      // apply migration
      if (n_mig > 0) {
        pair<int, int> deme_pair = mig_list[j].first;
        int k1 = deme_pair.first;
        int k2 = deme_pair.second;
        
        for (int l = 0; l < L; ++l) {
          int tmp1 = n_mig;
          int tmp2 = n_mig;
          int tmp3 = N;
          int tmp4 = N;
          for (int i = 0; i < alleles[l]; ++i) {
           
            // subtract migrants from deme1
            int n_mig_i1 = rhyper1(allele_counts[k1][l][i], tmp3 - allele_counts[k1][l][i], tmp1);
            tmp1 -= n_mig_i1;
            tmp3 -= allele_counts[k1][l][i];
            allele_counts[k1][l][i] -= n_mig_i1;
            
            // subtract migrants from deme2
            int n_mig_i2 = rhyper1(allele_counts[k2][l][i], tmp4 - allele_counts[k2][l][i], tmp2);
            tmp2 -= n_mig_i2;
            tmp4 -= allele_counts[k2][l][i];
            allele_counts[k2][l][i] -= n_mig_i2;
            
            // add migrants
            allele_counts[k1][l][i] += n_mig_i2;
            allele_counts[k2][l][i] += n_mig_i1;
            
          }
        }
        
      }
    }
    
    // loop through demes and loci
    for (int k = 0; k < K; ++k) {
      for (int l = 0; l < L; ++l) {
        
        // apply drift by drawing from multinomial distribution with probabilities
        // given by previous generation frequencies
        int draws_remaining = N;
        int N_remaining = N;
        for (int i = 0; i < (alleles[l] - 1); ++i) {
          int count_new = rbinom1(draws_remaining, allele_counts[k][l][i] / (double)N_remaining);
          N_remaining -= allele_counts[k][l][i];
          draws_remaining -= count_new;
          allele_counts[k][l][i] = count_new;
          
          if (N_remaining == 0) {
            break;
          }
        }
        allele_counts[k][l][alleles[l] - 1] = draws_remaining;
        
        // apply mutation
        if (mu > 0) {
          int n_mut = rbinom1(N, mu);
          if (n_mut > 0) {
            
            // subtract mutants from existing counts
            int draws_remaining = n_mut;
            int N_remaining = N;
            for (int i = 0; i < (alleles[l] - 1); ++i) {
              int n_mut_i = rbinom1(draws_remaining, allele_counts[k][l][i] / (double)N_remaining);
              N_remaining -= allele_counts[k][l][i];
              draws_remaining -= n_mut_i;
              allele_counts[k][l][i] -= n_mut_i;
              
              if (N_remaining == 0) {
                break;
              }
            }
            allele_counts[k][l][alleles[l] - 1] -= draws_remaining;
            
            // add new mutants to allele counts
            draws_remaining = n_mut;
            for (int i = 0; i < (alleles[l] - 1); ++i) {
              int n_mut_i = rbinom1(draws_remaining, 1.0 / double(alleles[l] - i));
              draws_remaining -= n_mut_i;
              allele_counts[k][l][i] += n_mut_i;
            }
            allele_counts[k][l][alleles[l] - 1] += draws_remaining;
            
          }
        }
        
      }  // end loop through loci
    }  // end loop through demes
    
    // store result
    if (t == t_out[t_out_next]) {
      allele_counts_store[t_out_next] = allele_counts;
      if (t_out_next < (n_t_out - 1)) {
        t_out_next++;
      }
    }
    
  }  // end simulation loop
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2 - t1);
  print("simulation completed in", time_span.count(), "seconds\n");
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("allele_counts") = allele_counts_store);
  return ret;
}

//------------------------------------------------
// apply box blur to matrix
Rcpp::NumericMatrix box_blur_cpp(Rcpp::NumericMatrix m, int d) {
  
  // initialise return matrix
  Rcpp::NumericMatrix ret(m.nrow(), m.ncol());
  int nr = ret.nrow();
  int nc = ret.ncol();
  
  // apply box blur to ret matrix
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      
      // get limits of box
      int i_min = (i-d > 1) ? i-d : 1;
      int i_max = (i+d < nr) ? i+d : nr;
      int j_min = (j-d > 1) ? j-d : 1;
      int j_max = (j+d < nc) ? j+d : nc;
      
      // sum values within box
      double v = 0;
      for (int a = i_min; a < i_max; ++a) {
        for (int b = j_min; b < j_max; ++b) {
          v += m(a,b);
        }
      }
      
      // store mean value
      v /= double((i_max - i_min)*(j_max - j_min));
      ret(i,j) = v;
    }
  }
  
  return ret;
}
