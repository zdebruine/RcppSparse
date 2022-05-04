#include <Rcpp.h>

// this function receives matrix coordinates "i" and "j" and returns
//  true/false with a 1/16th probability of true. The function is 
//  random uniform regardless of the values of i and j, but the values
//  of i and j are interchangeable, and thus this method is useful for
//  transpose-identical matrix construction
//
// seed must be >= 1. Seeds establish complete reproducibility and
//  different seeds yield unique solutions.
//
// modulus by 16 gives a density of 1/16th, or 6.25%, and is selected
//  due to efficiency
//
template <unsigned int inv_density>
inline bool index_status(unsigned int i, unsigned int j){
  // generate a non-periodic value from "i" and "j", inspired by FPE
  unsigned int x = (i * j) ^ (i & j);
  // Marsaglia's xorshift RNG, https://en.wikipedia.org/wiki/Xorshift
  x ^= x << 12;
  x ^= x >> 5;
  x ^= x >> 8;
  return (x % inv_density) == 0;
}

// simulate a random pattern matrix ("ngCMatrix")
//[[Rcpp::export]]
Rcpp::S4 random_ngCMatrix(const unsigned int nrow, const unsigned int ncol, const unsigned int inv_density = 16, const unsigned int seed = 123){
  if(seed < 1)
    Rcpp::stop("seed must be >= 1");
  if((inv_density & (inv_density - 1)) != 0 || inv_density > 128)
    Rcpp::stop("inverse density must be either 2, 4, 8, 16, 32, 64, or 128");

  std::vector<int> i;
  i.reserve(nrow * ncol / inv_density);
  Rcpp::IntegerVector p(ncol + 1);
  Rcpp::S4 s(std::string("ngCMatrix"));
  if(inv_density == 16){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<16>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  } else if(inv_density == 8){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<8>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  } else if(inv_density == 32){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<32>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  } else if(inv_density == 64){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<64>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  } else if(inv_density == 4){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<4>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  } else if(inv_density == 128){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<128>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  } else if(inv_density == 2){
    for(size_t col = seed; col < (ncol + seed); ++col){
      for(size_t row = seed; row < (nrow + seed); ++row)
        if(index_status<2>(row, col))
          i.push_back(row - seed);
      p[col + 1 - seed] = i.size();
    }
  }

  Rcpp::IntegerVector i_ = Rcpp::wrap(i);
  s.slot("i") = i_;
  s.slot("p") = p;
  s.slot("Dim") = Rcpp::IntegerVector::create(nrow, ncol);
  return s;
}

//[[Rcpp::export]]
Rcpp::S4 random_nsparseVector(const unsigned int length, const unsigned int inv_density = 16, const unsigned int seed = 123){
  if(seed < 1)
    Rcpp::stop("seed must be >= 1");
  if((inv_density & (inv_density - 1)) != 0 || inv_density > 128)
    Rcpp::stop("inverse density must be either 2, 4, 8, 16, 32, 64, or 128");
  
  std::vector<int> i;
  i.reserve(length / inv_density);
  Rcpp::S4 s(std::string("nsparseVector"));
  if(inv_density == 16){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<16>(pos, seed))
        i.push_back(pos);
  } else if(inv_density == 8){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<8>(pos, seed))
        i.push_back(pos);
  } else if(inv_density == 32){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<32>(pos, seed))
        i.push_back(pos);
  } else if(inv_density == 64){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<64>(pos, seed))
        i.push_back(pos);
  } else if(inv_density == 4){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<4>(pos, seed))
        i.push_back(pos);
  } else if(inv_density == 128){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<128>(pos, seed))
        i.push_back(pos);
  } else if(inv_density == 2){
    for(size_t pos = 0; pos < length; ++pos)
      if(index_status<2>(pos, seed))
        i.push_back(pos);
  }
  
  Rcpp::IntegerVector i_ = Rcpp::wrap(i);
  s.slot("i") = i_;
  s.slot("length") = length;
  return s;
}