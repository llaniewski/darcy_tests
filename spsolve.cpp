#define ARMA_USE_SUPERLU 1
#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
vec arma_spsolve(umat locations, vec values, vec rhs) {
  sp_mat A(locations, values);
  vec x=rhs;
  bool status = spsolve(x, A, rhs);
  if (status == false)  { Rcpp::Rcout << "no solution" << endl; }
  return x;
}

using namespace Rcpp;

class ADVarBool;

class ADVar {
  typedef std::pair< int, size_t > idx_t;
  typedef std::map< idx_t, double > diff_t;
  typedef std::map< std::string, int > name_idx_t;
public:
  typedef std::map< idx_t, ADVar > equations_t;
private:
  static name_idx_t name_idx;
  static int get_name_idx(const std::string& name){
    name_idx_t::iterator it = name_idx.find(name);
    if (it != name_idx.end()) return it->second;
    int k = (int) name_idx.size();
    name_idx[name] = k;
    return k;
  }
  diff_t diff;
  double val;

public:
  static idx_t idx(const std::string& name, const size_t& i) {
    return idx_t(get_name_idx(name), i);
  }
  ADVar(double val_=0): val(val_) {}
  ADVar(const std::string& name, const size_t& i,double val_=0) : val(val_) {
    diff[idx(name, i)] = 1;
  }
  double value() { return val; }
  ADVar operator+ (const ADVar& other) const {
    ADVar ret;
    ret.val = val + other.val;
    for (const auto& [key, value]: diff) ret.diff[key] += value;
    for (const auto& [key, value]: other.diff) ret.diff[key] += value;
    return ret;
  }
  ADVar operator- (const ADVar& other) const {
    ADVar ret;
    ret.val = val - other.val;
    for (const auto& [key, value]: diff) ret.diff[key] += value;
    for (const auto& [key, value]: other.diff) ret.diff[key] -= value;
    return ret;
  }
  ADVar operator* (const ADVar& other) const {
    ADVar ret;
    ret.val = val * other.val;
    for (const auto& [key, value]: diff) ret.diff[key] += value * other.val;
    for (const auto& [key, value]: other.diff) ret.diff[key] += value * val;
    return ret;
  }
  ADVarBool operator == (const ADVar& other) const;
};

class ADVarBool {
  ADVar val;
  int type;
public:
  ADVarBool(const ADVar& val_, const int& type_) : val(val_), type(type_) {};
  operator bool () {
    switch (type) {
    case 0: return val.value() == 0;
    case 1: return val.value() != 0;
    default:
      throw 0;
    }
  }
  operator ADVar() {
    if (type == 0) return val;
    throw 0;
  }
};

ADVarBool ADVar::operator == (const ADVar& other) const {
  return ADVarBool((*this) - other, 0);
}

ADVar::name_idx_t ADVar::name_idx;

std::vector<ADVar> ADvector(const std::string& name, const size_t& size) {
  std::vector<ADVar> ret;
  ret.reserve(size);
  for (size_t i=0; i<size; i++) {
    ret.emplace_back(name,i);
  }
  return ret;
}

class Equation {};

// [[Rcpp::export]]
void test(DataFrame face, DataFrame flux, DataFrame edge) {
  const size_t FACE = face.nrow();
  const size_t FLUX = flux.nrow();
  const size_t EDGE = edge.nrow();
  Rcpp::Rcout << "face: " << FACE << " flux: " << FLUX << " edge: " << EDGE << endl;
  NumericMatrix face_center = face["center"];
  Rcpp::Rcout << "center: " << face_center.nrow() << "x" << face_center.ncol() << endl; 
  auto fp = ADvector("face_pressure",FACE);
  auto ep = ADvector("edge_pressure",EDGE);
  auto fl = ADvector("flux",FLUX);
  ADVar::equations_t eq;
  eq[ADVar::idx("flux_balance",0)] = fl[0] == ep[0];
}
