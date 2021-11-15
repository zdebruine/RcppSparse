// 11/4/2021 Zach DeBruine (zacharydebruine@gmail.com)
// Please raise issues on github.com/zdebruine/RcppSparse/issues
//
// This header file contains the complete RcppSparse::Matrix class

#ifndef RcppSparse_h
#define RcppSparse_h

#include <RcppCommon.h>

// forward declare class
namespace RcppSparse {
    class Matrix;
    class tabVector;
    class tabMatrix;
}  // namespace RcppSparse

// forward declare Rcpp::as<> Exporter
namespace Rcpp {

    namespace traits {

        template <>
        class Exporter<RcppSparse::Matrix>;

        template <>
        class Exporter<RcppSparse::tabVector>;

        template <>
        class Exporter<RcppSparse::tabMatrix>;
    }  // namespace traits
}  // namespace Rcpp

#include <Rcpp.h>

//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

// ----------- helper function

template <typename T>
T reorder(std::vector<int> idx, T v) {
    T result(idx.size());
    for (size_t i = 0; i < idx.size(); ++i)
        result[i] = v[idx[i]];
    return result;
}

template <typename T>
T join(T a, T b) {
    T c(a.size() + b.size());
    for (int i = 0; i < a.size(); ++i)
        c(i) = a(i);
    for (int j = a.size(), k = 0; k < b.size(); ++j, ++k)
        c(j) = b(k);
    return c;
}

// ----------- RcppSparse::Matrix class

namespace RcppSparse {
    class Matrix {
    public:
        // public member objects
        Rcpp::NumericVector x;
        Rcpp::IntegerVector i, p, Dim;

        // constructors
        Matrix(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector Dim) : x(x), i(i), p(p), Dim(Dim) {}
        Matrix(const Rcpp::S4& s) {
            if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("i") || !s.hasSlot("Dim"))
                throw std::invalid_argument("Cannot construct RcppSparse::Matrix from this S4 object");
            x = s.slot("x");
            i = s.slot("i");
            p = s.slot("p");
            Dim = s.slot("Dim");
        }
        Matrix() {}

        unsigned int rows() { return Dim[0]; }
        unsigned int cols() { return Dim[1]; }
        unsigned int nrow() { return Dim[0]; }
        unsigned int ncol() { return Dim[1]; }
        unsigned int n_nonzero() { return x.size(); };
        Rcpp::NumericVector& nonzeros() { return x; };
        Rcpp::IntegerVector& innerIndexPtr() { return i; };
        Rcpp::IntegerVector& outerIndexPtr() { return p; };

        // create a deep copy of an R object
        Matrix clone() {
            Rcpp::NumericVector x_ = Rcpp::clone(x);
            Rcpp::IntegerVector i_ = Rcpp::clone(i);
            Rcpp::IntegerVector p_ = Rcpp::clone(p);
            Rcpp::IntegerVector Dim_ = Rcpp::clone(Dim);
            return Matrix(x_, i_, p_, Dim_);
        }

        // element lookup at specific index
        double at(int row, int col) const {
            for (int j = p[col]; j < p[col + 1]; ++j) {
                if (i[j] == row)
                    return x[j];
                else if (i[j] > row)
                    break;
            }
            return 0.0;
        }
        double operator()(int row, int col) const { return at(row, col); };
        double operator[](int index) const { return x[index]; };

        // subview clones
        Rcpp::NumericVector operator()(int row, Rcpp::IntegerVector& col) {
            Rcpp::NumericVector res(col.size());
            for (int j = 0; j < col.size(); ++j) res[j] = at(row, col[j]);
            return res;
        };
        Rcpp::NumericVector operator()(Rcpp::IntegerVector& row, int col) {
            Rcpp::NumericVector res(row.size());
            for (int j = 0; j < row.size(); ++j) res[j] = at(row[j], col);
            return res;
        };
        Rcpp::NumericMatrix operator()(Rcpp::IntegerVector& row, Rcpp::IntegerVector& col) {
            Rcpp::NumericMatrix res(row.size(), col.size());
            for (int j = 0; j < row.size(); ++j)
                for (int k = 0; k < col.size(); ++k)
                    res(j, k) = at(row[j], col[k]);
            return res;
        };

        // column access (copy)
        Rcpp::NumericVector col(int col) {
            Rcpp::NumericVector c(Dim[0], 0.0);
            for (int j = p[col]; j < p[col + 1]; ++j)
                c[i[j]] = x[j];
            return c;
        }
        Rcpp::NumericMatrix col(Rcpp::IntegerVector& c) {
            Rcpp::NumericMatrix res(Dim[0], c.size());
            for (int j = 0; j < c.size(); ++j) {
                res.column(j) = col(c[j]);
            }
            return res;
        }

        // row access (copy)
        Rcpp::NumericVector row(int row) {
            Rcpp::NumericVector r(Dim[1], 0.0);
            for (int col = 0; col < Dim[1]; ++col) {
                for (int j = p[col]; j < p[col + 1]; ++j) {
                    if (i[j] == row)
                        r[col] = x[j];
                    else if (i[j] > row)
                        break;
                }
            }
            return r;
        }
        Rcpp::NumericMatrix row(Rcpp::IntegerVector& r) {
            Rcpp::NumericMatrix res(r.size(), Dim[1]);
            for (int j = 0; j < r.size(); ++j) {
                res.row(j) = row(r[j]);
            }
            return res;
        }

        // colSums and rowSums family
        Rcpp::NumericVector colSums() {
            Rcpp::NumericVector sums(Dim[1]);
            for (int col = 0; col < Dim[1]; ++col)
                for (int j = p[col]; j < p[col + 1]; ++j)
                    sums(col) += x[j];
            return sums;
        }
        Rcpp::NumericVector rowSums() {
            Rcpp::NumericVector sums(Dim[0]);
            for (int col = 0; col < Dim[1]; ++col)
                for (int j = p[col]; j < p[col + 1]; ++j)
                    sums(i[j]) += x[j];
            return sums;
        }
        Rcpp::NumericVector colMeans() {
            Rcpp::NumericVector sums = colSums();
            for (int i = 0; i < sums.size(); ++i)
                sums[i] = sums[i] / Dim[0];
            return sums;
        };
        Rcpp::NumericVector rowMeans() {
            Rcpp::NumericVector sums = rowSums();
            for (int i = 0; i < sums.size(); ++i)
                sums[i] = sums[i] / Dim[1];
            return sums;
        };

        // crossprod
        Rcpp::NumericMatrix crossprod() {
            Rcpp::NumericMatrix res(Dim[1], Dim[1]);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int col1 = 0; col1 < Dim[1]; ++col1) {
                for (int col2 = col1; col2 < Dim[1]; ++col2) {
                    if (col1 == col2) {
                        for (int j = p[col1]; j < p[col1 + 1]; ++j)
                            res(col1, col1) += x[j] * x[j];
                    } else {
                        int col1_ind = p[col1], col1_max = p[col1 + 1];
                        int col2_ind = p[col2], col2_max = p[col2 + 1];
                        while (col1_ind < col1_max && col2_ind < col2_max) {
                            int row1 = i[col1_ind];
                            int row2 = i[col2_ind];
                            if (row1 == row2) {
                                res(col1, col2) += x[col1_ind] * x[col2_ind];
                                ++col1_ind;
                                ++col2_ind;
                            } else if (row1 < row2) {
                                do {
                                    ++col1_ind;
                                } while (i[col1_ind] < row2 && col1_ind < col1_max);
                            } else if (row2 < row1) {
                                do {
                                    ++col2_ind;
                                } while (i[col2_ind] < row1 && col2_ind < col2_max);
                            }
                        }
                        res(col2, col1) = res(col1, col2);
                    }
                }
            }
            return res;
        }

        // return indices of rows with nonzero values for a given column
        // this function is similar to Rcpp::Range, but unlike Rcpp::Range it is thread-safe
        std::vector<unsigned int> InnerIndices(int col) {
            std::vector<unsigned int> v(p[col + 1] - p[col]);
            for (int j = 0, it = p[col]; it < p[col + 1]; ++j, ++it)
                v[j] = (unsigned int)i[it];
            return v;
        }

        // return indices of rows with zeros values for a given column
        std::vector<unsigned int> emptyInnerIndices(int col) {
            // first get indices of non-zeros
            std::vector<unsigned int> nonzeros = InnerIndices(col);
            std::vector<unsigned int> all_vals(Dim[0]);
            std::iota(all_vals.begin(), all_vals.end(), 0);
            std::vector<unsigned int> zeros;
            std::set_difference(all_vals.begin(), all_vals.end(), nonzeros.begin(), nonzeros.end(),
                                std::inserter(zeros, zeros.begin()));
            return zeros;
        }

        // const column iterator
        class InnerIterator {
        public:
            InnerIterator(Matrix& ptr, int col) : ptr(ptr), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1]) {}
            operator bool() const { return (index < max_index); }
            InnerIterator& operator++() {
                ++index;
                return *this;
            }
            const double& value() const { return ptr.x[index]; }
            int row() const { return ptr.i[index]; }
            int col() const { return col_; }

        private:
            Matrix& ptr;
            int col_, index, max_index;
        };

        // equivalent to the "Forward Range" concept in two boost::ForwardTraversalIterator
        // iterates over non-zero values in `ptr.col(col)` at rows in `s`
        // `s` must be sorted in ascending order
        class InnerIteratorInRange {
        public:
            InnerIteratorInRange(Matrix& ptr, int col, std::vector<unsigned int>& s) : ptr(ptr), s(s), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1] - 1), s_max_index(s.size() - 1) {
                // decrement max_index and s_max_index to last case where ptr.i intersects with s
                while ((unsigned int)ptr.i[max_index] != s[s_max_index] && max_index >= index && s_max_index >= 0)
                    s[s_max_index] > (unsigned int)ptr.i[max_index] ? --s_max_index : --max_index;
                // increment index to the first case where ptr.i intersects with s
                while ((unsigned int)ptr.i[index] != s[s_index] && index <= max_index && s_index <= s_max_index)
                    s[s_index] < (unsigned int)ptr.i[index] ? ++s_index : ++index;
            }
            operator bool() const { return (index <= max_index && s_index <= s_max_index); }
            InnerIteratorInRange& operator++() {
                ++index;
                ++s_index;
                while (index <= max_index && s_index <= s_max_index && (unsigned int)ptr.i[index] != s[s_index])
                    s[s_index] < (unsigned int)ptr.i[index] ? ++s_index : ++index;
                return *this;
            }
            const double& value() const { return ptr.x[index]; }
            int row() const { return ptr.i[index]; }
            int col() const { return col_; }

        private:
            Matrix& ptr;
            const std::vector<unsigned int>& s;
            int col_, index, max_index, s_max_index, s_index = 0, s_size;
        };

        // iterates over non-zero values in ptr.col(col) not at rows in s_
        // basically, turn this into InnerIteratorInRange by computing a vector `s` of non-intersecting
        //    non-zero rows in ptr at time of initialization
        // s must be sorted in ascending order
        class InnerIteratorNotInRange {
        public:
            InnerIteratorNotInRange(Matrix& ptr, int col, std::vector<unsigned int>& s_) : ptr(ptr), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1] - 1) {
                s = std::vector<unsigned int>(ptr.p[col_ + 1] - ptr.p[col]);
                if (s.size() > 0) {
                    for (int j = 0, it = ptr.p[col_]; it < ptr.p[col_ + 1]; ++j, ++it)
                        s[j] = (unsigned int)ptr.i[it];
                    if (s_.size() > 0) {
                        // remove intersecting values in s_ from s
                        unsigned int si = 0, s_i = 0, z_i = 0;
                        std::vector<unsigned int> z = s;
                        while (si < s.size()) {
                            if (s_i > s_.size() || s_[s_i] > s[si]) {
                                z[z_i] = s[si];
                                ++si;
                                ++z_i;
                            } else if (s_[s_i] == s[si]) {
                                ++si;
                                ++s_i;
                            } else
                                ++s_i;
                        }
                        z.resize(z_i);
                        s = z;
                    }
                }
                s_max_index = s.size() - 1;

                // decrement max_index and s_max_index to last case where ptr.i intersects with s
                while ((unsigned int)ptr.i[max_index] != s[s_max_index] && max_index >= index && s_max_index >= 0)
                    s[s_max_index] > (unsigned int)ptr.i[max_index] ? --s_max_index : --max_index;
                // increment index to the first case where ptr.i intersects with s
                while ((unsigned int)ptr.i[index] != s[s_index] && index <= max_index && s_index <= s_max_index)
                    s[s_index] < (unsigned int)ptr.i[index] ? ++s_index : ++index;
            }
            operator bool() const { return (index <= max_index && s_index <= s_max_index); }
            InnerIteratorNotInRange& operator++() {
                ++index;
                ++s_index;
                while (index <= max_index && s_index <= s_max_index && (unsigned int)ptr.i[index] != s[s_index])
                    s[s_index] < (unsigned int)ptr.i[index] ? ++s_index : ++index;
                return *this;
            }
            const double& value() const { return ptr.x[index]; }
            int row() const { return ptr.i[index]; }
            int col() const { return col_; }

        private:
            Matrix& ptr;
            std::vector<unsigned int> s;
            int col_, index, max_index, s_max_index, s_index = 0, s_size;
        };

        // const row iterator
        class InnerRowIterator {
        public:
            InnerRowIterator(Matrix& ptr, int j) : ptr(ptr) {
                for (; index < ptr.Dim[1]; ++index) {
                    if (ptr.i[index] == j) break;
                }
                for (int r = 0; r < ptr.i.size(); ++r)
                    if (ptr.i[r] == j) max_index = r;
            }
            operator bool() const { return index <= max_index; };
            InnerRowIterator& operator++() {
                ++index;
                for (; index <= max_index; ++index) {
                    if (ptr.i[index] == row_) break;
                }
                return *this;
            };
            int col() {
                int j = 0;
                for (; j < ptr.p.size(); ++j) {
                    if (ptr.p[j] > index) break;
                }
                return j;
            };
            int row() { return row_; }
            double& value() const { return ptr.x[index]; };

        private:
            Matrix& ptr;
            int row_ = 0, index = 0, max_index = 0;
        };

        // number of nonzeros in a column
        unsigned int InnerNNZs(int col) {
            return p[col + 1] - p[col];
        }

        // is approximately symmetric
        bool isAppxSymmetric() {
            if (Dim[0] == Dim[1]) {
                InnerIterator col_it(*this, 0);
                InnerRowIterator row_it(*this, 0);
                while (++col_it && ++row_it) {
                    if (col_it.value() != row_it.value())
                        return false;
                }
                return true;
            }
            return false;
        }

        Matrix transpose() {
            Rcpp::S4 s(std::string("dgCMatrix"));
            s.slot("i") = i;
            s.slot("p") = p;
            s.slot("x") = x;
            s.slot("Dim") = Dim;
            Rcpp::Environment base = Rcpp::Environment::namespace_env("Matrix");
            Rcpp::Function t_r = base["t"];
            Rcpp::S4 At = t_r(Rcpp::_["x"] = s);
            return Matrix(At);
        };

        Rcpp::S4 wrap() {
            Rcpp::S4 s(std::string("dgCMatrix"));
            s.slot("x") = x;
            s.slot("i") = i;
            s.slot("p") = p;
            s.slot("Dim") = Dim;
            return s;
        }
    };

    // start writing R S4 methods
    // objective: R S4 tabMatrix class, with S4 tabVector class to power vectorized updates
    // write method for compressing tabVectors. R S4 does the subsetting work for insertions, and the bulk of the work for operations
    // this is an R package first and foremost, which happens to have an Rcpp

    // ----------- RcppSparse::tabMatrix class

    // constructed from std::vector<int> for thread safety (parallelized constructors are very useful)
    class tabVector {
    public:
        Rcpp::IntegerVector x;
        int length = 0;

        // constructors
        tabVector(Rcpp::IntegerVector x, int length) : x(x), length(length) {};
        tabVector(Rcpp::IntegerVector Ax, Rcpp::IntegerVector Ai, int length) : length(length) {
            if (Ax.size() > 0) {
                // table up A.x
                Rcpp::IntegerVector vals, counts;
                vals.push_back(Ax[0]);
                counts.push_back(1);
                for (int i = (0 + 1); i < Ax.size(); ++i) {
                    // check if Ax value is in array
                    bool add = true;
                    for (int j = 0; j < vals.size(); ++j) {
                        if (vals[j] == Ax[i]) {
                            ++counts[j];
                            add = false;
                            break;
                        }
                    }
                    if (add) {
                        vals.push_back(Ax[i]);
                        counts.push_back(1);
                    }
                }
                // get sort index for vals and counts
                std::vector<int> idx(vals.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::stable_sort(idx.begin(), idx.end(),
                                 [&vals](int i1, int i2) { return vals[i1] < vals[i2]; });

                // reorder vals and counts based on sort index
                vals = reorder(idx, vals);
                counts = reorder(idx, counts);

                int num_vals = std::accumulate(counts.begin(), counts.end(), (int)0);
                num_vals += vals.size();

                // return row indices in "i" ordered by values
                x = Rcpp::IntegerVector(num_vals);
                int ind = 0;
                for (auto value : vals) {
                    x[ind] = value * -1;
                    ++ind;
                    // look for value in column
                    for (int j = 0; j < Ax.size(); ++j) {
                        if (Ax[j] == value) {
                            x[ind] = Ai[j];
                            ++ind;
                        }
                    }
                }
            }
        }
        tabVector(tabMatrix A, int col) {
            x = A.x[Rcpp::Range(A.p[col], A.p[col + 1] - 1)];
            length = A.rows();
        }
        tabVector() {};

        // const iterator
        class Iterator {
        public:
            Iterator(tabVector& ptr) : ptr(ptr) {
                max_index_ = ptr.x.size();
                // catch for an empty vector
                if (index_ < max_index_) {
                    value_ = std::abs(ptr.x[index_]);
                    ++index_;
                }
            }

            operator bool() const {
                return (index_ < max_index_);
            }

            Iterator& operator++() {
                ++index;
                if (ptr.x[index_] < 0) {
                    value_ = std::abs(ptr.x[index_]);
                    ++index_;
                }
                return *this;
            }

            const int value() const {
                return value_;
            }

            void nextValue() const {
                ++index_;
                while (index_ < max_index_) {
                    if (ptr.x[index_] < 0) {
                        value_ = std::abs(ptr.x[index]);
                        ++index_;
                        return;
                    }
                    ++index_;
                }
            }

            void prevValue() const {
                --index_;
                while (index_ > 0) {
                    if (ptr.x[index_] < 0) {
                        value_ = std::abs(ptr.x[index]);
                        --index_;
                        return;
                    }
                }
            }

            const int row() const {
                return ptr.x[index_];
            }

            const int index() const {
                return index_;
            }

        private:
            tabVector& ptr;
            int row_, value_, index_ = 0, max_index_ = 0;
        };

        tabVector clone() {
            Rcpp::IntegerVector x_ = Rcpp::clone(x);
            int length_ = length;
            return tabVector(x_, length_);
        }

        int size() { return x.size(); }

        // COERCION METHODS

        template <>
        Rcpp::IntegerVector as<Rcpp::IntegerVector>() {
            Rcpp::IntegerVector v(length);
            int curr_value;
            for (int i = 0; i < x.size(); ++i) {
                if (x[i] < 0) {
                    curr_value = std::abs(x[i]);
                    ++i;
                }
                v[x[i]] = curr_value;
            }
            return v;
        }

        // wrap
        template <>
        Rcpp::S4 as<Rcpp::S4>() {
            Rcpp::S4 s(std::string("tabVector"));
            s.slot("x") = x;
            s.slot("length") = length;
            return s;
        }

        // SUBSETTING
        template <typename T>
        int operator[](T ind) {
            for (int j = 0; j < x.size(); ++j) {
                if (x[j] == (int)ind) {
                    do { --j; } while (x[j] > 0);
                    return std::abs(x[j]);
                }
            }
            return 0;
        }

        tabVector operator[](Rcpp::IntegerVector ind) {
            Rcpp::IntegerVector x_, p_;
            for (int i : ind) {
                for (int j = 0; j < x.size(); ++j) {
                    if (x[j] == i) {
                        p_.push_back(i);
                        // backtrack to last value
                        do { --j; } while (x[j] > 0);
                        x_.push_back(std::abs(x[j]));
                    }
                }
            }
            return tabVector(x_, p_, ind.size());
        }

        template <typename T>
        tabVector range(T start_row, T end_row) {
            Rcpp::IntegerVector range_((int)end_row - (int)start_row);
            for (int i = (int)start_row; i <= (int)end_row; ++i)
                range_(i) = i;
            return *this(range_);
        }

        tabVector operator()(Rcpp::IntegerVector ind) { return *this[](ind); }
        template <typename T> int operator()(T ind) { return *this[](ind); }

        // cross-product against dense vector
        int cross(Rcpp::IntegerVector y) {
            int value = 0;
            for (Iterator it(*this); it; ++it)
                value += it.value() * y(it.row());
            return value;
        }
        double cross(Rcpp::NumericVector y) {
            double value = 0;
            for (Iterator it(*this); it; ++it)
                value += (double)it.value() * y(it.row());
            return value;
        }
        int cross(tabVector y) {
            Rcpp::IntegerVector y_ = y.as_vector();
            return cross(y_);
        }
        int cross2(tabVector y) {
            if (x.size() == 0 || y.size() == 0) return 0;
            int v = 0, y_value = y.x[0];
            Rcpp::IntegerVector p_y(1), p_x(1); // value index pointers
            // within each value pair in "x" and "y", we do a forward traversal iterator
            // on the first pass, we index positions of new values in "y" as we go
            for (int i = 2; i < x.size(); ++i)
                if (x[i] < 0)
                    p_x.push_back(i);

            for (int it_x : p_x) {
                int x_value = x[it_x];
                ++it_x;
                int y_value = y[0];
                int it_y = 1;
                while (it_x < x.size() && it_y < y.size()) {
                    if (x[it_x] < y.x[it_y]) {
                        ++it_x;
                        if (x[it_x] < 0) {
                            it_x = 1;
                        }
                    } else if (y.x[it_y] < x[it_x]) {
                        ++it_y;
                        if (y.x[it_y] < 0) {
                            p_y.push_back(it_y);
                            y_value = y.x[it_y];
                            ++it_y;
                            it_x = 1;
                        }
                    } else if (x[it_x] == y[it_y]) {
                        v += x_value * y_value;
                        ++it_x;
                        ++it_y;
                        if (it_y < y.size()) {
                            if (y.x[it_y] < 0) {
                                p_y.push_back(it_y);
                                y_value = y.x[it_y];
                                ++it_y;
                            }
                            if (it_x < x.size()) {
                                if (x[it_x] < 0) {
                                    it_x = 1;
                                }
                            }
                        }
                    }
                }
            }
            // on the second...n pass, we use the previously determined index positions to avoid linear search
            //   and instead use the standard binary index search lookup
            // passes are repeated until we have reached the end of the x_vector

            // loop through each value pair to find rows in common
            for (Iterator it(*this); it; ++it) {
                for (Iterator it_y(*this); it; ++it) {
                    // 
                }
            }
            Rcpp::IntegerVector y_ = y.as_vector();
            return v;
        }

        // dot-product against dense vector
        tabVector dot(Rcpp::IntegerVector y) {
            Rcpp::IntegerVector v(length);
            for (Iterator it(*this); it; ++it)
                v(it.row()) = y(it.row()) * it.value();
            return tabVector(v, length);
        }
        tabVector dot(tabVector y) {
            return dot(y.as_vector());
        }
        Rcpp::NumericVector dot(Rcpp::NumericVector y) {
            Rcpp::NumericVector v(y.size());
            for (Iterator it(*this); it; ++it)
                v(it.row()) = y(it.row()) * (double)it.value();
            return v;
        }

        // cumulative functions
        int sum() {
            if (length == 0) return 0;
            int curr_val = x[0], v = 0, count = 1;
            for (int i = 2; i < x.size(); ++i) {
                if (x[i] < 1) {
                    v += count * curr_val;
                    curr_val = x[i];
                    ++i;
                    count = 1;
                } else ++count;
            }
            return std::abs(v);
        }

        int rows() { return length; }

        int max() {
            int max_val = x[0];
            for (int val : x) {
                if (val < max_val)
                    max_val = val;
            }
            return std::abs(max_val);
        }

        int min() {
            int min_val = x[0];
            for (int val : x) {
                if (val < 0) {
                    if (val == -1) return 1;
                    if (val > min_val) {
                        min_val = val;
                    }
                }
            }
            return std::abs(min_val);
        }

        // ------------ OPERATIONS

        // these operations coerce to dense
        Rcpp::IntegerVector operator+(int y) {
            Rcpp::IntegerVector v = as_vector();
            for (int i = 0; i < v.size(); ++i) v[i] += y;
            return v;
        }
        Rcpp::IntegerVector operator+(unsigned int y) { return (*this + (int)y); }
        Rcpp::NumericVector operator+(double y) {
            Rcpp::IntegerVector v = as_vector();
            Rcpp::NumericVector v_(v);
            for (int i = 0; i < v_.size(); ++i) v_[i] += y;
            return v_;
        }
        Rcpp::NumericVector operator+(Rcpp::NumericVector y) {
            Rcpp::IntegerVector v = as_vector();
            Rcpp::NumericVector v_(v);
            for (int i = 0; i < v_.size(); ++i) v_[i] += y[i];
            return v_;
        }
        Rcpp::IntegerVector operator+(Rcpp::IntegerVector y) {
            Rcpp::IntegerVector v = as_vector();
            for (int i = 0; i < v.size(); ++i) v[i] += y[i];
            return v;
        }
        Rcpp::NumericVector operator/(Rcpp::NumericVector y) {
            Rcpp::IntegerVector v = as_vector();
            Rcpp::NumericVector v_(v);
            for (int i = 0; i < v_.size(); ++i) v_[i] /= y[i];
            return v_;
        }
        Rcpp::NumericVector operator/(double y) {
            Rcpp::IntegerVector v = as_vector();
            Rcpp::NumericVector v_(v);
            for (int i = 0; i < v_.size(); ++i) v_[i] /= y;
            return v_;
        }
        Rcpp::IntegerVector operator-(int y) { return *this + (-y); }
        Rcpp::IntegerVector operator-(unsigned int y) { return *this + (-(int)y); }
        Rcpp::NumericVector operator-(double y) { return *this + (-y); }
        Rcpp::IntegerVector operator-(Rcpp::IntegerVector y) { return *this + (-y); }
        Rcpp::NumericVector operator-(Rcpp::NumericVector y) { return *this + (-y); }
        Rcpp::NumericVector operator/(double y) { Rcpp::NumericVector v_(as_vector()); return v_ / y; }

        // integer multiplication is special!
        void operator*=(int y) {
            for (int& x_ : x)
                if (x_ < 0)
                    x_ *= y;
        }
        void operator*=(unsigned int y) { return (*this *= (int)y); }
        tabVector operator*(int y) {
            tabVector x_ = clone();
            x_ *= y;
            return x_;
        }
        tabVector operator*(unsigned int y) { return (*this * (unsigned int)y); }
        tabVector operator*(Rcpp::IntegerVector y) { return dot(y); }
        Rcpp::NumericVector operator*(Rcpp::NumericVector y) { return dot(y); }

        // operations between two tabVectors
        // rhs is always coerced to dense
        Rcpp::IntegerVector operator+(tabVector y) { return *this + y.as_vector(); }
        Rcpp::IntegerVector operator-(tabVector y) { return *this - y.as_vector(); }
        Rcpp::NumericVector operator/(tabVector y) { Rcpp::NumericVector y_(y.as_vector()); return *this / y_; }
        tabVector operator*(tabVector y) { return dot(y.as_vector()); }
        void operator*=(tabVector y) { *this = dot(y); }
        void operator*=(Rcpp::IntegerVector y) { *this = dot(y); }

        // LOGICAL COMPARISONS
        // these operations coerce to dense
        Rcpp::LogicalVector operator!=(int y) {
            if (y == 0) {
                Rcpp::LogicalVector v(length, false);
                for (int x_ : x)
                    if (x_ > 0)
                        v(x_) = true;
                return v;
            } else {
                Rcpp::IntegerVector x_ = as_vector();
                return x_ != y;
            }
        }
        Rcpp::LogicalVector operator==(int y) {
            Rcpp::LogicalVector v = (*this != y);
            return !v;
        }

        Rcpp::LogicalVector operator<(int y) {
            if (y <= 0) {
                return Rcpp::LogicalVector(length, false);
            } else if (y == 1) {
                return (*this == 0);
            } else {
                Rcpp::IntegerVector x_ = as_vector();
                return x_ < y;
            }
        }

        Rcpp::LogicalVector operator>(int y) {
            if (y == 0) return (*this != 0);
            else return !(*this < y);
        }

        Rcpp::LogicalVector operator<=(int y) { return (*this < (y + 1)); }
        Rcpp::LogicalVector operator>=(int y) { return (*this > (y - 1)); }
        Rcpp::LogicalVector operator!=(unsigned int y) { return *this != (int)y; }
        Rcpp::LogicalVector operator==(unsigned int y) { return *this == (int)y; }
        Rcpp::LogicalVector operator<(unsigned int y) { return *this < (int)y; }
        Rcpp::LogicalVector operator>(unsigned int y) { return *this > (int)y; }
        Rcpp::LogicalVector operator<=(unsigned int y) { return *this <= (int)y; }
        Rcpp::LogicalVector operator>=(unsigned int y) { return *this >= (int)y; }

        int n_nonzeros() {
            int v = 0;
            for (int val : x) {
                if (val >= 0) ++v;
            }
            return v;
        }
    };

    class tabMatrix {
    public:
        // public member objects
        Rcpp::IntegerVector x, p, Dim;

        // constructors
        tabMatrix(Rcpp::IntegerVector x, Rcpp::IntegerVector p, Rcpp::IntegerVector Dim) : x(x), p(p), Dim(Dim) {}
        tabMatrix(const Rcpp::S4& s) {
            if (!s.hasSlot("x") || !s.hasSlot("i") || !s.hasSlot("Dim"))
                throw std::invalid_argument("Cannot construct RcppSparse::tabMatrix from this S4 object");
            x = s.slot("x");
            p = s.slot("p");
            Dim = s.slot("Dim");
        }
        tabMatrix(Matrix A) {
            std::vector<std::vector<int>> x_(A.cols());
            std::vector<int> p__ = Rcpp::as<std::vector<int>>(A.p);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (unsigned int i = 0; i < A.cols(); ++i) {
                int start = A.p[i], stop = A.p[i + 1];
                if (stop > start) {
                    // table up A.x
                    std::vector<int> vals, counts;
                    vals.push_back(A.x[start]);
                    counts.push_back(1);
                    for (int i = (start + 1); i < stop; ++i) {
                        // check if A.x value is in array
                        bool add = true;
                        for (size_t j = 0; j < vals.size(); ++j) {
                            if (vals[j] == A.x[i]) {
                                ++counts[j];
                                add = false;
                                break;
                            }
                        }
                        if (add) {
                            vals.push_back(A.x[i]);
                            counts.push_back(1);
                        }
                    }
                    // get sort index for vals and counts
                    std::vector<int> idx(vals.size());
                    std::iota(idx.begin(), idx.end(), 0);
                    std::stable_sort(idx.begin(), idx.end(),
                                     [&vals](int i1, int i2) { return vals[i1] < vals[i2]; });

                    // reorder vals and counts based on sort index
                    vals = reorder(idx, vals);
                    counts = reorder(idx, counts);

                    int num_vals = std::accumulate(counts.begin(), counts.end(), (int)0);
                    num_vals += vals.size();

                    // return row indices in "i" ordered by values
                    std::vector<int> result(num_vals);
                    int ind = 0;
                    for (auto value : vals) {
                        result[ind] = value * -1;
                        ++ind;
                        // look for value in column
                        for (int j = start; j < stop; ++j) {
                            if (A.x[j] == value) {
                                result[ind] = A.i[j];
                                ++ind;
                            }
                        }
                    }
                    x_[i] = result;
                }
            }
            p = Rcpp::IntegerVector(A.cols() + 1);
            for (size_t i = 0; i < x_.size(); ++i)
                p[i + 1] = x_[i].size() + p[i];
            x = Rcpp::IntegerVector(p[A.cols()]);
            for (unsigned int i = 0, idx = 0; i < A.cols(); ++i)
                for (unsigned int j = 0; j < x_[i].size(); ++j, ++idx)
                    x[idx] = x_[i][j];
            Dim = A.Dim;
        }
        tabMatrix(Rcpp::NumericMatrix A) {
            Dim = Rcpp::IntegerVector(2);
            p = Rcpp::IntegerVector(1);
            Dim[0] = A.rows();
            Dim[1] = A.cols();
            for (int i = 0; i < A.cols(); ++i) {
                x.push_back(tabVector(A.column(i)).x);
                p.push_back(x.size());
            }
        }
        tabMatrix() {}

        unsigned int rows() { return Dim[0]; }
        unsigned int cols() { return Dim[1]; }
        unsigned int nrow() { return Dim[0]; }
        unsigned int ncol() { return Dim[1]; }

        // const column iterator
        class InnerIterator {
        public:
            InnerIterator(tabMatrix& ptr, int col) : ptr(ptr), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1]) {
                // catch for empty columns
                if (index < max_index) {
                    value_ = std::abs(ptr.x[index]);
                    ++index;
                }
            }

            operator bool() const {
                return (index < max_index);
            }

            InnerIterator& operator++() {
                ++index;
                if (ptr.x[index] < 0) {
                    value_ = std::abs(ptr.x[index]);
                    ++index;
                }
                return *this;
            }

            const int value() const {
                return value_;
            }

            const int row() const {
                return ptr.x[index];
            }

            const int col() const {
                return col_;
            }

        private:
            tabMatrix& ptr;
            int col_, row_, value_, index, max_index;
        };

        // const row iterator
        // terribly slow, only used for printing objects to R console
        // this is so bad that it is not even publicly documented to avoid usage
        class RowIterator {
        public:
            RowIterator(tabMatrix& ptr, int row) : ptr(ptr), row_(row) {
                // find first index in row
                max_index = ptr.x.size();
                while (index < max_index && ptr.x[index] != row_)
                    ++index;

                if (index < max_index) {
                    // there are non-zero values in this row

                    // find current value
                    int val_index = index - 1;
                    while (ptr.x[val_index] >= 0)
                        --val_index;
                    value_ = std::abs(ptr.x[val_index]);

                    // find column in which first non-zero resides
                    while (ptr.p[col_ + 1] < index)
                        ++col_;

                    // find last index in row
                    if (ptr.x[max_index] != row_) {
                        while (--max_index) {
                            if (ptr.x[max_index] == row_)
                                break;
                        }
                        ++max_index;
                    }
                }
            }

            operator bool() const {
                return (index < max_index);
            }

            RowIterator& operator++() {
                ++index;
                if (index < max_index) {
                    while (ptr.x[index] != row_) ++index;
                    // now find the corresponding value
                    int val_index = index - 1;
                    while (ptr.x[val_index] >= 0)
                        --val_index;
                    value_ = std::abs(ptr.x[val_index]);

                    // find corresponding column
                    while (ptr.p[col_ + 1] < index) ++col_;
                }
                return *this;
            }

            const int value() const {
                return value_;
            }

            const int row() const {
                return row_;
            }

            const int col() const {
                return col_;
            }

        private:
            tabMatrix& ptr;
            int col_ = 0, row_, value_ = -1, index = 0, max_index = 0;
        };

        tabMatrix clone() {
            Rcpp::IntegerVector x_ = Rcpp::clone(x);
            Rcpp::IntegerVector p_ = Rcpp::clone(p);
            Rcpp::IntegerVector Dim_ = Rcpp::clone(Dim);
            return tabMatrix(x_, p_, Dim_);
        }

        tabVector col(int col) { return tabVector(A, col); }

        template <typename T>
        tabVector row(T i) {
            Rcpp::IntegerVector x_, i_;
            for (RowIterator it(*this, (int)i); it; ++it) {
                x_.push_back(it.value());
                i_.push_back(it.col());
            }
            return tabVector(x_, i_, Dim[0]);
        }

        tabMatrix col(Rcpp::IntegerVector cols) {
            Rcpp::IntegerVector x_, p_(1);
            for (int col : cols) {
                Rcpp::IntegerVector x__ = x[Rcpp::Range(p[col], p[col + 1] - 1)];
                x_ = join(x_, x__);
                p_.push_back(x_.size());
            }
            Rcpp::IntegerVector Dim_ = Dim;
            Dim_(1) = cols.size();
            return tabMatrix(x_, p_, Dim_);
        }

        template <typename T>
        int at(T row, T col) {
            for (int i = p[(int)col]; i < p[(int)col + 1]; ++i) {
                if (x[i] == (int)row) {
                    // backtrack to last value
                    do { --i; } while (x[i] > 0);
                    return std::abs(x[i]);
                }
            }
            return 0;
        }

        template <typename T>
        int operator()(T row, T col) {
            return at(row, col);
        }

        template <typename T>
        tabVector operator()(T row, Rcpp::IntegerVector cols) {
            Rcpp::IntegerVector x_, p_;
            for (int col : cols) {
                // find index of row in col
                for (int j = p[col]; j < p[col + 1]; ++j) {
                    if (x[j] == (int)row) {
                        p_.push_back(col);
                        // backtrack to last value
                        do { --j; } while (x[j] > 0);
                        x_.push_back(std::abs(x[j]));
                        break;
                    }
                }
            }
            return tabVector(x_, p_, cols.size());
        }

        template <typename T>
        tabVector operator()(Rcpp::IntegerVector rows, T col_) {
            tabVector m = col((int)col_);
            return m(rows);
        }

        tabMatrix operator()(Rcpp::IntegerVector rows, Rcpp::IntegerVector cols) {
            Rcpp::IntegerVector x_, p_(1);
            for (int col_ : cols) {
                tabVector m = col(col_);
                m(rows);
                x_ = join(x_, m.x);
                p_.push_back(x_.size());
            }
            Rcpp::IntegerVector Dim_(2);
            Dim_[0] = rows.size();
            Dim_[1] = cols.size();
            return tabMatrix(x_, p_, Dim_);
        }

        // COERCION METHODS
        template <>
        Rcpp::S4 as<Rcpp::S4>() {
            Rcpp::S4 s(std::string("tabMatrix"));
            s.slot("x") = x;
            s.slot("p") = p;
            s.slot("Dim") = Dim;
            return s;
        }

        template <>
        RcppSparse::Matrix as<RcppSparse::Matrix>() {
            // convert to dgCMatrix
            std::vector<std::vector<int>> x_(Dim[1]);
            std::vector<std::vector<int>> i_(Dim[1]);

            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int i = 0; i < Dim[1]; ++i) {
                std::vector<int> rows, values;
                for (InnerIterator it(*this, i); it; ++it) {
                    rows.push_back(it.row());
                    values.push_back(it.value());
                }

                // get sort_index of row vector
                std::vector<int> idx(rows.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::stable_sort(idx.begin(), idx.end(),
                                 [&rows](int i1, int i2) { return rows[i1] < rows[i2]; });

                // reorder vals and counts based on sort index
                rows = reorder(idx, rows);
                values = reorder(idx, values);
                x_[i] = values;
                i_[i] = rows;
            }
            // calculate "p" vector
            Rcpp::IntegerVector p_(Dim[1] + 1), i__(p[Dim[1]]);
            Rcpp::NumericVector x__(p[Dim[1]]);
            for (size_t i = 0; i < x_.size(); ++i)
                p_[i + 1] = x_[i].size() + p_[i];
            for (int i = 0, idx = 0; i < Dim[1]; ++i) {
                for (size_t j = 0; j < x_[i].size(); ++j, ++idx) {
                    x__[idx] = x_[i][j];
                    i__[idx] = i_[i][j];
                }
            }
            return Matrix(x__, i__, p_, Dim);
        }

        template <>
        Rcpp::IntegerMatrix as<Rcpp::IntegerMatrix>() {
            Rcpp::IntegerMatrix v(Dim[0], Dim[1]);
            for (int col = 0; col < Dim[1]; ++col) {
                int curr_value;
                for (int ind = p[col]; ind < p[col + 1]; ++ind) {
                    if (x[ind] < 0) {
                        curr_value = std::abs(x[ind]);
                        ++ind;
                    }
                    v(x[ind], col) = curr_value;
                }
            }
            return v;
        }

        // cross-product against dense vector
        Rcpp::NumericVector cross(Rcpp::NumericVector y) {
            Rcpp::NumericVector result(Dim[1]);
            for (int i = 0; i < Dim[1]; ++i)
                for (InnerIterator it(*this, i); it; ++it)
                    result(i) += it.value() * y(it.row());
            return result;
        }

        // cross-product against dense vector
        tabVector cross(Rcpp::IntegerVector y) {
            Rcpp::IntegerVector result(Dim[1]);
            for (int i = 0; i < Dim[1]; ++i)
                for (InnerIterator it(*this, i); it; ++it)
                    result(i) += it.value() * y(it.row());
            return tabVector(y, Dim[0]);
        }

        // cross-product against dense matrix
        Rcpp::NumericMatrix cross(Rcpp::NumericMatrix y) {
            Rcpp::NumericMatrix result(Dim[1], y.cols());
            for (int j = 0; j < y.cols(); ++j) {
                for (int i = 0; i < Dim[1]; ++i)
                    for (InnerIterator it(*this, i); it; ++it)
                        result(i, j) += it.value() * y(it.row(), j);
            }
            return result;
        }

        // cross-product against dense matrix
        Rcpp::IntegerMatrix cross(Rcpp::IntegerMatrix y) {
            Rcpp::IntegerMatrix result(Dim[1], y.cols());
            for (int j = 0; j < y.cols(); ++j) {
                for (int i = 0; i < Dim[1]; ++i)
                    for (InnerIterator it(*this, i); it; ++it)
                        result(i, j) += it.value() * y(it.row(), j);
            }
            return tabMatrix(y);
        }

        Rcpp::IntegerVector colSums() {
            Rcpp::IntegerVector result(Dim[1]);
            for (int col = 0; col < Dim[1]; ++col) {
                if (p[col] != p[col + 1]) {
                    int curr_val = x[p[col]], count = 1;
                    for (int ind = p[col] + 2; ind < p[col + 1]; ++ind) {
                        if (x[ind] < 1) {
                            result(col) += count * curr_val;
                            curr_val = x[ind];
                            ++ind;
                            count = 1;
                        } else ++count;
                    }
                }
            }
            return result;
        }

        int min() {
            return tabVector(x, Dim[0] * Dim[1]).min();
        }

        int max() {
            return tabVector(x, Dim[0] * Dim[1]).max();
        }

        int n_nonzeros() {
            int v = 0;
            for (int val : x)
                if (val >= 0) ++v;
            return v;
        }

        // ------------- OPERATIONS
        // only support operations that maintain sparsity

        // "*"
        tabMatrix operator*(unsigned int y) {}
        tabMatrix operator*(Rcpp::IntegerVector y) {}
        tabMatrix operator*(Rcpp::IntegerMatrix y) {}
        tabMatrix operator*(tabVector y) {}
        tabMatrix operator*(tabMatrix y) {}
        tabMatrix operator*=(unsigned int y) {}
        tabMatrix operator*=(Rcpp::IntegerVector y) {}
        tabMatrix operator*=(Rcpp::IntegerMatrix y) {}
        tabMatrix operator*=(tabVector y) {}
        tabMatrix operator*=(tabMatrix y) {}

    };

}  // namespace RcppSparse

namespace Rcpp {
    namespace traits {

        template <>
        class Exporter<RcppSparse::Matrix> {
            Rcpp::NumericVector x_;
            Rcpp::IntegerVector i, p, Dim;

        public:
            Exporter(SEXP x) {
                Rcpp::S4 s(x);
                if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("i") || !s.hasSlot("Dim"))
                    throw std::invalid_argument("Cannot construct RcppSparse::Matrix from this S4 object");
                x_ = s.slot("x");
                i = s.slot("i");
                p = s.slot("p");
                Dim = s.slot("Dim");
            }

            RcppSparse::Matrix get() {
                return RcppSparse::Matrix(x_, i, p, Dim);
            }
        };

        template <>
        class Exporter<RcppSparse::tabVector> {
            Rcpp::IntegerVector x_;
            int length;

        public:
            Exporter(SEXP x) {
                Rcpp::S4 s(x);
                if (!s.hasSlot("x") || !s.hasSlot("length"))
                    throw std::invalid_argument("Cannot construct RcppSparse::tabVector from this S4 object");
                x_ = s.slot("x");
                length = s.slot("length");
            }

            RcppSparse::tabVector get() {
                return RcppSparse::tabVector(x_, length);
            }
        };

        template <>
        class Exporter<RcppSparse::tabMatrix> {
            Rcpp::IntegerVector Dim, p, x_;

        public:
            Exporter(SEXP x) {
                Rcpp::S4 s(x);
                if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("Dim"))
                    throw std::invalid_argument("Cannot construct RcppSparse::tabMatrix from this S4 object");
                x_ = s.slot("x");
                p = s.slot("p");
                Dim = s.slot("Dim");
            }

            RcppSparse::tabMatrix get() {
                return RcppSparse::tabMatrix(x_, p, Dim);
            }
        };
    }  // namespace traits
}  // namespace Rcpp



#endif
