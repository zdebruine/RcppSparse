#ifndef TABMATRIX_TABMATRIX_H
#define TABMATRIX_TABMATRIX_H

#ifndef TABMATRIX_TABVECTOR_H
#include "tabVector.h"
#endif

namespace tabMatrix {
template <int RTYPE>
class tabMatrix {
   public:
    // public member objects
    Rcpp::IntegerVector x, p, Dim;

    // constructors
    tabMatrix(Rcpp::IntegerVector x, Rcpp::IntegerVector p, Rcpp::IntegerVector Dim) : x(x), p(p), Dim(Dim) {}
    tabMatrix(const Rcpp::S4& s) {
        if (!s.hasSlot("x") || !s.hasSlot("i") || !s.hasSlot("Dim"))
            throw std::invalid_argument("Cannot construct tabMatrix::tabMatrix from this S4 object");
        x = s.slot("x");
        p = s.slot("p");
        Dim = s.slot("Dim");
    }
    tabMatrix(cscMatrix A) {
        std::vector<std::vector<int>> x_(A.cols());
        std::vector<int> p__ = Rcpp::as<std::vector<int>>(A.p);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned int i = 0; i < A.cols(); ++i) {
            int start = A.p[i], stop = A.p[i + 1];
            if (stop > start) {
                // table up A.x
                std::vector<int> counts;
                Rcpp::Vector<RTYPE> vals;
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
                reorder(idx, vals);
                reorder(idx, counts);

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
    tabMatrix(Rcpp::IntegerMatrix A) {
        Dim = Rcpp::IntegerVector(2);
        p = Rcpp::IntegerVector(1);
        Dim[0] = A.rows();
        Dim[1] = A.cols();
        for (int i = 0; i < A.cols(); ++i) {
            Rcpp::IntegerVector A_col_i = A.column(i);
            tabVector<RTYPE> x_ = tabVector<RTYPE>(A_col_i, Dim[0]);
            x = join(x, x_.x);
            p.push_back(x_.size());
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

    tabVector<RTYPE> col(int col) {
        return tabVector<RTYPE>(x[Rcpp::Range(p[col], p[col + 1] - 1)], rows());
    }

    template <typename T>
    tabVector<RTYPE> row(T i) {
        Rcpp::IntegerVector x_, i_;
        for (RowIterator it(*this, (int)i); it; ++it) {
            x_.push_back(it.value());
            i_.push_back(it.col());
        }
        return tabVector<RTYPE>(x_, i_, Dim[0]);
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
                do {
                    --i;
                } while (x[i] > 0);
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
    tabVector<RTYPE> operator()(T row, Rcpp::IntegerVector cols) {
        Rcpp::IntegerVector x_, p_;
        for (int col : cols) {
            // find index of row in col
            for (int j = p[col]; j < p[col + 1]; ++j) {
                if (x[j] == (int)row) {
                    p_.push_back(col);
                    // backtrack to last value
                    do {
                        --j;
                    } while (x[j] > 0);
                    x_.push_back(std::abs(x[j]));
                    break;
                }
            }
        }
        return tabVector<RTYPE>(x_, p_, cols.size());
    }

    template <typename T>
    tabVector<RTYPE> operator()(Rcpp::IntegerVector rows, T col_) {
        tabVector<RTYPE> m = col((int)col_);
        return m(rows);
    }

    tabMatrix operator()(Rcpp::IntegerVector rows, Rcpp::IntegerVector cols) {
        Rcpp::IntegerVector x_, p_(1);
        for (int col_ : cols) {
            tabVector<RTYPE> m = col(col_);
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
    Rcpp::S4 as_S4() {
        Rcpp::S4 s(std::string("tabMatrix"));
        s.slot("x") = x;
        s.slot("p") = p;
        s.slot("Dim") = Dim;
        return s;
    }

    cscMatrix as_cscMatrix() {
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
            reorder(idx, rows);
            reorder(idx, values);
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
        return cscMatrix(x__, i__, p_, Dim);
    }

    Rcpp::IntegerMatrix as_IntegerMatrix() {
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
    tabVector<RTYPE> cross(Rcpp::IntegerVector y) {
        Rcpp::IntegerVector result(Dim[1]);
        for (int i = 0; i < Dim[1]; ++i)
            for (InnerIterator it(*this, i); it; ++it)
                result(i) += it.value() * y(it.row());
        return tabVector<RTYPE>(y, Dim[0]);
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

    // cross-product against dense integer matrix
    tabMatrix cross(Rcpp::IntegerMatrix y) {
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
                    } else
                        ++count;
                }
            }
        }
        return result;
    }

    int min() {
        return tabVector<RTYPE>(x, Dim[0] * Dim[1]).min();
    }

    int max() {
        return tabVector<RTYPE>(x, Dim[0] * Dim[1]).max();
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
    tabMatrix operator*(tabVector<RTYPE> y) {}
    tabMatrix operator*(tabMatrix y) {}
    tabMatrix operator*=(unsigned int y) {}
    tabMatrix operator*=(Rcpp::IntegerVector y) {}
    tabMatrix operator*=(Rcpp::IntegerMatrix y) {}
    tabMatrix operator*=(tabVector<RTYPE> y) {}
    tabMatrix operator*=(tabMatrix y) {}
};

}  // namespace tabMatrix

#endif