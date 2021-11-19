#ifndef TABMATRIX_TABVECTOR_H
#define TABMATRIX_TABVECTOR_H

#ifndef TABMATRIX_CSCMATRIX_H
#include "cscMatrix.h"
#endif

namespace tabMatrix {

// intended typenames are Rcpp::IntegerVector or Rcpp::NumericVector
template <int RTYPE>
class tabVector {
   public:
    typedef Rcpp::Vector<RTYPE> vec_t;
    typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;

    // PUBLIC MEMBERS
    vec_t x;
    int length = 0;

    // CONSTRUCTORS
    tabVector(){};

    // construct directly
    template <typename int_type>
    tabVector(vec_t x, int_type length) : x(x), length(length){};

    // interchange between integer and double tabVector types
    tabVector(tabVector<INTSXP> vec) : x(vec.x), length(vec.length){};
    tabVector(tabVector<REALSXP> vec) : x(vec.x), length(vec.length){};

    // construct from sparseVector
    template <typename int_type>
    tabVector(vec_t Ax, Rcpp::IntegerVector Ai, int_type length) : length(length) {
        // TO DO:  catch for zeros in the vector. Possibly implement separate method for converting dense vectors to sparse vectors
        if (Ax.size() > 0) {
            // table up A.x
            Rcpp::IntegerVector counts;
            vec_t vals;
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
            reorder(idx, vals);
            reorder(idx, counts);

            int num_vals = std::accumulate(counts.begin(), counts.end(), 0);
            num_vals += vals.size();

            // return row indices in "i" ordered by values
            x = vec_t(num_vals);
            int ind = 0;
            for (auto value : vals) {
                x[ind] = -std::abs(value);
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

    // BASIC MEMBER FUNCTIONS
    tabVector clone() { return tabVector(Rcpp::clone(x), length); }
    int size() { return x.size(); }
    int length() { return length; }
    int rows() { return length; }

    // ITERATORS

    template <typename int_type>
    tabVector(int_type length) : length(length){};

    // const iterator
    class Iterator {
       public:
        Iterator(tabVector& ptr) : ptr(ptr) {
            max_index_ = ptr.x.size();
            if (index_ < max_index_) {
                value_ = std::abs(ptr.x[index_]);
                ++index_;
            }
        }

        const storage_t value() const { return value_; }
        operator bool() const { return (index_ < max_index_); }
        const int row() const { return ptr.x[index_]; }
        const int index() const { return index_; }

        Iterator& operator++() {
            ++index_;
            if (ptr.x[index_] < 0) {
                value_ = std::abs(ptr.x[index_]);
                ++index_;
            }
            return *this;
        }

        Iterator& operator--() {
            --index_;
            if (index_ > 1) {
                if (ptr.x[index_] < 0) {
                    --index_;
                    int index2 = index_;
                    while (--index2 && index2 >= 0) {
                        if (ptr.x[index2] < 0) {
                            value_ = std::abs(ptr.x[index2]);
                            return *this;
                        }
                    }
                }
            } else {
                index_ = max_index_;
            }
            return *this;
        }

        void nextValue() {
            ++index_;
            while (index_ < max_index_) {
                if (ptr.x[index_] < 0) {
                    value_ = std::abs(ptr.x[index_]);
                    ++index_;
                }
                ++index_;
            }
        }

        void restartValue() {
            // backtrack index to the first instance of the current value
            while (index_ > 0 && ptr.x[index_ - 1] != -value_)
                --index_;
        }

        void prevValue() {
            if (index_ < 2) {
                // do not call prevValue if there is no possibility of a previous value
                index_ = max_index_;
            } else {
                // find index of current value definition
                do --index_;
                while (ptr.x[index_] > 0);
                if (index_ < 2) {
                    // oh drat, there was no previous value
                    index_ = max_index_;
                } else {
                    // find first index of preceding value
                    do --index_;
                    while (ptr.x[index_ - 1] > 0);
                    value_ = std::abs(ptr.x[index_ - 1]);
                }
            }
        }

        void begin() {
            index_ = 1;
            if (1 < max_index_)
                value_ = std::abs(ptr.x[0]);
        }

        void end() {
            index_ = max_index_;
            if (max_index_ > 0) {
                int index2 = index_ - 1;
                while (--index2) {
                    if (ptr.x[index2] < 0) {
                        value_ = std::abs(ptr.x[index2]);
                    }
                }
            }
        }

       private:
        tabVector& ptr;
        storage_t value_;
        int index_ = 0, max_index_ = 0;
    };

    class IntersectionIterator {
       public:
        IntersectionIterator(tabVector<RTYPE>& v1, tabVector<RTYPE>& v2) : v1(v1), v2(v2) {
            if (v1.length != v2.length)
                Rcpp::stop("cannot construct IntersectionIterator from tabVectors of different lengths");
            i1_max = v1.size();
            i2_max = v2.size();
            if (i1_max > 1 && i2_max > 1) {
                x1 = v1.x[0];
                x2 = v2.x[0];
                ++(*this);
            } else {
                is_valid = false;
            }
        }

        operator bool() { return is_valid; }
        const int row() { return v1.x[i1]; }
        const storage_t value1() { return std::abs(x1); }
        const storage_t value2() { return std::abs(x2); }

        IntersectionIterator& operator++() {
            // starting from a value where the index is shared
            ++i1;
            ++i2;
            if (i1 >= i1_max || i2 >= i2_max || v1.x[i1] < 0 || v2.x[i2] < 0) {
                // can we advance v2 to the next value?
                while (i2 < i2_max) {
                    if (v2.x[i2] < 0) break;
                    ++i2;
                }
                if (i2 < i2_max) {
                    x2 = v2.x[i2];
                    ++i2;
                    while (v1.x[i1 - 1] >= 0) --i1;
                } else {
                    // can we advance v1 to the next value?
                    while (i1 < i1_max) {
                        if (v1.x[i1] < 0) break;
                        ++i1;
                    }
                    if (i1 < i1_max) {
                        x1 = v1.x[i1];
                        ++i1;
                        i2 = 1;
                        x2 = v2.x[0];
                    } else {
                        is_valid = false;
                        return *this;  // we're done!
                    }
                }
            }
            while (v1.x[i1] != v2.x[i2]) {
                if (i1 >= i1_max || i2 >= i2_max || v1.x[i1] < 0 || v2.x[i2] < 0) {
                    // can we advance v2 to the next value?
                    while (i2 < i2_max) {
                        if (v2.x[i2] < 0) break;
                        ++i2;
                    }
                    if (i2 < i2_max) {
                        x2 = v2.x[i2];
                        ++i2;
                        while (v1.x[i1 - 1] >= 0) --i1;
                    } else {
                        // can we advance v1 to the next value?
                        while (i1 < i1_max) {
                            if (v1.x[i1] < 0) break;
                            ++i1;
                        }
                        if (i1 < i1_max) {
                            x1 = v1.x[i1];
                            ++i1;
                            i2 = 1;
                            x2 = v2.x[0];
                        } else {
                            is_valid = false;
                            return *this;  // we're done!
                        }
                    }
                } else if (v1.x[i1] < v2.x[i2]) {
                    ++i1;
                } else {
                    ++i2;
                }
            }
            return *this;
        }

       private:
        tabVector &v1, &v2;
        storage_t x1, x2;
        int i1 = 0, i2 = 0, i1_max, i2_max;
        bool is_valid = true;
    };

    // COERCION METHODS

    // as dense vector
    vec_t as_Vector() {
        vec_t v(length);
        for (Iterator it(*this); it; ++it)
            v[it.row()] = it.value();
        return v;
    }

    // as sparse vector
    // see section for explicit specialization for either isparseVector or dsparseVector
    Rcpp::S4 as_sparseVector() {}

    // SUBSETTING

    // scalar
    template <typename int_type>
    storage_t operator()(int_type ind) {
        static_assert(std::is_scalar<int_type>::value), "Cannot subset a tabVector using an object of this type/class");
        for (Iterator it(*this); it; ++it) {
            if (it.row() == ind)
                return it.value();
        }
        return 0;
    }

    // range
    template <typename int_type>
    tabVector operator()(int_type start, int_type stop) {
        static_assert(std::is_scalar<int_type>::value), "Cannot subset a tabVector using an object of this type/class");
        vec_t x_;
        if (length > 0) {
            int curr_val = x[0];
            bool pushed_value = true;
            for (int i = 1; i < x.size(); ++i) {
                if (x[i] < 0) {
                    curr_val = std::abs(x[i]);
                    pushed_value = false;
                    ++i;
                }
                if (x[i] >= start && x[i] <= stop) {
                    if (!pushed_value) {
                        x_.push_back(x[i]);
                        pushed_value = true;
                    }
                    x_.push_back(x[i] - start);
                }
            }
        }
        return tabVector<RTYPE>(x_, stop - start + 1);
    }

    // subset operator, send out to various specializations
    // check if length is 1, if it is a range, otherwise basic subsetting

    // cannot assume that the vector is sorted
    // check if this is a range

    // dense non-negative integer vector
    tabVector operator()(Rcpp::IntegerVector ind) {
        // is this a single index?
        if (ind.size() == 1) {
            return *this(ind[0]);
        }

        // is this a range?
        bool is_range = true;
        for (int i = 0; i < ind.size() - 1; ++i) {
            if (ind[i] + 1 != ind[i + 1]) {
                is_range = false;
                break;
            }
        }
        if (is_range)
            return *this(ind[0], ind[ind.size() - 1]);

        // this is not a range, brute-force subsetting
        vec_t x_;
        Rcpp::IntegerVector i_;
        for (int idx = 0; idx < ind.size(); ++idx) {
            for (Iterator it(*this); it; ++it) {
                if (it.row() == ind[idx]) {
                    x_.push_back(it.value());
                    i_.push_back(idx);
                    break;
                }
            }
        }
        return tabVector<RTYPE>(x_, i_, ind.size());
    }

    // dense logical vector
    tabVector operator()(Rcpp::LogicalVector ind) {
        Rcpp::IntegerVector v = find_true(ind);
        return *this(v);
    }

    // WRITING

    // erase a single index
    template <typename int_type>
    void erase(int_type ind) {
        static_assert(std::is_scalar<int_type>::value), "Cannot subset a tabVector using an object of this type/class");
        if (ind < length) {
            for (int i = 0; i < x.size(); ++i) {
                if (x[i] == ind) {
                    x.erase(ind);
                    if (i == x.size() || (x[i - 1] < 0 && x[i] < 0))
                        x.erase(i - 1);
                }
            }
            length -= 1;
        } else
            Rcpp::stop("cannot erase an index >= length of tabVector")
    }

    // range-based erasure where we ignore values in the range [start, stop]
    template <typename int_type>
    void erase(int_type start, int_type stop) {
        static_assert(std::is_scalar<int_type>::value), "Cannot subset a tabVector using an object of this type/class");
        vec_t x_;
        int diff = stop - start + 1;
        if (length > 0) {
            int curr_val = x[0];
            bool pushed_value = true;
            for (int i = 1; i < x.size(); ++i) {
                if (x[i] < 0) {
                    curr_val = std::abs(x[i]);
                    pushed_value = false;
                    ++i;
                }
                if (x[i] < start) {
                    if (!pushed_value) {
                        x_.push_back(x[i]);
                        pushed_value = true;
                    }
                    x_.push_back(x[i]);
                } else if (x[i] > stop) {
                    if (!pushed_value) {
                        x_.push_back(x[i]);
                        pushed_value = true;
                    }
                    x_.push_back(x[i] - diff);
                }
            }
        }
        x = x_;
        length -= diff;
    }

    // erase a bunch of non-contiguous indices, which is a horrible idea
    void erase(Rcpp::IntegerVector ind) {
        for (int idx : ind)
            erase(idx);
    }

    // ASSIGNMENT

    // assign a single value
    template <typename int_type, typename T>
    void assign(int_type idx, T value) {
        // see if the row already exists, if so, erase it
        for (int i = 0; i < x.size(); ++i) {
            if (x[i] == idx) {
                x.erase(i);
                ++x.length;
                break;
            }
        }
        // now see if the value exists
        int val_ind = x.size();
        for (int i = 0; i < x.size(); ++i) {
            if (-x[i] == (storage_t)value) {
                val_ind = i + 1;
                break;
            }
        }

        // if value exists
        if (val_ind < x.size()) {
            while (x[val_ind] < idx) ++val_ind;
            // insert value at val_ind
            x.insert(idx, val_ind);
        } else {
            x.push_back((storage_t)-value);
            x.push_back(idx);
        }
    }

    // assign multiple values, which is a horrible idea. Should probably construct using sparseVector or other methods.
    void assign(Rcpp::IntegerVector idx, Rcpp::Vector<RTYPE> value) {
        for (int i = 0; i < idx.size(); ++i)
            assign(idx[i], value[i]);
    }

    // ARITHMETIC
    Rcpp::LogicalVector is_zero() {
        Rcpp::LogicalVector l(length, false);
        for (Iterator it(*this); it; ++it)
            l[it.row()] = true;
        return l;
    }

    // cumulative functions
    storage_t sum() {
        storage_t v = 0;
        for (Iterator it(*this); it; ++it)
            v += it.value();
        return v;
    }

    storage_t max() {
        storage_t max_val = x[0];
        for (storage_t val : x) {
            if (val < max_val)
                max_val = val;
        }
        return std::abs(max_val);
    }

    storage_t min() {
        storage_t min_val = x[0];
        for (storage_t val : x) {
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
    // arithmetic against scalars or dense vectors gives dense outputs

    // +-/ with scalars or dense vectors returns dense
    vec_t operator+(storage_t y) { return as_Vector() + y; }
    vec_t operator-(storage_t y) { return as_Vector() - y; }
    vec_t operator/(storage_t y) { return as_Vector() / y; }
    vec_t operator*(storage_t y) { return as_Vector() / y; }

    vec_t operator+(vec_t y) {
        vec_t v = y.clone();
        for (Iterator it(*this); it; ++it)
            v[it.row()] += it.value();
        return v;
    }
    vec_t operator-(vec_t y) {
        vec_t v = y.clone();
        for (Iterator it(*this); it; ++it)
            v[it.row()] -= it.value();
        return v;
    }
    vec_t operator/(vec_t y) {
        vec_t v = y.clone();
        for (Iterator it(*this); it; ++it)
            v[it.row()] /= it.value();
        return v;
    }
    // multiply operator for vectors
    vec_t operator*(vec_t y) {
        vec_t v = y.clone();
        for (Iterator it(*this); it; ++it)
            v[it.row()] *= it.value();
        return v;
    }

    // tabVectors may be multiplied against one another, but division by zero, addition of negative values, or subtraction could cause undefined behavior or introduce zeros

    // multiply in place operator for scalar
    void operator*=(storage_t y) {
        if (y == 0) return tabVector<RTYPE>();
        y = std::abs(y);
        for (int& x_ : x)
            if (x_ < 0)
                x_ *= y;
    }
    tabVector operator*(tabVector<RTYPE> y) {
        vec_t x_;
        Rcpp::IntegerVector i_;
        for (IntersectionIterator it(*this, y); it; ++it) {
            x_.push_back(it.value1() * it.value2());
            i_.push_back(it.row());
        }
        return tabVector<RTYPE>(x_, i_, y.length);
    }
    tabVector operator*=(tabVector<RTYPE> y) {
        *this = (*this * y);
    }

    /*
    addition of two tabVectors is possible, but it requires a UnionIterator of sorts

    tabVector operator+(tabVector<RTYPE> y){
        // TO DO
    }
    tabVector operator+=(tabVector<RTYPE> y){
        // TO DO
    }
    */

    // crossproduct
    storage_t crossprod(vec_t y) {
        storage_t value = 0;
        for (Iterator it(*this); it; ++it)
            value += it.value() * y(it.row());
        return value;
    }

    storage_t crossprod(tabVector<RTYPE> y) {
        storage_t value = 0;
        for (IntersectionIterator it(*this, y); it; ++it)
            value += it.value1() * it.value2();
        return value;
    }

    storage_t crossprod() {
        storage_t value = 0;
        for (Iterator it(*this); it; ++it)
            value += it.value() * it.value();
        return value;
    }

    // tcrossprod
    Rcpp::Matrix<RTYPE> tcrossprod() {
        Rcpp::Matrix<RTYPE> res(length, length);
        for (Iterator it(*this); it; ++it)
            for (Iterator it2(*this); it2; ++it2)
                res(it.row(), it2.row()) = it.value() * it2.value();
        return res;
    }
    Rcpp::Matrix<RTYPE> tcrossprod(tabVector<RTYPE> y) {
        if (length != y.length) Rcpp::stop("length of parent object and 'y' for 'tcrossprod' operation are not compatible");
        Rcpp::Matrix<RTYPE> res(length, y.length);
        for (Iterator it(*this); it; ++it)
            for (Iterator it2(*this); it2; ++it2)
                res(it.row(), it2.row()) = it.value() * it2.value();
        return res;
    }
    Rcpp::Matrix<RTYPE> tcrossprod(vec_t y) {
        if (length != y.size()) Rcpp::stop("length of parent object and 'y' for 'tcrossprod' operation are not compatible");
        Rcpp::Matrix<RTYPE> res(length, y.size());
        for (Iterator it(*this); it; ++it)
            for (int i = 0; i < y.size(); ++i)
                res(it.row(), i) = it.value() * y[i];
        return res;
    }

    void square() {
        for (auto& x_ : x)
            if (x_ < 0)
                x_ *= -x_;
    }

    // TO DO: dist could really benefit from a UnionIterator
    storage_t dist(tabVector<RTYPE> y, std::string method = "euclidean") {
        if (method == "euclidean" || method == "square") {
            storage_t v = y.square().sum() + x.square().sum();
            for (IntersectionIterator it(*this, y); it; ++it) {
                v -= it.value1() * it.value1() + it.value2() * it.value2();
                v += std::pow(it.value1() - it.value2(), 2);
            }
            if (method == "square") return v;
            return std::sqrt(v);
        } else if (method == "cosine") {
            return std::sqrt(crossprod(x)) / (square().sum() * y.square().sum());
        } else if (method == "jaccard") {
            storage_t intersection = crossprod(y);
            return intersection / (sum() + y.sum() - intersection);
        } else if (method == "abs") {
            storage_t v = y.sum() + x.sum();
            for (IntersectionIterator it(*this, y); it; ++it) {
                v -= it.value1() + it.value2();
                v += std::abs(it.value1() - it.value2());
            }
        }
        Rcpp::stop("distance method type is not currently implemented.");
    }

    // LOGICAL COMPARISONS
    // scalar comparizons against either "int" or "double"
    Rcpp::IntegerVector operator==(storage_t y) {
        Rcpp::IntegerVector result(length);
        int idx = 0;
        if (y == 0) {
            for (Iterator it(*this); it; ++it) {
                result(idx) = it.row();
                ++idx;
            }
        } else {
            for (Iterator it(*this); it; ++it) {
                if (it.value() == y) {
                    result(idx) = it.row();
                    ++idx;
                }
            }
        }
        return result[Rcpp::Range(0, idx)];
    }

    // TO DO
    // == != > < <= >= for tabVector vs. tabVector
    // all other comparisons vs storage_t
    // all other comparisons vs. vec_t

    int n_nonzeros() {
        int v = 0;
        for (int i = 0; i < x.size(); ++i)
            if (x[i] < 0)
                ++v;
        return v;
    }

    // member functions that require explicit specialization
    // specializations needed either to indicate preservation of storage type (i.e. int vs. int operations)
    //    or to indicate conversion of INTSXP tabVector to REALSXP tabVector

    Rcpp::S4 as_S4() {}

    void pow(storage_t y) {
        for (auto& x_ : x)
            if (x_ < 0)
                x_ = std::abs(x_) ^ y;
    }
};

// EXPLICIT SPECIALIZATIONS
// member functions that must be specialized by template (i.e. INTSXP or REALSXP)

// Rcpp::as to S4 objects
template <>
Rcpp::S4 tabVector<INTSXP>::as_S4() {
    Rcpp::S4 s(std::string("itabVector"));
    s.slot("x") = x;
    s.slot("length") = length;
    return s;
}

template <>
Rcpp::S4 tabVector<REALSXP>::as_S4() {
    Rcpp::S4 s(std::string("dtabVector"));
    s.slot("x") = x;
    s.slot("length") = length;
    return s;
}

// Coerce to sparseVector
template <>
Rcpp::S4 tabVector<INTSXP>::as_sparseVector() {
    Rcpp::S4 s(std::string("isparseVector"));
    Rcpp::IntegerVector x_, i_;
    for (Iterator it(*this); it; ++it) {
        x_.push_back(it.value);
        i_.push_back(it.row());
    }
    s.slot("x") = x_;
    s.slot("i") = i_;
    s.slot("length") = length();
    return s;
}

template <>
Rcpp::S4 tabVector<REALSXP>::as_sparseVector() {
    Rcpp::S4 s(std::string("dsparseVector"));
    Rcpp::IntegerVector i_;
    Rcpp::NumericVector x_;
    for (Iterator it(*this); it; ++it) {
        x_.push_back(it.value);
        i_.push_back(it.row());
    }
    s.slot("x") = x_;
    s.slot("i") = i_;
    s.slot("length") = length();
    return s;
}

}  // namespace tabMatrix

// in-place dense arithmetic
template <int RTYPE>
void operator+=(Rcpp::Vector<RTYPE> y, tabMatrix::tabVector<RTYPE> x) {
    for (typename tabMatrix::tabVector<RTYPE>::Iterator it(x); it; ++it)
        y[it.row()] += it.value();
}
template <int RTYPE>
void operator-=(Rcpp::Vector<RTYPE> y, tabMatrix::tabVector<RTYPE> x) {
    for (typename tabMatrix::tabVector<RTYPE>::Iterator it(x); it; ++it)
        y[it.row()] -= it.value();
}
template <int RTYPE>
void operator*=(Rcpp::Vector<RTYPE> y, tabMatrix::tabVector<RTYPE> x) {
    for (typename tabMatrix::tabVector<RTYPE>::Iterator it(x); it; ++it)
        y[it.row()] *= it.value();
}
template <int RTYPE>
void operator/=(Rcpp::Vector<RTYPE> y, tabMatrix::tabVector<RTYPE> x) {
    for (typename tabMatrix::tabVector<RTYPE>::Iterator it(x); it; ++it)
        y[it.row()] /= it.value();
}

// TO DO. Overload any functions or operators where the RHS is a tabVector and the LHS is dense or scalar

// TO DO. Overload any member functions as functional interface where one or both of the LHS or RHS include tabVector

#endif