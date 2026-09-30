// Copyright (c) 2026 Teodor Iličin and Rok Žitko. BSD-3-Clause.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace py = pybind11;
using Word = std::uint64_t;
using Complex = std::complex<double>;
using Index = std::int64_t;

static int popcount(Word x) {
    int count = 0;
    while (x) { x &= x - 1; ++count; }
    return count;
}

class Basis {
public:
    int nimp, nmodes, nwords, cutoff, parity, target_sz, target_eta, target_number;
    std::vector<int> spins, eta;
    std::vector<Word> states;
    std::vector<int> qp_counts;
    std::vector<Index> table;
    std::size_t max_dimension;
    static constexpr int no_sz = 1000000000;

    Basis(int nimp_, const std::vector<int>& spins_, int cutoff_, int parity_,
          int target_sz_, const std::vector<int>& eta_, int target_eta_,
          std::size_t max_dimension_, int target_number_ = -1)
        : nimp(nimp_), nmodes(static_cast<int>(spins_.size())),
          nwords((nmodes + 63) / 64), cutoff(cutoff_), parity(parity_),
          target_sz(target_sz_), target_eta(target_eta_), target_number(target_number_), spins(spins_), eta(eta_),
          max_dimension(max_dimension_) {
        if (nimp < 0 || nimp > nmodes || nmodes == 0 || cutoff < 0 ||
            cutoff > nmodes - nimp || (parity != 0 && parity != 1) ||
            max_dimension == 0 || target_number < -1 || target_number > nmodes)
            throw std::invalid_argument("invalid Fock-basis dimensions or sector");
        for (int s : spins) if (s != 1 && s != -1)
            throw std::invalid_argument("mode spin labels must be +1 or -1");
        if (target_eta != 0 && target_eta != 1 && target_eta != -1)
            throw std::invalid_argument("eta must be 0, +1 or -1");
        if (target_eta && eta.size() != spins.size())
            throw std::invalid_argument("missing eta labels");
        for (int e : eta) if (e != 1 && e != -1)
            throw std::invalid_argument("mode eta labels must be +1 or -1");
        std::vector<std::vector<int>> groups(4);
        for (int i = 0; i < nmodes; ++i)
            groups[(i < nimp ? 0 : 2) + (spins[i] == -1)].push_back(i);
        std::vector<Word> state(nwords, 0);
        std::vector<int> counts(4);
        std::function<void(int)> choose_group;
        std::function<void(int,int,int)> choose;
        int current_qp = 0;
        choose_group = [&](int g) {
            if (g == 4) {
                if (target_eta) {
                    int e = 1;
                    for (int i = 0; i < nmodes; ++i)
                        if ((state[i / 64] >> (i % 64)) & 1) e *= eta[i];
                    if (e != target_eta) return;
                }
                if (qp_counts.size() >= max_dimension)
                    throw std::length_error("basis exceeds max_dimension");
                states.insert(states.end(), state.begin(), state.end());
                qp_counts.push_back(current_qp);
            } else choose(g, 0, counts[g]);
        };
        choose = [&](int g, int start, int left) {
            if (!left) { choose_group(g + 1); return; }
            const auto& modes = groups[g];
            for (int j = start; j <= static_cast<int>(modes.size()) - left; ++j) {
                int i = modes[j];
                state[i / 64] |= Word(1) << (i % 64);
                choose(g, j + 1, left - 1);
                state[i / 64] &= ~(Word(1) << (i % 64));
            }
        };
        for (int q = 0; q <= cutoff; ++q) {
            current_qp = q;
            for (int iu = 0; iu <= static_cast<int>(groups[0].size()); ++iu)
                for (int id = 0; id <= static_cast<int>(groups[1].size()); ++id) {
                    if ((iu + id + q) % 2 != parity) continue;
                    if (target_number >= 0 && iu + id + q != target_number) continue;
                    for (int bu = 0; bu <= q; ++bu) {
                        int bd = q - bu;
                        if (bu > static_cast<int>(groups[2].size()) ||
                            bd > static_cast<int>(groups[3].size())) continue;
                        if (target_sz != no_sz && iu - id + bu - bd != target_sz) continue;
                        counts = {iu, id, bu, bd};
                        choose_group(0);
                    }
                }
        }
        if (qp_counts.empty()) throw std::invalid_argument("empty symmetry sector");
        std::size_t capacity = 1;
        while (capacity < 2 * size()) {
            if (capacity > std::numeric_limits<std::size_t>::max() / 2)
                throw std::length_error("basis index overflow");
            capacity *= 2;
        }
        table.assign(capacity, -1);
        for (std::size_t i = 0; i < size(); ++i) {
            const Word* s = data(i);
            auto pos = hash(s) & (table.size() - 1);
            while (table[pos] != -1) pos = (pos + 1) & (table.size() - 1);
            table[pos] = static_cast<Index>(i);
        }
    }

    std::size_t size() const { return qp_counts.size(); }
    const Word* data(std::size_t i) const { return states.data() + i * nwords; }
    std::size_t hash(const Word* s) const {
        Word h = 0x9e3779b97f4a7c15ULL;
        for (int w = 0; w < nwords; ++w) {
            Word x = s[w] + 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            h ^= (x ^ (x >> 31)) + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
    Index lookup(const Word* s) const {
        auto pos = hash(s) & (table.size() - 1);
        while (table[pos] != -1) {
            auto i = table[pos];
            if (std::equal(s, s + nwords, data(i))) return i;
            pos = (pos + 1) & (table.size() - 1);
        }
        return -1;
    }
    py::array_t<Word> occupations() const {
        py::array_t<Word> a({static_cast<py::ssize_t>(size()), static_cast<py::ssize_t>(nwords)});
        std::copy(states.begin(), states.end(), a.mutable_data());
        return a;
    }
    py::array_t<int> counts() const {
        py::array_t<int> a(size());
        std::copy(qp_counts.begin(), qp_counts.end(), a.mutable_data());
        return a;
    }
    std::size_t bytes() const {
        return states.capacity() * sizeof(Word) + qp_counts.capacity() * sizeof(int)
             + table.capacity() * sizeof(Index);
    }
};

struct Term {
    Complex coefficient;
    std::vector<int> ops;
    int dq = 0;
};

class Operator {
    std::shared_ptr<Basis> basis;
    std::vector<Term> creation;
    std::vector<std::vector<Term>> anchored;
    std::vector<Complex> diagonal;
public:
    bool real = true;
    Operator(std::shared_ptr<Basis> b, const std::vector<Complex>& coefficients,
             const std::vector<std::vector<int>>& strings) : basis(std::move(b)) {
        if (!basis) throw std::invalid_argument("operator requires a non-null basis");
        if (coefficients.size() != strings.size())
            throw std::invalid_argument("operator coefficient/string count mismatch");
        anchored.resize(basis->nmodes);
        std::vector<Term> diagonal_terms;
        Complex constant = 0;
        std::vector<Complex> onebody(basis->nmodes, 0);
        for (std::size_t k = 0; k < strings.size(); ++k) {
            Term term{coefficients[k], strings[k], 0};
            if (!std::isfinite(term.coefficient.real()) || !std::isfinite(term.coefficient.imag()))
                throw std::invalid_argument("non-finite operator coefficient");
            if (term.coefficient.imag() != 0) real = false;
            std::vector<int> cre, ann;
            for (int op : term.ops) {
                // Check signed bounds before abs: abs(INT_MIN) is undefined.
                if (!op || op < -basis->nmodes || op > basis->nmodes)
                    throw std::invalid_argument("fermionic mode outside basis");
                int i = std::abs(op) - 1;
                // The diagonal and anchor optimizations require normal ordering,
                // without repeated creators/annihilators of the same mode.
                if ((op > 0 && (!ann.empty() || (!cre.empty() && i <= cre.back()))) ||
                    (op < 0 && !ann.empty() && i >= ann.back()))
                    throw std::invalid_argument("fermionic strings must be canonical: ascending creators then descending annihilators");
                if (i >= basis->nimp) term.dq += op > 0 ? 1 : -1;
                (op > 0 ? cre : ann).push_back(i);
            }
            std::reverse(ann.begin(), ann.end());
            if (cre == ann) {
                if (cre.empty()) constant += term.coefficient;
                else if (cre.size() == 1) onebody[cre[0]] += term.coefficient;
                else diagonal_terms.push_back(term);
            } else {
                int anchor = -1;
                for (auto it = term.ops.rbegin(); it != term.ops.rend(); ++it)
                    if (*it < 0) { anchor = -*it - 1; break; }
                if (anchor < 0) creation.push_back(term);
                else anchored[anchor].push_back(term);
            }
        }
        diagonal.assign(basis->size(), constant);
        for (std::size_t col = 0; col < basis->size(); ++col) {
            const Word* s = basis->data(col);
            for (int i = 0; i < basis->nmodes; ++i)
                if ((s[i / 64] >> (i % 64)) & 1) diagonal[col] += onebody[i];
            for (const auto& term : diagonal_terms) {
                bool occupied = true;
                for (int op : term.ops) {
                    if (op < 0) break;
                    int i = op - 1;
                    if (!((s[i / 64] >> (i % 64)) & 1)) { occupied = false; break; }
                }
                if (occupied) diagonal[col] += term.coefficient;
            }
        }
    }

    Index apply(const Term& term, const Word* source, std::vector<Word>& target, int& sign,
                const Basis& destination) const {
        std::copy(source, source + basis->nwords, target.begin());
        sign = 1;
        for (auto it = term.ops.rbegin(); it != term.ops.rend(); ++it) {
            const int i = std::abs(*it) - 1, w = i / 64, bit = i % 64;
            const bool occupied = (target[w] >> bit) & 1;
            if (occupied == (*it > 0)) return -1;
            int n = 0;
            for (int a = 0; a < w; ++a) n += popcount(target[a]);
            n += popcount(target[w] & ((Word(1) << bit) - 1));
            if (n % 2) sign = -sign;
            target[w] ^= Word(1) << bit;
        }
        return destination.lookup(target.data());
    }

    template<class Function>
    void connections(std::size_t col, std::vector<Word>& target, Function&& emit,
                     const Basis* other = nullptr) const {
        const Basis& destination = other ? *other : *basis;
        const Word* s = basis->data(col);
        const int q = basis->qp_counts[col];
        auto process = [&](const Term& term) {
            if (q + term.dq < 0 || q + term.dq > destination.cutoff) return;
            int sign;
            Index row = apply(term, s, target, sign, destination);
            if (row >= 0) emit(row, static_cast<double>(sign) * term.coefficient);
        };
        for (const auto& term : creation) process(term);
        for (int w = 0; w < basis->nwords; ++w) {
            Word bits = s[w];
            int bit = 0;
            while (bits) {
                if (bits & 1) for (const auto& term : anchored[64*w + bit]) process(term);
                bits >>= 1; ++bit;
            }
        }
    }

    template<class T>
    py::array_t<T> matvec(py::array_t<T, py::array::c_style | py::array::forcecast> input) const {
        if (input.ndim() != 1 || static_cast<std::size_t>(input.size()) != basis->size())
            throw std::invalid_argument("state-vector dimension mismatch");
        if constexpr (std::is_same_v<T, double>)
            if (!real) throw std::invalid_argument("real action requested for complex operator");
        py::array_t<T> output(basis->size());
        const T* x = input.data(); T* y = output.mutable_data();
        py::gil_scoped_release release;
        std::fill(y, y + basis->size(), T(0));
        std::vector<Word> target(basis->nwords);
        for (std::size_t col = 0; col < basis->size(); ++col) {
            if (x[col] == T(0)) continue;
            if constexpr (std::is_same_v<T, double>) y[col] += diagonal[col].real() * x[col];
            else y[col] += diagonal[col] * x[col];
            connections(col, target, [&](Index row, Complex v) {
                if constexpr (std::is_same_v<T, double>) y[row] += v.real() * x[col];
                else y[row] += v * x[col];
            });
        }
        return output;
    }

    py::array_t<Complex> diag() const {
        py::array_t<Complex> a(diagonal.size());
        std::copy(diagonal.begin(), diagonal.end(), a.mutable_data());
        return a;
    }

    py::array_t<Complex> between(py::array_t<Complex, py::array::c_style | py::array::forcecast> input,
                                std::shared_ptr<Basis> destination) const {
        if (!destination) throw std::invalid_argument("transition requires a non-null destination basis");
        if (basis->nmodes != destination->nmodes || basis->nimp != destination->nimp)
            throw std::invalid_argument("transition bases use different fermionic modes");
        if (input.ndim() != 1 || static_cast<std::size_t>(input.size()) != basis->size())
            throw std::invalid_argument("state-vector dimension mismatch");
        py::array_t<Complex> output(destination->size());
        const Complex* x = input.data(); Complex* y = output.mutable_data();
        py::gil_scoped_release release;
        std::fill(y, y + destination->size(), Complex(0));
        std::vector<Word> target(basis->nwords);
        for (std::size_t col = 0; col < basis->size(); ++col) {
            if (x[col] == Complex(0)) continue;
            Index row = destination->lookup(basis->data(col));
            if (row >= 0) y[row] += diagonal[col] * x[col];
            connections(col, target, [&](Index r, Complex value) { y[r] += value*x[col]; }, destination.get());
        }
        return output;
    }

    py::tuple sparse(std::size_t max_nnz) const {
        std::vector<Index> pointers(basis->size() + 1, 0), rows;
        std::vector<Complex> values;
        {
            py::gil_scoped_release release;
            std::vector<Word> target(basis->nwords);
            std::vector<std::pair<Index,Complex>> column;
            for (std::size_t col = 0; col < basis->size(); ++col) {
                column.clear();
                if (diagonal[col] != Complex(0)) column.emplace_back(col, diagonal[col]);
                connections(col, target, [&](Index row, Complex v) { column.emplace_back(row, v); });
                std::sort(column.begin(), column.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                for (std::size_t i = 0; i < column.size();) {
                    Index row = column[i].first; Complex v = 0;
                    do { v += column[i++].second; } while (i < column.size() && column[i].first == row);
                    if (v == Complex(0)) continue;
                    if (values.size() >= max_nnz) throw std::length_error("sparse matrix exceeds max_nnz; use matrix-free action");
                    rows.push_back(row); values.push_back(v);
                }
                pointers[col + 1] = values.size();
            }
        }
        py::array_t<Complex> v(values.size());
        py::array_t<Index> r(rows.size()), p(pointers.size());
        std::copy(values.begin(), values.end(), v.mutable_data());
        std::copy(rows.begin(), rows.end(), r.mutable_data());
        std::copy(pointers.begin(), pointers.end(), p.mutable_data());
        return py::make_tuple(v, r, p);
    }
};

PYBIND11_MODULE(_core, m) {
    m.doc() = "Fermionic Fock bases and projected operator actions (C++17)";
    py::class_<Basis, std::shared_ptr<Basis>>(m, "Basis")
        .def(py::init<int,const std::vector<int>&,int,int,int,const std::vector<int>&,int,std::size_t>())
        .def(py::init<int,const std::vector<int>&,int,int,int,const std::vector<int>&,int,std::size_t,int>())
        .def_property_readonly("dimension", &Basis::size)
        .def_property_readonly("bytes", &Basis::bytes)
        .def("occupations", &Basis::occupations)
        .def("qp_counts", &Basis::counts);
    py::class_<Operator>(m, "Operator")
        .def(py::init<std::shared_ptr<Basis>,const std::vector<Complex>&,const std::vector<std::vector<int>>&>(),
             "Construct from canonical strings: ascending creators followed by descending annihilators (no repetitions).")
        .def_readonly("real", &Operator::real)
        .def("matvec_real", &Operator::matvec<double>)
        .def("matvec_complex", &Operator::matvec<Complex>)
        .def("between", &Operator::between)
        .def("diagonal", &Operator::diag)
        .def("sparse", &Operator::sparse);
}
