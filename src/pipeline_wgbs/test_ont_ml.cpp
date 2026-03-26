/**
 * test_ont_ml.cpp
 * Unit tests for subset_to_Cm_section ML slicing when C+h and C+m have
 * different numbers of positions.
 *
 * Compile alongside patter object files:
 *   g++ -std=c++11 -o test_ont_ml test_ont_ml.cpp ont.o patter_utils.o
 *
 * Or via setup.py:
 *   python setup.py --targets test_ont_ml
 */
#include "patter.h"
#include <cassert>
#include <iostream>
#include <sstream>

// ---- helpers ---------------------------------------------------------------

static std::string vec2str(const std::vector<int> &v) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < v.size(); i++) {
        if (i) oss << ",";
        oss << v[i];
    }
    oss << "]";
    return oss.str();
}

static void check_eq(const std::string &label,
                     const std::string &got,
                     const std::string &expected) {
    if (got != expected) {
        std::cerr << "FAIL  " << label << "\n"
                  << "      got:      " << got << "\n"
                  << "      expected: " << expected << "\n";
        std::exit(1);
    }
    std::cout << "PASS  " << label << "\n";
}

// ---- tests -----------------------------------------------------------------

/* Test 1: C+h (1 position) comes before C+m (2 positions).
 * ML has 3 values: 1 for C+h, then 2 for C+m.
 *
 * Bug: current code uses MM_pos * nr_MM_vals as offset, so for C+m it
 * checks  ML_vec.size() >= (1+1)*2 = 4  →  3 >= 4  is FALSE,
 * leaving ML_str un-trimmed (all 3 values), which then triggers a
 * size-mismatch throw in parse_np_fields_by_mod.
 */
static void test_unequal_h_before_m() {
    // C+h has 1 MM entry; C+m has 2 MM entries.
    // ML: first 1 value for C+h, next 2 for C+m.
    std::string MM_h = "C+h,1;C+m,0,3";
    std::string ML_h = ",100,200,180";
    bool dot_h = false;
    subset_to_Cm_section(MM_h, ML_h, dot_h, "h");
    check_eq("test_unequal_h_before_m / MM for h", MM_h, "C+h,1");
    check_eq("test_unequal_h_before_m / ML for h", ML_h, ",100");

    std::string MM_m = "C+h,1;C+m,0,3";
    std::string ML_m = ",100,200,180";
    bool dot_m = false;
    subset_to_Cm_section(MM_m, ML_m, dot_m, "m");
    check_eq("test_unequal_h_before_m / MM for m", MM_m, "C+m,0,3");
    check_eq("test_unequal_h_before_m / ML for m", ML_m, ",200,180");
}

/* Test 2: C+m (3 positions) comes before C+h (1 position) — also common. */
static void test_unequal_m_before_h() {
    std::string MM_m = "C+m,0,3,2;C+h,1";
    std::string ML_m = ",200,180,160,220";
    bool dot_m = false;
    subset_to_Cm_section(MM_m, ML_m, dot_m, "m");
    check_eq("test_unequal_m_before_h / MM for m", MM_m, "C+m,0,3,2");
    check_eq("test_unequal_m_before_h / ML for m", ML_m, ",200,180,160");

    std::string MM_h = "C+m,0,3,2;C+h,1";
    std::string ML_h = ",200,180,160,220";
    bool dot_h = false;
    subset_to_Cm_section(MM_h, ML_h, dot_h, "h");
    check_eq("test_unequal_m_before_h / MM for h", MM_h, "C+h,1");
    check_eq("test_unequal_m_before_h / ML for h", ML_h, ",220");
}

/* Test 3: Equal counts — regression, should still work. */
static void test_equal_counts() {
    std::string MM_h = "C+m,0,3;C+h,0,2";
    std::string ML_h = ",200,180,120,100";
    bool dot_h = false;
    subset_to_Cm_section(MM_h, ML_h, dot_h, "h");
    check_eq("test_equal_counts / MM for h", MM_h, "C+h,0,2");
    check_eq("test_equal_counts / ML for h", ML_h, ",120,100");

    std::string MM_m = "C+m,0,3;C+h,0,2";
    std::string ML_m = ",200,180,120,100";
    bool dot_m = false;
    subset_to_Cm_section(MM_m, ML_m, dot_m, "m");
    check_eq("test_equal_counts / MM for m", MM_m, "C+m,0,3");
    check_eq("test_equal_counts / ML for m", ML_m, ",200,180");
}

/* Test 4: Only C+m (no C+h) — regression. */
static void test_only_m() {
    std::string MM_m = "C+m,0,3,2";
    std::string ML_m = ",200,180,160";
    bool dot_m = false;
    subset_to_Cm_section(MM_m, ML_m, dot_m, "m");
    check_eq("test_only_m / MM", MM_m, "C+m,0,3,2");
    check_eq("test_only_m / ML", ML_m, ",200,180,160");

    std::string MM_h = "C+m,0,3,2";
    std::string ML_h = ",200,180,160";
    bool dot_h = false;
    subset_to_Cm_section(MM_h, ML_h, dot_h, "h");
    check_eq("test_only_m / MM for h (should be empty)", MM_h, "");
    check_eq("test_only_m / ML for h (should be empty)", ML_h, "");
}

/* Test 5: Large unequal counts matching the pattern seen in tiny_test.bam
 * (C+h first with fewer positions, C+m second with many more). */
static void test_large_unequal_real_pattern() {
    // Simulate: C+h has 3 positions, C+m has 5 positions
    // ML: 3 values for h (10,20,30), then 5 for m (40,50,60,70,80)
    std::string MM = "C+h,1,2,3;C+m,0,1,1,1,1";
    std::string ML = ",10,20,30,40,50,60,70,80";

    {
        std::string MM_h = MM, ML_h = ML;
        bool dot = false;
        subset_to_Cm_section(MM_h, ML_h, dot, "h");
        check_eq("test_large_unequal / ML for h", ML_h, ",10,20,30");
    }
    {
        std::string MM_m = MM, ML_m = ML;
        bool dot = false;
        subset_to_Cm_section(MM_m, ML_m, dot, "m");
        check_eq("test_large_unequal / ML for m", ML_m, ",40,50,60,70,80");
    }
}

// ---- main ------------------------------------------------------------------

int main() {
    std::cout << "=== test_ont_ml: ML slicing in subset_to_Cm_section ===\n";
    test_unequal_h_before_m();
    test_unequal_m_before_h();
    test_equal_counts();
    test_only_m();
    test_large_unequal_real_pattern();
    std::cout << "=== All tests passed ===\n";
    return 0;
}
