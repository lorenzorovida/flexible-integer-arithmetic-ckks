#include <iostream>
#include "openfhe.h"
#include "CKKSController.h"
#include "chrono"
#include <functional>
#include "Utils.h"
#include "Logger.h"

using namespace lbcrypto;
using namespace std;
using namespace chrono;

CKKSController cc;
int ring_size = 12;
int verbose = 2;
int wordsize = 32;

bool test = false;
bool input_mode = false;

void read_arguments(int argc, char* argv[]);
void random_operations(int bits);
void random_operations_batched(int bits);

int main(int argc, char* argv[]) {
    read_arguments(argc, argv);

    // Con 13 levels 256-bits
    cc.generate_context_for_bootstrapping(1 << ring_size, 13);
    cc.generate_rotations_for_additions(wordsize * 2);
    cc.generate_rotations_for_multiplications(wordsize * 2);
    cc.generate_rotations_for_bit_length(wordsize);

    /*
     * Experiments
     */
    Ctxt numerator = cc.encrypt_multi_int({3200232430, 0, 0, 0}, 32, 12);
    Ctxt inpt = cc.encrypt_multi_int({10000003,200,300,400}, 32, 12);
    Ctxt b = cc.inverse_bit_length(inpt, 32);
    cc.print(b, 32);


    vector<vector<double>> coeffs;
    coeffs.push_back(read_vector_file("../coeffs/p1-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p2-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p3-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p4-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p5-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p6-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p7-norm-247-LUT-DIVISION.txt"));

    for (int i = 0; i < 32-7; i++) {
        // Garbage
        coeffs.push_back(read_vector_file("../coeffs/p1-norm-247-LUT-DIVISION.txt"));
    }


    Ctxt s = cc.get_context()->EvalChebyshevSeriesPSBatchRepeated(b, coeffs, -1, 1, (int)(b->GetSlots() / coeffs.size()));

    s = cc.binboot(s);

    cout << "s_bin: ";
    cc.print(s, 32);


    Ctxt b_norm = cc.blind_rotation(inpt, s, 32);
    cout << "b_norm: ";
    cc.print(b_norm, 32*2);

    b_norm = cc.binboot(b_norm);

    int LUT_BITS = 8;
    Ctxt b_norm_rot = cc.rot(b_norm, 32 - 1 - LUT_BITS);

    cout << "b_norm_rot: ";
    cc.print(b_norm_rot, 32*2);




    vector<double> mask_idx(b_norm_rot->GetSlots());
    for (int i = 0; i < LUT_BITS; i++) mask_idx[i] = pow(2, i);
    Ctxt idx = cc.mult(b_norm_rot, mask_idx);

    //Idx is in decimal
    for (int i = 0; i < log2(256); i++) {
        idx = cc.add(idx, cc.rot(idx, pow(2, i)));
    }

    vector<double> mask_first(idx->GetSlots());
    mask_first[0] = 1;
    idx = cc.mult(idx, mask_first);

    Ctxt idx_masked_clone = idx->Clone();

    for (int i = 0; i < log2(32); i++) {
        idx = cc.add(idx, cc.rot(idx, -pow(2, i)));
    }

    idx = cc.add(idx, cc.rot(idx_masked_clone, -32));
    idx = cc.add(idx, cc.rot(idx_masked_clone, -32-1));

    cout << "idx (decimal): ";
    cc.print(idx, 32+2);
    cout << idx->GetLevel() << endl;

    coeffs.clear();
    for (int i = 0; i < 32 + 2; i++) {
        coeffs.push_back(read_vector_file("../coeffs/LUTs/32 bits/division/LUT-DIVISION-32-bits-" + to_string(i) + ".txt"));
    }
    for (int i = 0; i < (32 * 2) - (32 + 2); i++) {
        //Garbage
        coeffs.push_back(read_vector_file("../coeffs/LUTs/32 bits/division/LUT-DIVISION-32-bits-0.txt"));
    }



    Ctxt binx = cc.get_context()->EvalChebyshevSeriesPSBatchRepeated(idx, coeffs, 0, 256, (int)(idx->GetSlots() / coeffs.size()));
    binx = cc.binboot(binx);

    cout << "b_norm: " << cc.print_ints(b_norm, 64, 1) << endl;
    cout << "x: " << cc.print_ints(binx, 64, 1) << endl;

    for (int iter = 0; iter < 3; iter++) {

        Ctxt term = cc.mul_integer(b_norm, binx, 64, 64, 1, 1, true);

        cout << "term1: " << cc.print_ints(term, 128, 1) << endl;


        vector<double> vl(term->GetSlots());
        for (int i = 0; i < 32*2+1; i++) vl[i] = 1;
        term = cc.sub(cc.encode(vl, term->GetLevel()), term);

        // here we should have a +1, but for convergence is not required i guess? it will get canceled after anyways

        cout << "term2: " << cc.print_ints(term, 32 * 4, 1) << endl;

        term = cc.rot(term, 32); //Tolgo i 32 bit meno significativi..

        binx = cc.mul_integer(binx, term, 32 * 2, 32 * 2, 1, 1, true);

        cout << "x1: " << cc.print_ints(binx, 32 * 2, 1) << endl;

        binx = cc.rot(binx, 32);

        cout << "x2: " << cc.print_ints(binx, 32 * 2, 1) << endl;

    }

    Ctxt result = cc.mul_integer(numerator, binx, 64, 64, 1, 1, true);

    cout << "result before blind: " << cc.print_ints(result, 32 * 4, 1) << endl;

    result = cc.rot(result, 32);
    vector<double> firstnmask(result->GetSlots());
    for (int i = 0; i < 32 * 4; i++) firstnmask[i] = 1;
    result = cc.mult(result, firstnmask);

    cout << "result before blind: " << cc.print_ints(result, 32 * 4, 1) << endl;


    cc.print(s, 32 * 2);

    result = cc.blind_rotation(result, s, 32 * 2);

    cout << "result before blind: " << cc.print_ints(result, 32 * 2, 1) << endl;

    result = cc.rot(result, 32);
    result = cc.mult(result, firstnmask);


    cout << "result before blind: ";
    cout << cc.print_ints(result, 9, 1);


    exit(0);

    /*
     * End experiments
     */

    if (test) {
        cout << "Keygen works, you are good to go to use the program :)" << endl;
        return 0;
    }

    if (input_mode) {
        random_operations(0);
    } else {
        random_operations_batched(wordsize);
    }
}

void random_operations_batched(int bits) {
    int slots = cc.get_context()->GetRingDimension() / (bits * bits);

    vector<uint128_t> a;
    vector<uint128_t> b;

    Logger log(verbose);

    log.yellow_bold(1) << endl << "Running batched " << bits << "-bits operations experiment!" << endl;

    for (int i = 0; i < slots; i++) {
        a.push_back(random_number(bits));
        b.push_back(random_number(bits));
    }

    log(1) << "a: " << to_string_uint128(a) << endl << "b: " << to_string_uint128(b) << endl << endl;

    Ctxt c1 = cc.encrypt_multi_int(a, bits, 12);
    Ctxt c2 = cc.encrypt_multi_int(b, bits, 12);

    auto time = steady_clock::now();

    Ctxt csum = cc.binboot(cc.add_integer(c1, c2, bits));
    log.info(1) << "Addition (a + b)" << endl;
    log(2) << "Expected: " << to_string_uint128(add_simd(a, b)) << endl;
    log(2) << "Obtained: " << cc.print_ints(csum, bits + 1, slots) << endl;
    if (verbose >= 1) print_duration(time, "Addition took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt csub = cc.binboot(cc.sub_integer(c1, c2, bits));
    log.info(1) << "Comparison (a ≤ b)" << endl;
    log(2) << "Expected: " << comp_simd(a, b) << endl;
    log(2) << "Obtained: " << last_bits(cc.decode(cc.decrypt(csub)), slots, bits) << endl;
    if (verbose >= 1) print_duration(time, "Comparison took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt cmultmod = cc.mul_integer(c1, c2, bits, bits, slots, slots, false);

    log.info(1) << "Multiplication (a * b) % 2^n" << endl;
    log(2) << "Expected: " << to_string_uint128(mul_simd(a, b, bits)) << endl;
    log(2) << "Obtained: " << cc.print_ints(cmultmod, bits, slots) << endl;
    if (verbose >= 1) print_duration(time, "Multiplication took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt cmult = cc.mul_integer(c1, c2, bits, bits, slots, slots, true);

    if (bits > 64) {
        log.info(1) << "Multiplication with overflow (a * b), Warning: results will be inaccurate as we only have access to uint128_t :-(" << endl;
    } else {
        log.info(1) << "Multiplication with overflow (a * b)" << endl;
    }
    log(2) << "Expected: " << to_string_uint128(mul_simd(a, b, 2 * bits)) << endl;
    log(2) << "Obtained: " << cc.print_ints(cmult, bits, slots, true) << endl;

    if (verbose >= 1) print_duration(time, "Multiplication took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt cshift = cc.rot(cmult, -2);
    log.info(1) << "Logical shift (a * b << 2)" << endl;
    log(2) << "Expected: " << to_string_uint128(shift_simd(mul_simd(a, b, bits), 2)) << endl;
    log(2) << "Obtained: " << cc.print_ints(cshift, bits + 2, slots) << endl;
    if (verbose >= 1) print_duration(time, "Logical shift took: ");
    log(1) << endl;
}

void random_operations(int bits) {
    //uint128_t a = random_number(bits);
    //uint128_t b = random_number(bits);

    cout << "Insert the desired number of bits (8, 16, 32, 64, 128, 256): " << endl;
    cin >> bits;

    if (bits != 8 && bits != 16 && bits != 32 && bits != 64 && bits != 128 && bits != 256) {
        cerr << "The amount of bits (" << wordsize << ") is not supported. Pick one out of (8, 16, 32, 64, 128, 256)" << endl;
        return;
    }

    string a_str, b_str;
    uint128_t a = 0, b = 0;

    cout << "Insert A: ";
    cin >> a_str;
    for (char c : a_str) {
        if (c >= '0' && c <= '9') {
            a = a * 10 + (c - '0');
        }
    }

    cout << "Insert B: ";
    cin >> b_str;
    for (char c : b_str) {
        if (c >= '0' && c <= '9') {
            b = b * 10 + (c - '0');
        }
    }

    Logger log(verbose);

    log(1) << endl
           << "Running single " << bits << "-bits operations experiment!"
           << endl;

    log(1) << "a: " << to_string_uint128(a) << ", b: " << to_string_uint128(b) << endl << endl;

    Ctxt c1 = cc.encrypt_single_int(a, bits, 12);
    Ctxt c2 = cc.encrypt_single_int(b, bits, 12);

    auto time = steady_clock::now();

    Ctxt csum = cc.binboot(cc.add_integer(c1, c2, bits));
    log(1) << "Addition (a + b)" << endl;
    log(2) << "Expected: " << to_string_uint128(a + b) << endl;
    log(2) << "Obtained: " << to_string_uint128(bits_to_int128(cc.decode(cc.decrypt(csum)), bits + 1)) << endl;
    if (verbose >= 1) print_duration(time, "Addition took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt csub = cc.binboot(cc.sub_integer(c1, c2, bits));
    log(1) << "Comparison (a ≤ b)" << endl;
    log(2) << "Expected: " << (a <= b) << endl;
    log(2) << "Obtained: " << cc.decode(cc.decrypt(csub))[bits] << endl;
    if (verbose >= 1) print_duration(time, "Comparison took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt cmultmod = cc.mul_integer(c1, c2, bits, bits, 1, 1, false);

    log(1) << "Multiplication (a * b) % 2^n" << endl;
    log(2) << "Expected (" << bits << " bits): " << to_string_uint128((a * b) & ((uint128_t(1) << bits) - 1)) << endl;
    log(2) << "Obtained (" << bits << " bits): " << to_string_uint128(bits_to_int128(cc.decode(cc.decrypt(cmultmod)), bits)) << endl;
    if (verbose >= 1) print_duration(time, "Multiplication took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt cmult = cc.mul_integer(c1, c2, bits, bits, 1, 1, true);

    if (bits > 64) {
        log(1) << "Multiplication with overflow (a * b), Warning: results will be inaccurate as we only have access to uint128_t :-(" << endl;
    } else {
        log(1) << "Multiplication with overflow (a * b)" << endl;
    }
    log(2) << "Expected (" << bits * 2 << " bits): " << to_string_uint128((uint128_t)(a * b)) << endl;
    log(2) << "Obtained (" << bits * 2 << " bits): " << to_string_uint128(bits_to_int128(cc.decode(cc.decrypt(cmult)), bits * 2)) << endl;
    if (verbose >= 1) print_duration(time, "Multiplication took: ");
    log(1) << "-----" << endl;

    time = steady_clock::now();

    Ctxt cshift = cc.rot(cmult, -2);

    log(1) << "Logical shift (a * b << 2)" << endl;
    log(2) << "Expected: " << to_string_uint128(((uint128_t)(a * b) << 2) & ((uint128_t(1) << bits) - 1)) << endl;
    log(2) << "Obtained: " << to_string_uint128(bits_to_int128(cc.decode(cc.decrypt(cshift)), bits + 2)) << endl;
    if (verbose >= 1) print_duration(time, "Logical shift took: ");

    log(1) << endl << endl;
}

void read_arguments(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--ring" && i + 1 < argc) {
            ring_size = stoi(argv[i + 1]);
            ++i;
        }
        if (arg == "--verbose" && i + 1 < argc) {
            verbose = stoi(argv[i + 1]);
            ++i;
        }
        if (arg == "--bits" && i + 1 < argc) {
            wordsize = stoi(argv[i + 1]);

            if (wordsize != 8 && wordsize != 16 && wordsize != 32 && wordsize != 64 && wordsize != 128 &&
                wordsize != 256) {
                cerr << "The amount of bits (" << wordsize
                     << ") is not supported. Pick one out of (8, 16, 32, 64, 128, 256)" << endl;
            }

            ++i;
        }
        if (arg == "--test") {
            cout << "The program has been compiled and linked successfully, now checking if keygen works..." << endl;
            test = true;
        }
        if (arg == "--input") {
            input_mode = true;
        }
    }
}