#include "CKKSController.h"

void CKKSController::generate_context_for_bootstrapping(int ring, int levels) {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::SPARSE_ENCAPSULATED);

    int dcrtBits = 36;
    int firstMod = 37;

    depth = levels + FHECKKSRNS::GetBootstrapDepth({4, 3}, lbcrypto::SPARSE_ENCAPSULATED);

    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(ring);
    parameters.SetNumLargeDigits(4);

    parameters.SetBatchSize(ring / 2);

    this->slots = ring / 2;

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    parameters.SetMultiplicativeDepth(depth);

    context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);
    context->Enable(FHE);

    key_pair = context->KeyGen();

    context->EvalMultKeyGen(key_pair.secretKey);
    print_moduli_chain(key_pair.publicKey->GetPublicElements()[0]);
    generate_bootstrapping();
}

void CKKSController::generate_bootstrapping() {
    int slots_bootstrapping = slots;

    context->EvalBootstrapSetup({4, 3}, {0, 0}, slots_bootstrapping, 0, true, true);
    context->EvalBootstrapKeyGen(key_pair.secretKey, slots_bootstrapping);
}

void CKKSController::generate_rotation_key(int index) {
    vector<int> rotations;

    rotations.push_back(index);

    context->EvalRotateKeyGen(key_pair.secretKey, rotations);
}

void CKKSController::generate_rotation_keys(vector<int> indexes) {
    context->EvalRotateKeyGen(key_pair.secretKey, indexes);
}

void CKKSController::generate_rotation_keys_inverse(vector<int> indexes) {
    vector<int> inverted;
    for (size_t i = 0; i < indexes.size(); i++) {
        inverted.push_back(-indexes[i]);
    }
    context->EvalRotateKeyGen(key_pair.secretKey, inverted);
}

void CKKSController::generate_rotations_for_additions(int bits) {
    for (int i = 1; i < bits; i *= 2)
        generate_rotation_key(-i);
}

void CKKSController::generate_rotations_for_bit_length(int bits) {
    for (int i = 1; i < bits; i *= 2)
        generate_rotation_key(i);

    for (int i = 1; i < 256; i *= 2)
        generate_rotation_key(i);

    generate_rotation_key(7);

    generate_rotation_key(1);
    generate_rotation_key(2);
    generate_rotation_key(3);
    generate_rotation_key(4);
    generate_rotation_key(5);
    generate_rotation_key(6);
    generate_rotation_key(bits - 1 - 8);
    generate_rotation_key(bits - 7);
    generate_rotation_key(-bits);
    generate_rotation_key(-bits-1);
}

void CKKSController::generate_rotations_for_multiplications(int bits) {
    if (bits == 8) {
        generate_rotation_keys_inverse(rot_8_bits);
    } else if (bits == 16) {
        generate_rotation_keys_inverse(rot_16_bits);
    } else if (bits == 32) {
        generate_rotation_keys_inverse(rot_32_bits);
    } else if (bits == 64) {
        generate_rotation_keys_inverse(rot_64_bits);
    } else if (bits == 128) {
        generate_rotation_keys_inverse(rot_128_bits);
    } else if (bits == 256) {
        generate_rotation_keys_inverse(rot_256_bits);
    } else {
        cerr << "Unsupported number of bits (" << bits << "), use 16, 32 or 64" << endl;
        return;
    }
}


Ptxt CKKSController::encode(const vector<double> &vec, int lvl) {
    Ptxt p = context->MakeCKKSPackedPlaintext(vec, 1, lvl, nullptr, vec.size());
    p->SetLength(vec.size());

    return p;
}

Ptxt CKKSController::encode(const vector<int> &vec, int lvl) {
    std::vector<std::complex<double>> complex_vec;
    complex_vec.reserve(vec.size());

    std::transform(vec.begin(), vec.end(),
                   std::back_inserter(complex_vec),
                   [](int x) { return std::complex<double>(x, 0.0); });


    Ptxt p = context->MakeCKKSPackedPlaintext(complex_vec, 1, lvl, nullptr, vec.size());
    p->SetLength(vec.size());

    return p;
}

Ptxt CKKSController::encode(double value, int lvl) {
    vector<double> repeated_value;
    for (uint32_t i = 0; i < slots; i++) repeated_value.push_back(value);

    return encode(repeated_value, lvl);
}

Ctxt CKKSController::encrypt(const vector<double> &vec) {
    return encrypt(encode(vec));
}

Ctxt CKKSController::encrypt(const vector<int> &vec) {
    return encrypt(encode(vec));
}

Ctxt CKKSController::encrypt(const vector<double> &vec, int lvl) {
    return encrypt(encode(vec, lvl));
}

Ctxt CKKSController::encrypt(const vector<int> &vec, int lvl) {
    return encrypt(encode(vec, lvl));
}

Ctxt CKKSController::encrypt(const Ptxt &p) {
    return context->Encrypt(p, key_pair.publicKey);
}

Ctxt CKKSController::encrypt_single_int(uint128_t val, int bits, int lvl) {
    vector<int> vec = intToBitsLSB(val, bits);
    append_zeros(vec, context->GetRingDimension() / 2 - bits);
    return encrypt(vec, lvl);
}

Ctxt CKKSController::encrypt_multi_int(vector<uint128_t> val, int bits, int lvl) {
    vector<int> toBeEncoded;

    for (size_t i = 0; i < val.size(); i++) {
        vector<int> vec = intToBitsLSB(val[i], bits);
        append_zeros(vec, bits * bits / 2 - bits);

        toBeEncoded.insert(toBeEncoded.end(), vec.begin(), vec.end());
    }

    return encrypt(toBeEncoded, lvl);
}

vector<double> CKKSController::decode(const Ptxt& p) {
    return p->GetRealPackedValue();
}

Ptxt CKKSController::decrypt(const Ctxt &c) {
    Ptxt p;
    context->Decrypt(key_pair.secretKey, c, &p);

    return p;
}

Ctxt CKKSController::add(const Ctxt &a, const Ctxt &b) {
    return context->EvalAdd(a, b);
}

Ctxt CKKSController::add(const Ctxt &a, const Ptxt &b) {
    Ptxt temp(b);
    return context->EvalAdd(a, temp);
}

Ctxt CKKSController::add(const Ctxt &a, double d) {
    Ptxt temp = encode(d);
    return context->EvalAdd(a, temp);
}

Ctxt CKKSController::add_tree(vector<Ctxt> v) {
    return context->EvalAddMany(v);
}

Ctxt CKKSController::sub(const Ctxt &a, const Ctxt &b) {
    return context->EvalSub(a, b);
}

Ctxt CKKSController::sub(const Ctxt &c, const Ptxt &p) {
    Ptxt temp(p);
    return context->EvalSub(c, temp);
}

Ctxt CKKSController::sub(const Ptxt &p, const Ctxt &c) {
    Ptxt temp(p);
    return context->EvalSub(temp, c);
}

Ctxt CKKSController::sub(double d, const Ctxt &c) {
    //Ptxt temp(p);
    return context->EvalSub(d, c);
}

Ctxt CKKSController::mult(const Ctxt &c, const Ptxt& p) {
    return context->EvalMult(c, p);
}

Ctxt CKKSController::mult(const Ctxt &c, vector<double> p) {
    Ptxt ptxt = encode(p, c->GetLevel());
    return context->EvalMult(c, ptxt);
}

Ctxt CKKSController::mult(const Ctxt &c, vector<int> p) {
    Ptxt ptxt = encode(p, c->GetLevel());
    return context->EvalMult(c, ptxt);
}

Ctxt CKKSController::mult(const Ctxt &c1, const Ctxt &c2) {
    return context->EvalMult(c1, c2);
}

Ctxt CKKSController::mult(const Ctxt &c, double v) {
    return context->EvalMult(c, encode(v));
}

Ctxt CKKSController::square(const Ctxt &c) {
    return context->EvalSquare(c);
}

Ctxt CKKSController::rot(const Ctxt &c, int rotIndex) {
    if (rotIndex == 0) return c;
    return context->EvalRotate(c, rotIndex);
}
Ctxt CKKSController::rot_fast(const Ctxt &c, int rotIndex, shared_ptr<vector<lbcrypto::DCRTPoly>> precomputations) {
    if (rotIndex == 0) return c;
    return context->EvalFastRotation(c, rotIndex, context->GetCyclotomicOrder(), precomputations);
}

vector<int> CKKSController::rot(const vector<int>& vec, int rotIndex) {
    if (vec.empty()) return {};

    int n = vec.size();
    rotIndex = rotIndex % n;
    if (rotIndex == 0) return vec;

    std::vector<int> result = vec;

    if (rotIndex > 0) {
        // Positive shift → left rotation
        std::rotate(result.begin(), result.begin() + rotIndex, result.end());
    } else {
        // Negative shift → right rotation
        rotIndex = -rotIndex;
        std::rotate(result.rbegin(), result.rbegin() + rotIndex, result.rend());
    }

    return result;
}

Ctxt CKKSController::clean_and_reduce(const Ctxt &c) {
    return context->EvalMult(context->EvalSquare(c), context->EvalSquare(context->EvalSub(c, 2)));
}

Ctxt CKKSController::clean(const Ctxt &c) {
    Ctxt sq = context->EvalSquare(c);
    Ctxt t1 = context->EvalMult(c, -2);

    return context->EvalAdd(context->EvalMult(sq, t1), context->EvalMult(sq, 3));
}

Ctxt CKKSController::mod2shallow(const Ctxt &c) {
    return context->EvalSub(context->EvalMult(c, 2), context->EvalSquare(c));
}

Ctxt CKKSController::bintodec(const Ctxt &c, int repetitions) {
    vector<double> mask;

    for (int i = 0; i < repetitions; i++) {
        mask.insert(mask.end(), {1, 2, 4, 8, 0, 0, 0, 0});
    }

    Ctxt res = mult(c, mask);
    res = add(res, rot(res, 1));
    res = add(res, rot(res, 2));

    vector<double> mask2;

    for (int i = 0; i < repetitions; i++) {
        mask2.insert(mask2.end(), {sqrt(1.0 / (225.0 / 2.0)), 0, 0, 0, 0, 0, 0, 0});
    }

    res = mult(res, mask2);
    res = add(res, rot(res, -1));
    res = add(res, rot(res, -2));
    res = add(res, rot(res, -4));

    return res;
}

std::pair<Ctxt, Ctxt> CKKSController::csa3(const Ctxt &a, const Ctxt &b, const Ctxt &c, bool clean_vals) {
    Ctxt S;

    if (clean_vals) {
        S = clean_and_reduce(add(a, b));
        S = mod2shallow(add(S, clean(c)));
    } else {
        S = mod2shallow(add(a, b));
        S = mod2shallow(add(S, c));
    }

    Ctxt C = majoritybit(a, b, c);

    return {S, C};
}

Ctxt CKKSController::multiplier4bits(const Ctxt &a, const Ctxt &b, int repetitions) {
    Ctxt result = mult(a, b);
    result = add(result, -1);

    vector<vector<double>> coeffs;

    coeffs.push_back(read_vector_file("../coeffs/p1-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p2-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p3-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p4-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p5-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p6-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p7-norm-369.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p8-norm-369.txt"));

    Ctxt resultpoly = context->EvalChebyshevSeriesPSBatchRepeated(result, coeffs, -1, 1, repetitions);

    resultpoly = binboot(resultpoly);

    return resultpoly;
}

Ctxt CKKSController::majoritybit(const Ctxt &a, const Ctxt &b, const Ctxt &c) {
    Ctxt total = add(add(a, b), c);
    Ctxt sq = context->EvalSquare(total);
    Ctxt t1 = mult(total, -1.0/3.0);

    return add(add(mult(t1, sq), mult(sq, 3.0/2.0)), mult(total, -7.0/6.0));

}

Ctxt CKKSController::csa4(const Ctxt &a, const Ctxt &b, const Ctxt &c, const Ctxt &d, int bits) {
    Ctxt s1, c1;
    std::tie(s1, c1) = csa3(a, b, c, false);
    c1 = rot(c1, -1);

    Ctxt s2, c2;
    std::tie(s2, c2) = csa3(s1, c1, d, false);
    c2 = rot(c2, -1);

    Ctxt result = add_integer(s2, c2, bits, false);

    return result;
}

Ctxt CKKSController::binary_or(const Ctxt &a, const Ctxt &b) {
    // a+b - a*b
    return sub(add(a, b), mult(a, b));
}


Ctxt CKKSController::add_integer(const Ctxt &a, const Ctxt &b, int bits, bool clean_first) {
    Ctxt p;

    if (clean_first) {
        p = clean_and_reduce(add(a,b));
    } else {
        p = square(sub(a, b));
    }

    Ctxt absum = p->Clone();

    Ctxt g = mult(a, b);

    for (int i = 1; i < bits; i *= 2) {
        Ctxt p_shift = rot(p, -i);
        Ctxt g_shift = rot(g, -i);

        Ctxt pg = mult(p, g_shift);
        g = sub(add(g, pg), mult(p, g));

        if (i < bits - 1) p = mult(p, p_shift);
    }

    g = rot(g, -1);
    Ctxt s = square(sub(absum, g));

    return s;
}

Ctxt CKKSController::sub_integer(const Ctxt &a, const Ctxt &b, int bits, bool clean_first) {
    vector<double> ones;

    int s = context->GetRingDimension() / (bits * bits);

    for (int i = 0; i < s; i++) {
        for (int j = 0; j < bits + 1; j++) {
            ones.push_back(1);
        }
        for (int j = 0; j < bits * bits / 2 - bits - 1; j++) {
            ones.push_back(0);
        }
    }

    Ctxt inverted = context->EvalSub(encrypt(encode(ones, b->GetLevel())), b);
    return add_integer(a, inverted, bits, clean_first);
}

Ctxt CKKSController::binary_mult(const Ctxt &a, const Ctxt &b, int bits, int repetitions) {
    Ctxt result;

    int rep_size = bits * bits / 2;
    int dunn = bits * bits / 8;

    Ctxt a_processed, b_processed;

    if (bits > 8) {
        vector<double> mask_low(a->GetSlots(), 0);
        vector<double> mask_high(a->GetSlots(), 0);

        for (int j = 0; j < repetitions; ++j) {
            for (int i = 0; i < bits / 2; ++i) {
                mask_low[(rep_size * j) + i] = 1;
            }
        }

        for (int j = 0; j < repetitions; ++j) {
            for (int i = bits / 2; i < bits; ++i) {
                mask_high[(rep_size * j) + i] = 1;
            }
        }

        a_processed = mult(a, mask_low);
        Ctxt a_processed_high = mult(a, mask_high);

        a_processed = add(a_processed, rot(a_processed, -(dunn)));
        a_processed = add(a_processed, rot(a_processed_high, -(dunn * 2 - bits / 2)));
        a_processed = add(a_processed, rot(a_processed_high, -(dunn * 3 - bits / 2)));

        b_processed = mult(b, mask_low);
        Ctxt b_processed_high = mult(b, mask_high);
        b_processed = add(b_processed, rot(b_processed, -(dunn * 2)));
        b_processed = add(b_processed, rot(b_processed_high, -(dunn - bits / 2)));
        b_processed = add(b_processed, rot(b_processed_high, -(dunn * 3 - bits / 2)));

    } else {
        a_processed = a->Clone();
        a_processed = add(a_processed, rot(a, -8));
        a_processed = add(a_processed, rot(a, -(8 + 4)));
        a_processed = add(a_processed, rot(a, -(16 + 4)));

        b_processed = b->Clone();
        b_processed = add(b_processed, rot(b, -16));
        b_processed = add(b_processed, rot(b, -4));
        b_processed = add(b_processed, rot(b, -20));
    }

    if (bits == 8) {
        result = multiplier4bits(bintodec(a_processed, repetitions * 4),
                                    bintodec(b_processed, repetitions * 4),
                                    repetitions * 4);
    } else {
        result = binary_mult(a_processed, b_processed, bits / 2, 4 * repetitions);
    }

    int dunn2 = dunn * 2;

    vector<double> mask1(a->GetSlots(), 0.0);

    for (int j = 0; j < repetitions; ++j) {
        for (int i = 0; i < bits; ++i) {
            mask1[(j * rep_size) + i] = 1.0;
            mask1[(j * rep_size) + i + dunn2] = 1.0;
        }
    }

    vector<double> mask2(a->GetSlots(), 0.0);

    for (int j = 0; j < repetitions; ++j) {
        for (int i = 0; i < bits; ++i) {
            mask2[(j * rep_size) + rep_size / 4 + i] = 1.0;
            mask2[(j * rep_size) + rep_size / 4 + i + dunn2] = 1.0;
        }
    }

    Ctxt p1 = mult(result, mask1);
    Ctxt p2 = rot(p1, -(-rep_size / 2 + bits / 2));

    Ctxt p3, p4;

    if (bits == 8) {
        p3 = rot(mult(result, mask2), 16);
    } else {
        p3 = rot(mult(result, mask2), -(-rep_size / 4 + bits / 2));
    }

    if (bits == 8) {
        p4 = rot(p3, -12);
    } else if (bits == 16) {
        p4 = rot(mult(result, mask2), 80);
    } else if (bits == 32) {
        p4 = rot(mult(result, mask2), 352);
    } else if (bits == 64) {
        p4 = rot(mult(result, mask2), 1472);
    }

    result = csa4(p1, p2, p3, p4, bits);

    result = binboot(result);


    return result;
}

Ctxt CKKSController::process_array(const Ctxt& c, const Ctxt& c_processed, const std::vector<std::pair<int,int>>& mask_roll_pairs, int mask_size, int rep, shared_ptr<vector<DCRTPoly>> rot_precomputations) {
    Ctxt c_processed_clone = c_processed->Clone();

    for (auto [start, roll_base] : mask_roll_pairs) {
        int total_size = mask_size * rep;

        vector<int> mask(total_size, 0);

        for (int i = 0; i < rep; ++i) {
            for (int j = 0; j < 4; ++j) {
                mask[i * mask_size + start + j] = 1;
            }
        }

        int shift = roll_base - start;
        Ctxt rolled_ctxt = rot_fast(c, -shift, rot_precomputations);
        vector<int> rolled_mask = rot(mask, -shift);

        c_processed_clone = add(c_processed_clone, mult(rolled_ctxt, rolled_mask));
    }

    return c_processed_clone;
}

Ctxt CKKSController::mul_integer(const Ctxt &a, const Ctxt &b, int bits, int bits_original, int repetitions, int repetitions_original, bool overflow) {
    Ctxt result;

    int rep_size = bits * bits / 2;

    // The size of the basic multiplicator (8 bits)
    int base_mult = 8;

    Ctxt a_processed, b_processed;

    shared_ptr<vector<DCRTPoly>> a_precomputations = context->EvalFastRotationPrecompute(a);
    shared_ptr<vector<DCRTPoly>> b_precomputations = context->EvalFastRotationPrecompute(b);

    if (bits == 8) {
        int mask_size = bits_original * (bits_original / 2);

        vector<int> masklow(a->GetSlots(), 0);
        for (int j = 0; j < repetitions_original; j++) {
            masklow[0 + j * mask_size] = 1;
            masklow[1 + j * mask_size] = 1;
            masklow[2 + j * mask_size] = 1;
            masklow[3 + j * mask_size] = 1;
        }

        a_processed = mult(a, masklow);

        vector<int> maskhigh(a->GetSlots(), 0);

        for (int j = 0; j < repetitions_original; j++) {
            maskhigh[4 + j * mask_size] = 1;
            maskhigh[5 + j * mask_size] = 1;
            maskhigh[6 + j * mask_size] = 1;
            maskhigh[7 + j * mask_size] = 1;
        }

        a_processed = add(a_processed, mult(rot(a, -(16 - 4)), rot(maskhigh, -(16 - 4))));

        if (bits_original > 8) {
            a_processed = process_array(a, a_processed, {{8, 64}, {12, 80}}, mask_size, repetitions_original, a_precomputations);
        }

        if (bits_original > 16) {
            a_processed = process_array(a, a_processed, {{16, 256}, {20, 272}, {24, 320}, {28, 336}}, mask_size, repetitions_original, a_precomputations);
        }

        if (bits_original > 32) {
            a_processed = process_array(a, a_processed,{{32, 1024}, {36, 1040}, {40, 1088}, {44, 1104}, {48, 1280}, {52, 1296}, {56, 1344}, {60, 1360}}, mask_size, repetitions_original, a_precomputations);
        }

        if (bits_original > 64) {
            a_processed = process_array(a, a_processed, {{64, 4096}, {68, 4112}, {72, 4160}, {76, 4176}, {80, 4352}, {84, 4368}, {88, 4416}, {92, 4432}, {96, 5120}, {100, 5136}, {104, 5184}, {108, 5200}, {112, 5376}, {116, 5392}, {120, 5440}, {124, 5456}},mask_size,repetitions_original,a_precomputations);
        }

        if (bits_original > 128) {
            a_processed = process_array(a,a_processed,{{128, 16384}, {132, 16400}, {136, 16448}, {140, 16464},{144, 16640}, {148, 16656}, {152, 16704}, {156, 16720},{160, 17408}, {164, 17424}, {168, 17472}, {172, 17488},{176, 17664}, {180, 17680}, {184, 17728}, {188, 17744},{192, 20480}, {196, 20496}, {200, 20544}, {204, 20560},{208, 20736}, {212, 20752}, {216, 20800}, {220, 20816},{224, 21504}, {228, 21520}, {232, 21568}, {236, 21584},{240, 21760}, {244, 21776}, {248, 21824}, {252, 21840}},mask_size,repetitions_original,a_precomputations);
        }

        if (bits_original > 4) a_processed = add(a_processed, rot(a_processed, -8));
        if (bits_original > 8) a_processed = add(a_processed, rot(a_processed, -32));
        if (bits_original > 16) a_processed = add(a_processed, rot(a_processed, -128));
        if (bits_original > 32) a_processed = add(a_processed, rot(a_processed, -512));
        if (bits_original > 64) a_processed = add(a_processed, rot(a_processed, -2048));
        if (bits_original > 128) a_processed = add(a_processed, rot(a_processed, -8192));

        // B //

        b_processed = mult(b, masklow);
        b_processed = add(b_processed, mult(rot(b, - 4), rot(maskhigh,  - 4)));

        if (bits_original > 8) {
            b_processed = process_array(b, b_processed, {{8, 32}, {12, 40}}, mask_size, repetitions_original, b_precomputations);
        }

        if (bits_original > 16) {
            b_processed = process_array(b, b_processed, {{16, 128}, {20, 136}, {24, 160}, {28, 168}}, mask_size, repetitions_original, b_precomputations);
        }

        if (bits_original > 32) {
            b_processed = process_array(b, b_processed, {{32, 512}, {36, 520}, {40, 544}, {44, 552}, {48, 640}, {52, 648}, {56, 672}, {60, 680}}, mask_size, repetitions_original, b_precomputations);
        }

        if (bits_original > 64) {
            b_processed = process_array(b, b_processed, {{64, 2048}, {68, 2056}, {72, 2080}, {76, 2088}, {80, 2176}, {84, 2184}, {88, 2208}, {92, 2216}, {96, 2560}, {100, 2568}, {104, 2592}, {108, 2600}, {112, 2688}, {116, 2696}, {120, 2720}, {124, 2728}}, mask_size, repetitions_original, b_precomputations);
        }

        if (bits_original > 128) {
            b_processed = process_array(b, b_processed, {{128, 8192}, {132, 8200}, {136, 8224}, {140, 8232}, {144, 8320}, {148, 8328}, {152, 8352}, {156, 8360}, {160, 8704}, {164, 8712}, {168, 8736}, {172, 8744}, {176, 8832}, {180, 8840}, {184, 8864}, {188, 8872}, {192, 10240}, {196, 10248}, {200, 10272}, {204, 10280}, {208, 10368}, {212, 10376}, {216, 10400}, {220, 10408}, {224, 10752}, {228, 10760}, {232, 10784}, {236, 10792}, {240, 10880}, {244, 10888}, {248, 10912}, {252, 10920}}, mask_size, repetitions_original, b_precomputations);
        }


        if (bits_original > 4) b_processed = add(b_processed, rot(b_processed, -16));
        if (bits_original > 8) b_processed = add(b_processed, rot(b_processed, -64));
        if (bits_original > 16) b_processed = add(b_processed, rot(b_processed, -256));
        if (bits_original > 32) b_processed = add(b_processed, rot(b_processed, -1024));
        if (bits_original > 64) b_processed = add(b_processed, rot(b_processed, -4096));
        if (bits_original > 128) b_processed = add(b_processed, rot(b_processed, -16384));


        //auto t = steady_clock::now();
        result = multiplier4bits(bintodec(a_processed, repetitions * 4),
                                 bintodec(b_processed, repetitions * 4), repetitions * 4);
        //print_duration(t, "4 bits multiplier took: ");
    } else {
        result = mul_integer(a, b, bits / 2, bits_original, 4 * repetitions, repetitions_original, overflow);
    }

    int dunn = (bits * bits / base_mult) * 2;

    vector<double> mask1(a->GetSlots(), 0.0);

    for (int j = 0; j < repetitions; ++j) {
        for (int i = 0; i < bits; ++i) {
            mask1[(j * rep_size) + i] = 1.0;
            mask1[(j * rep_size) + i + dunn] = 1.0;
        }
    }

    vector<double> mask2(a->GetSlots(), 0.0);

    for (int j = 0; j < repetitions; ++j) {
        for (int i = 0; i < bits; ++i) {
            mask2[(j * rep_size) + rep_size / 4 + i] = 1.0;
            mask2[(j * rep_size) + rep_size / 4 + i + dunn] = 1.0;
        }
    }

    Ctxt p1 = mult(result, mask1);
    Ctxt p2 = rot(p1, -(-rep_size/2 + bits/2));

    Ctxt p3, p4;

    if (bits == 8) {
        p3 = rot(mult(result, mask2), 16);
    } else {
        p3 = rot(mult(result, mask2), -(-rep_size / 4 + bits / 2));
    }


    if (bits == 8) {
        p4 = rot(p3, -12);
    } else {
        p4 = rot(mult(result, mask2), (((bits - 2) * (3 * bits - 2)) / 8));
    }

    if (!overflow && bits == bits_original) {
        pair<Ctxt, Ctxt> out = csa3(p1, p2, p3);
        result = binboot(add_integer(out.first, rot(out.second, -1), bits));
    } else {
        result = binboot(csa4(p1, p2, p3, p4, bits));
    }

    return result;
}

Ctxt CKKSController::shf_integer(const Ctxt &a, int shift, int bits) {
    vector<double> mask;

    int int_slots = context->GetRingDimension() / (bits * bits);

    for (int i = 0; i < int_slots; i++) {
        for (int j = 0; j < bits; j++) {
            mask.push_back(1);
        }
        for (int j = 0; j < (bits / 2) * (bits / 2) / 2 - bits; j++) {
            mask.push_back(0);
        }
    }

    return rot(mult(a, mask), -shift);
}

Ctxt CKKSController::bit_length(const Ctxt &a, int bits) {
    //TODO write it (similar as below)
}

Ctxt CKKSController::inverse_bit_length(const Ctxt &a, int bits) {
    int step = 1;
    Ctxt result = a->Clone();
    while (step < bits) {
        result = binary_or(result, rot(result, step));
        step *= 2;
    }

    // Now we sum all the ones
    for (int i = 0; i < log2(bits); i++) {
        result = add(result, rot(result, -pow(2, i)));
    }

    //The result is in the last (partial) slot
    vector<double> mask(result->GetSlots()); //TODO mettere batched

    mask[bits - 1] = -1/(bits / 2.0);
    result = mult(result, mask);

    mask[bits - 1] = 1.0;
    result = add(result, encode(mask, result->GetLevel()));

    Ctxt resultclone = result->Clone();

    result = add(result, rot(result, 1));
    result = add(result, rot(result, 2));
    result = add(result, rot(result, 4));
    result = sub(result, rot(resultclone, 7));

    result = rot(result, bits - 7);

    return result;
}

/*
 * Takes as input a 'bits' size input a, and a number of at most 7 bits (128 at most), LSB to MSB in binary, to the RIGHT
 */
Ctxt CKKSController::blind_rotation(const Ctxt &a, const Ctxt &index, int bits) {
    Ctxt result = a->Clone();

    for (int i = 0; i < 7; i++) {
        vector<double> mask(a->GetSlots());
        mask[i] = 1;
        Ctxt current_index = mult(index, encode(mask, index->GetLevel()));
        current_index = rot(current_index, i);
        for (int j = 0; j < log2(bits); j++) {
            //Filling
            current_index = add(current_index, rot(current_index, -pow(2, j)));
        }

        //If condition
        result = add(mult(result, sub(1, current_index)), mult(rot(result, -pow(2, i)), current_index));
    }

    return result;
}

Ctxt CKKSController::div_integer(const Ctxt &num, const Ctxt &den, int bits) {
    int LUT_BITS = 8;

    //This is bits - den.bit_length()
    Ctxt b = inverse_bit_length(den, bits);

    vector<vector<double>> coeffs;
    coeffs.push_back(read_vector_file("../coeffs/p1-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p2-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p3-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p4-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p5-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p6-norm-247-LUT-DIVISION.txt"));
    coeffs.push_back(read_vector_file("../coeffs/p7-norm-247-LUT-DIVISION.txt"));

    for (int i = 0; i < 32-7; i++) coeffs.push_back(read_vector_file("../coeffs/p1-norm-247-LUT-DIVISION.txt")); //Garbage slots

    Ctxt s = context->EvalChebyshevSeriesPSBatchRepeated(b, coeffs, -1, 1, (int)(b->GetSlots() / coeffs.size()));

    s = binboot(s);

    Ctxt den_norm = blind_rotation(den, s, bits); //This performs the den << (bits - s.bit_length()) step

    den_norm = binboot(den_norm);

    Ctxt den_norm_rot = rot(den_norm, bits - 1 - LUT_BITS);

    // Now the index must be in decimal to be given as input to the Chebyshev-LUT
    vector<double> mask(den_norm_rot->GetSlots());
    for (int i = 0; i < LUT_BITS; i++) mask[i] = pow(2, i);
    Ctxt idx = mult(den_norm_rot, mask);

    //N.b. idx \in [0, 256]
    for (int i = 0; i < log2(256); i++) {
        idx = add(idx, rot(idx, pow(2, i)));
    }

    fill(mask.begin(), mask.end(), 0.0);
    mask[0] = 1;
    idx = mult(idx, mask);

    Ctxt idx_masked_clone = idx->Clone();

    for (int i = 0; i < log2(bits); i++) {
        idx = add(idx, rot(idx, -pow(2, i)));
    }

    idx = add(idx, rot(idx_masked_clone, -bits));
    idx = add(idx, rot(idx_masked_clone, -bits-1));

    //LUT output occupies bits + 2 slots, so we repeat idx 'bits + 2' times

    coeffs.clear();
    for (int i = 0; i < bits + 2; i++) coeffs.push_back(read_vector_file("../coeffs/LUTs/" + to_string(bits) + " bits/division/LUT-DIVISION-32-bits-" + to_string(i) + ".txt"));
    for (int i = 0; i < (bits * 2) - (bits + 2); i++) coeffs.push_back(read_vector_file("../coeffs/LUTs/" + to_string(bits) + " bits/division/LUT-DIVISION-32-bits-0.txt")); //Garbage

    Ctxt x = get_context()->EvalChebyshevSeriesPSBatchRepeated(idx, coeffs, 0, 256, (int)(idx->GetSlots() / coeffs.size()));
    x = binboot(x);

    //x is the actual hint, let's go

    for (int iter = 0; iter < 3; iter++) {
        Ctxt term = mul_integer(x, den_norm, bits * 2, bits * 2, 1, 1, true);

        fill(mask.begin(), mask.end(), 0.0);
        for (int i = 0; i < bits*2+1; i++) mask[i] = 1;
        term = sub(encode(mask, term->GetLevel()), term);

        // here we should have a +1, but for convergence is not required i guess? it will get canceled after anyways

        term = rot(term, bits); //Removing the first 'bits' least significative

        x = mul_integer(x, term, bits * 2, bits * 2, 1, 1, true);

        x = rot(x, bits);

    }

    Ctxt result = mul_integer(num, x, bits * 2, bits * 2, 1, 1, true);

    cout << "result before blind: " << cc.print_ints(result, 32 * 4, 1) << endl;

    result = rot(result, bits);

    fill(mask.begin(), mask.end(), 0.0);
    for (int i = 0; i < 32 * 4; i++) mask[i] = 1;
    result = mult(result, mask);

    result = blind_rotation(result, s, bits * 2);

    result = rot(result, bits);
    result = mult(result, mask);

    return result;
}


Ctxt CKKSController::reduce(const Ctxt &c) {
    return context->EvalSub(context->EvalMult(c, 2), context->EvalSquare(c));
}

Ctxt CKKSController::bootstrap(const Ctxt &c) {
    //cout << "Lv boot : " << c->GetLevel() << endl;
    Ctxt cboot = context->EvalBootstrap(c);
    //cout << "Lv after: " << cboot->GetLevel() << endl;
    return cboot;
}

Ctxt CKKSController::binboot(const Ctxt &c) {
    //cout << "Input level : " << c->GetLevel() << endl;
    Ctxt cboot = context->EvalBootstrapStCFirstBits(c);
    cout << "X" << endl;
    //cout << "Output level: " << cboot->GetLevel() << endl;
    return cboot;
}

Ctxt CKKSController::chebyshev(const Ctxt &c, vector<double> coeffs, int a, int b) {
    return context->EvalChebyshevSeries(c, coeffs, a, b);
}

Ctxt CKKSController::chebyshev_batch(const Ctxt &c, vector<vector<double>> coeffs, int a, int b) {
    return context->EvalChebyshevSeriesPSBatch(c, coeffs, a, b);
}

Ctxt CKKSController::chebyshev_batch_rep(const Ctxt &c, vector<vector<double>> coeffs, int a, int b, int repetitions) {
    return context->EvalChebyshevSeriesPSBatchRepeated(c, coeffs, a, b, repetitions);
}

void CKKSController::print(const Ctxt &c) {
    int s = c->GetSlots();

    Ptxt result;
    context->Decrypt(key_pair.secretKey, c, &result);
    result->SetSlots(s);
    vector<double> v = result->GetRealPackedValue();

    cout << "(Level: " << c->GetLevel() << ") [ ";

    for (int i = 0; i < s; i += 1) {
        string segno = "";
        if (v[i] > 0) {
            segno = "";
        } else {
            segno = "-";
            v[i] = -v[i];
        }


        if (static_cast<uint32_t>(i) == slots - 1) {
            if (abs(v[i]) <= 0.0001)
                cout << "0 ]";
            else {
                cout << segno << v[i] << " ]";
            }
        } else {
            if (abs(v[i]) <= 0.0001)
                cout << "0" << " ";
            else
                cout << segno << v[i] << " ";
        }
    }

    cout << endl;
}

void CKKSController::print(const Ctxt &c, int slots) {
    int s = slots;

    Ptxt result;
    context->Decrypt(key_pair.secretKey, c, &result);
    result->SetSlots(s);
    vector<double> v = result->GetRealPackedValue();

    cout << "(Level: " << c->GetLevel() << ") [ ";

    for (int i = 0; i < s; i += 1) {
        string segno = "";
        if (v[i] > 0) {
            segno = "";
        } else {
            segno = "-";
            v[i] = -v[i];
        }


        if (static_cast<uint32_t>(i) == static_cast<uint32_t>(slots - 1)) {
            if (abs(v[i]) <= 0.0001)
                cout << "0 ]";
            else {
                cout << segno << v[i] << " ]";
            }
        } else {
            if (abs(v[i]) <= 0.0001)
                cout << "0" << " ";
            else
                cout << segno << v[i] << " ";
        }
    }

    cout << endl;
}

string CKKSController::print_ints(const Ctxt &c, int bits, int s, bool doublebits) {
    return to_string_uint128(ptxt_to_vec(decode(decrypt(c)), s, bits, doublebits));
}

void CKKSController::print_moduli_chain(const lbcrypto::DCRTPoly &poly) {
    cout << "log(N): " << std::round(log(context->GetRingDimension())/log(2)) << ", circuit depth: " << depth << endl;
}

CryptoContext<DCRTPoly> CKKSController::get_context() {
    return context;
}
