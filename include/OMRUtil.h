#pragma once

#include "PVWToBFVSeal.h"
#include "SealUtils.h"
#include "retrieval.h"
#include "client.h"
#include "LoadAndSaveUtils.h"
#include "MathUtil.h"
#include "global.h"
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <thread>

using namespace seal;


void choosePertinentMsg(int numOfTransactions, int pertinentMsgNum, vector<int>& pertinentMsgIndices, prng_seed_type& seed) {
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    RandomToStandardAdapter engine(rng->create());
    uniform_int_distribution<uint64_t> dist(0, numOfTransactions - 1);
    for (int i = 0; i < pertinentMsgNum; i++) {
        auto temp = dist(engine);
        while(find(pertinentMsgIndices.begin(), pertinentMsgIndices.end(), temp) != pertinentMsgIndices.end()){
            temp = dist(engine);
        }
        pertinentMsgIndices.push_back(temp);
    }
    sort(pertinentMsgIndices.begin(), pertinentMsgIndices.end());
    // pertinentMsgIndices.push_back(0);

    cout << "Expected Message Indices: " << pertinentMsgIndices << endl;
}


vector<vector<uint64_t>> preparingTransactionsFormal(vector<int>& pertinentMsgIndices, PVWpk& pk, int numOfTransactions, int pertinentMsgNum,
                                                      const PVWParam& params, int partySize = 1) {


    vector<vector<uint64_t>> ret;
    vector<int> zeros(params.ell, 0);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    choosePertinentMsg(numOfTransactions, pertinentMsgNum, pertinentMsgIndices, seed);

    for(int i = 0; i < numOfTransactions; i++){
        PVWCiphertext tempclue;

        // create clues with new SK for the rest of messages in the same group
        for (int p = 0; p < partySize - 1; p++) {
            PVWCiphertext tempclue;
            auto sk2 = PVWGenerateSecretKey(params);
            PVWEncSK(tempclue, zeros, sk2, params);
            saveClues(tempclue, i*partySize + p);
        }

        // w.l.o.g assume the index of recipient within party is |partySize - 1|, i.e., the last in the group
        if(find(pertinentMsgIndices.begin(), pertinentMsgIndices.end(), i) != pertinentMsgIndices.end()) {
            PVWEncPK(tempclue, zeros, pk, params);
            ret.push_back(loadDataSingle(i));
            expectedIndices.push_back(uint64_t(i));
        }
        else
        {
            auto sk2 = PVWGenerateSecretKey(params);
            PVWEncSK(tempclue, zeros, sk2, params);
        }
        saveClues(tempclue, i*partySize + partySize - 1);
    }
    return ret;
}

// Phase 1, obtaining PV's
Ciphertext (vector<PVWCiphertext>& SICPVW, vector<Ciphertext> switchingKey, const RelinKeys& relin_keys,
                            const GaloisKeys& gal_keys, const size_t& degree, const SEALContext& context, const PVWParam& params,
                            const int numOfTransactions, const int partialSize = 0) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    time_start = chrono::high_resolution_clock::now();

    Evaluator evaluator(context);
    
    vector<Ciphertext> packedSIC(params.ell);
    computeBplusASPVWOptimized(packedSIC, SICPVW, switchingKey, gal_keys, context, params, partialSize);

    int rangeToCheck = 850; // range check is from [-rangeToCheck, rangeToCheck-1]
    newRangeCheckPVW(packedSIC, rangeToCheck, relin_keys, degree, context, params);

    time_end = chrono::high_resolution_clock::now();
    cout << "\n Phase1 time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start) << "\n";

    return packedSIC[0];
}


// Phase 1, obtaining PV's based on encrypted targetId
// used in GOMR1/2_ObliviousMultiplexer_BFV
Ciphertext serverOperations1obtainPackedSICWithCluePoly(vector<Ciphertext> switchingKey, const RelinKeys& relin_keys, const GaloisKeys& gal_keys,
                                                        const size_t& degree, const SEALContext& context, const PVWParam& params,
                                                        const int numOfTransactions, uint64_t *total_plain_ntt, uint64_t *total_load) {
    Evaluator evaluator(context);
    
    vector<Ciphertext> packedSIC(params.ell);
    computeBplusASPVWOptimizedWithCluePoly(packedSIC, switchingKey, relin_keys, gal_keys, context, params, total_plain_ntt, total_load);

    int rangeToCheck = 850; // range check is from [-rangeToCheck, rangeToCheck-1]
    newRangeCheckPVW(packedSIC, rangeToCheck, relin_keys, degree, context, params);

    return packedSIC[0];
}

// Phase 1, obtaining PV's based on encrypted secret SK and shared SK
// used in GOMR1/2_FG
Ciphertext serverOperations1obtainPackedSICWithFixedGroupClue(vector<vector<int>>& fgClues, vector<Ciphertext> switchingKey, const RelinKeys& relin_keys,
                            const GaloisKeys& gal_keys, const size_t& degree, const SEALContext& context, const PVWParam& params,
                            const int numOfTransactions, const int partialSize = partial_size_glb) {
    Evaluator evaluator(context);
    
    vector<Ciphertext> packedSIC(params.ell);
    computeBplusASPVWOptimizedWithFixedGroupClue(packedSIC, fgClues, switchingKey, gal_keys, context, params, partialSize);

    int rangeToCheck = 850; // range check is from [-rangeToCheck, rangeToCheck-1]
    newRangeCheckPVW(packedSIC, rangeToCheck, relin_keys, degree, context, params);

    return packedSIC[0];
}

// Phase 2, retrieving
void serverOperations2therest(Ciphertext& lhs, vector<vector<int>>& bipartite_map, Ciphertext& rhs,
                        Ciphertext& packedSIC, const vector<vector<uint64_t>>& payload, const RelinKeys& relin_keys, const GaloisKeys& gal_keys,
                        const size_t& degree, const SEALContext& context, const SEALContext& context2, const PVWParam& params, const int numOfTransactions, 
                        int& counter, int partySize = 1, const int payloadSize = 306){

    Evaluator evaluator(context);
    int step = 32; // simply to save memory so process 32 msgs at a time
    
    bool expandAlter = true;
    
    for(int i = counter; i < counter+numOfTransactions; i += step){
        vector<Ciphertext> expandedSIC;
        // step 1. expand PV
        if(expandAlter)
            expandSIC_Alt(expandedSIC, packedSIC, gal_keys, gal_keys_last, int(degree), context, context2, step, i-counter);
        else
            expandSIC(expandedSIC, packedSIC, gal_keys, gal_keys_last, int(degree), context, context2, step, i-counter);

        // transform to ntt form for better efficiency especially for the last two steps
        for(size_t j = 0; j < expandedSIC.size(); j++)
            if(!expandedSIC[j].is_ntt_form())
                evaluator.transform_to_ntt_inplace(expandedSIC[j]);

        // step 2. deterministic retrieval
        deterministicIndexRetrieval(lhs, expandedSIC, context, degree, i, partySize);

        // step 3-4. multiply weights and pack them
        // The following two steps are for streaming updates
        vector<vector<Ciphertext>> payloadUnpacked;
        payloadRetrievalOptimizedwithWeights(payloadUnpacked, payload, bipartite_map_glb, weights_glb, expandedSIC, context, degree, i, i - counter);
        // Note that if number of repeatitions is already set, this is the only step needed for streaming updates
        payloadPackingOptimized(rhs, payloadUnpacked, bipartite_map_glb, degree, context, gal_keys, i);   
    }
    if(lhs.is_ntt_form())
        evaluator.transform_from_ntt_inplace(lhs);
    if(rhs.is_ntt_form())
        evaluator.transform_from_ntt_inplace(rhs);

    counter += numOfTransactions;
}

// Phase 2, retrieving for OMR3
void serverOperations3therest(vector<Ciphertext>& lhsCounter, vector<vector<int>>& bipartite_map, Ciphertext& rhs,
                        Ciphertext& packedSIC, const vector<vector<uint64_t>>& payload, const RelinKeys& relin_keys, const GaloisKeys& gal_keys, const PublicKey& public_key,
                        const size_t& degree, const SEALContext& context, const SEALContext& context2, const PVWParam& params, const int numOfTransactions, 
                        int& counter, int numberOfCt = 1, int partySize = 1, int slotPerBucket = 3, const int payloadSize = 306){

    Evaluator evaluator(context);

    int step = 32;
    for(int i = counter; i < counter+numOfTransactions; i += step){
        // step 1. expand PV
        vector<Ciphertext> expandedSIC;
        expandSIC_Alt(expandedSIC, packedSIC, gal_keys, gal_keys_last, int(degree), context, context2, step, i-counter);
        // transform to ntt form for better efficiency for all of the following steps
        for(size_t j = 0; j < expandedSIC.size(); j++)
            if(!expandedSIC[j].is_ntt_form())
                evaluator.transform_to_ntt_inplace(expandedSIC[j]);

        // step 2. randomized retrieval
        randomizedIndexRetrieval_opt(lhsCounter, expandedSIC, context, public_key, i, degree,
                                     repetition_glb, numberOfCt, num_bucket_glb, partySize, slotPerBucket);
        // step 3-4. multiply weights and pack them
        // The following two steps are for streaming updates
        vector<vector<Ciphertext>> payloadUnpacked;
        payloadRetrievalOptimizedwithWeights(payloadUnpacked, payload, bipartite_map_glb, weights_glb, expandedSIC, context, degree, i, i-counter);
        // Note that if number of repeatitions is already set, this is the only step needed for streaming updates
        payloadPackingOptimized(rhs, payloadUnpacked, bipartite_map_glb, degree, context, gal_keys, i);
    }
    for(size_t i = 0; i < lhsCounter.size(); i++){
            evaluator.transform_from_ntt_inplace(lhsCounter[i]);
    }
    if(rhs.is_ntt_form())
        evaluator.transform_from_ntt_inplace(rhs);
    
    counter += numOfTransactions;
}


vector<vector<long>> receiverDecoding(Ciphertext& lhsEnc, vector<vector<int>>& bipartite_map, Ciphertext& rhsEnc,
                        const size_t& degree, const SecretKey& secret_key, const SEALContext& context, const int numOfTransactions, int partySize = 1,
                        int seed = 3, const int payloadUpperBound = 306, const int payloadSize = 306){

    // 1. find pertinent indices
    map<int, pair<int, int>> pertinentIndices;
    decodeIndices(pertinentIndices, lhsEnc, numOfTransactions, degree, secret_key, context, partySize);
    cout << "Pertinent message indices found with its group PV value: " << endl;
    for (map<int, pair<int, int>>::iterator it = pertinentIndices.begin(); it != pertinentIndices.end(); it++)
    {
        cout << it->first << "," << it->second.second << "  ";
    }
    cout << endl;

    // 2. forming rhs
    vector<vector<int>> rhs;
    vector<Ciphertext> rhsEncVec{rhsEnc};
    formRhs(rhs, rhsEncVec, secret_key, degree, context, OMRtwoM);

    // 3. forming lhs
    vector<vector<int>> lhs;
    formLhsWeights(lhs, pertinentIndices, bipartite_map_glb, weights_glb, 0, OMRtwoM);

    // 4. solving equation
    auto newrhs = equationSolving(lhs, rhs, payloadSize);

    return newrhs;
}

vector<vector<long>> receiverDecodingOMR3(vector<Ciphertext>& lhsCounter, vector<vector<int>>& bipartite_map, Ciphertext& rhsEnc,
                        const size_t& degree, const SecretKey& secret_key, const SEALContext& context, const int numOfTransactions,
                        int partySize = 1, int slot_per_bucket = 3, int seed = 3, const int payloadUpperBound = 306, const int payloadSize = 306){
    // 1. find pertinent indices
    map<int, pair<int, int>> pertinentIndices;
    decodeIndicesRandom_opt(pertinentIndices, lhsCounter, secret_key, context, partySize, slot_per_bucket);
    for (map<int, pair<int, int>>::iterator it = pertinentIndices.begin(); it != pertinentIndices.end(); it++)
    {
        cout << it->first << "," << it->second.second << "  ";
    }
    cout << std::endl;

    // 2. forming rhs
    vector<vector<int>> rhs;
    vector<Ciphertext> rhsEncVec{rhsEnc};
    formRhs(rhs, rhsEncVec, secret_key, degree, context, OMRthreeM);

    // 3. forming lhs
    vector<vector<int>> lhs;
    formLhsWeights(lhs, pertinentIndices, bipartite_map_glb, weights_glb, 0, OMRthreeM);

    // 4. solving equation
    auto newrhs = equationSolving(lhs, rhs, payloadSize);

    return newrhs;
}

// to check whether the result is as expected
bool checkRes(vector<vector<uint64_t>> expected, vector<vector<long>> res){
    for(size_t i = 0; i < expected.size(); i++){
        bool flag = false;
        for(size_t j = 0; j < res.size(); j++){
            if(expected[i][0] == uint64_t(res[j][0])){
                if(expected[i].size() != res[j].size())
                {
                    cerr << "expected and res length not the same" << endl;
                    return false;
                }
                for(size_t k = 1; k < res[j].size(); k++){
                    if(expected[i][k] != uint64_t(res[j][k]))
                        break;
                    if(k == res[j].size() - 1){
                        flag = true;
                    }
                }
            }
        }
        if(!flag)
            return false;
    }
    return true;
}

// check OMD detection key size
// We are:
//      1. packing PVW sk into ell ciphertexts
//      2. using seed mode in SEAL
void OMDlevelspecificDetectKeySize(){
    auto params = PVWParam(450, 65537, 1.3, 16000, 4); 
    auto sk = PVWGenerateSecretKey(params);
    cout << "Finishing generating sk for PVW cts\n";
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = poly_modulus_degree_glb;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    auto coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 28, 
                                                                            39, 60, 60, 60, 
                                                                            60, 60, 60, 60, 60, 60,
                                                                            32, 30, 60 });
    parms.set_coeff_modulus(coeff_modulus);
    parms.set_plain_modulus(65537);

	prng_seed_type seed;
    for (auto &i : seed)
    {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    parms.set_random_generator(rng);

    SEALContext context(parms, true, sec_level_type::none);
    print_parameters(context); 
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    BatchEncoder batch_encoder(context);
    GaloisKeys gal_keys;

    seal::Serializable<PublicKey> pk = keygen.create_public_key();
	seal::Serializable<RelinKeys> rlk = keygen.create_relin_keys();
	stringstream streamPK, streamRLK, streamRTK;
    auto reskeysize = pk.save(streamPK);
	reskeysize += rlk.save(streamRLK);
	reskeysize += keygen.create_galois_keys(vector<int>({1})).save(streamRTK);

    public_key.load(context, streamPK);
    relin_keys.load(context, streamRLK);
    gal_keys.load(context, streamRTK); 
	vector<Ciphertext> switchingKeypacked = omr::generateDetectionKey(context, poly_modulus_degree, public_key, secret_key, sk, params);
	stringstream data_stream;
    for(size_t i = 0; i < switchingKeypacked.size(); i++){
        reskeysize += switchingKeypacked[i].save(data_stream);
    }
    cout << "Detection Key Size: " << reskeysize << " bytes" << endl;
}

// check OMR detection key size
// We are:
//      1. packing PVW sk into ell ciphertexts
//      2. use level-specific rot keys
//      3. using seed mode in SEAL
void levelspecificDetectKeySize(){
    auto params = PVWParam(450, 65537, 1.3, 16000, 4); 
    auto sk = PVWGenerateSecretKey(params);
    cout << "Finishing generating sk for PVW cts\n";

    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = poly_modulus_degree_glb;
    auto degree = poly_modulus_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    auto coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 28, 
                                                                            60, 60, 60, 60, 60, 
                                                                            60, 60, 60, 60, 60, 60,
                                                                            60, 30, 60 });
    parms.set_coeff_modulus(coeff_modulus);
    parms.set_plain_modulus(65537);


	prng_seed_type seed;
    for (auto &i : seed)
    {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    parms.set_random_generator(rng);

    SEALContext context(parms, true, sec_level_type::none);
    print_parameters(context); 
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    BatchEncoder batch_encoder(context);
    GaloisKeys gal_keys;

    vector<int> steps = {0};
    for(int i = 1; i < int(poly_modulus_degree/2); i *= 2){
	    steps.push_back(i);
    }

    stringstream lvlRTK, lvlRTK2;
    /////////////////////////////////////// Level specific keys
    vector<Modulus> coeff_modulus_next = coeff_modulus;
    coeff_modulus_next.erase(coeff_modulus_next.begin() + 4, coeff_modulus_next.end()-1);
    EncryptionParameters parms_next = parms;
    parms_next.set_coeff_modulus(coeff_modulus_next);
    parms_next.set_random_generator(rng);
    SEALContext context_next = SEALContext(parms_next, true, sec_level_type::none);

    SecretKey sk_next;
    sk_next.data().resize(coeff_modulus_next.size() * degree);
    sk_next.parms_id() = context_next.key_parms_id();
    util::set_poly(secret_key.data().data(), degree, coeff_modulus_next.size() - 1, sk_next.data().data());
    util::set_poly(
        secret_key.data().data() + degree * (coeff_modulus.size() - 1), degree, 1,
        sk_next.data().data() + degree * (coeff_modulus_next.size() - 1));
    KeyGenerator keygen_next(context_next, sk_next); 
    vector<int> steps_next = {0,32,64,128,256,512,1024,2048,4096,8192};
    auto reskeysize = keygen_next.create_galois_keys(steps_next).save(lvlRTK);
        //////////////////////////////////////
    vector<Modulus> coeff_modulus_last = coeff_modulus;
    coeff_modulus_last.erase(coeff_modulus_last.begin() + 3, coeff_modulus_last.end()-1);
    EncryptionParameters parms_last = parms;
    parms_last.set_coeff_modulus(coeff_modulus_last);
    parms_last.set_random_generator(rng);
    SEALContext context_last = SEALContext(parms_last, true, sec_level_type::none);

    SecretKey sk_last;
    sk_last.data().resize(coeff_modulus_last.size() * degree);
    sk_last.parms_id() = context_last.key_parms_id();
    util::set_poly(secret_key.data().data(), degree, coeff_modulus_last.size() - 1, sk_last.data().data());
    util::set_poly(
        secret_key.data().data() + degree * (coeff_modulus.size() - 1), degree, 1,
        sk_last.data().data() + degree * (coeff_modulus_last.size() - 1));
    KeyGenerator keygen_last(context_last, sk_last); 
    vector<int> steps_last = {1,2,4,8,16};
    reskeysize += keygen_last.create_galois_keys(steps_last).save(lvlRTK2);
    //////////////////////////////////////

    seal::Serializable<PublicKey> pk = keygen.create_public_key();
	seal::Serializable<RelinKeys> rlk = keygen.create_relin_keys();
	stringstream streamPK, streamRLK, streamRTK;
    reskeysize += pk.save(streamPK);
	reskeysize += rlk.save(streamRLK);
	reskeysize += keygen.create_galois_keys(vector<int>({1})).save(streamRTK);

    public_key.load(context, streamPK);
    relin_keys.load(context, streamRLK);
    gal_keys.load(context, streamRTK); 
	vector<Ciphertext>  switchingKeypacked = omr::generateDetectionKey(context, poly_modulus_degree, public_key, secret_key, sk, params);
	stringstream data_stream;
    for(size_t i = 0; i < switchingKeypacked.size(); i++){
        reskeysize += switchingKeypacked[i].save(data_stream);
    }
    cout << "Detection Key Size: " << reskeysize << " bytes" << endl;
}


//////////////////////////////////////////////////// For Oblivious Multiplexer ////////////////////////////////////////////////////

// Pick random Zq elements as ID of recipients, in form of a (partySize x idSize) matrix.
vector<vector<int>> initializeRecipientId(const PVWParam& params, int partySize, int idSize) {
    vector<vector<int>> ids(partySize, vector<int> (idSize, -1)); 

    lbcrypto::DiscreteUniformGeneratorImpl<regevSK> dug;
    dug.SetModulus(params.q);

    for (int i = 0; i < (int)ids.size(); i++) {
        NativeVector temp = dug.GenerateVector(idSize);
        for (int j = 0; j < (int)ids[0].size(); j++) {
            ids[i][j] = temp[j].ConvertToInt();
        }
    }

    return ids;
}

/**
 * @brief Serve as sanity check, does not deal with GOMR logic; verify that the pertinent message clue can be accurately recovered
 * from the cluePoly based on the recipient's target ID.
 * 
 * @param params PVW parameter
 * @param extended_id exponential-extended version of a recipient ID, size = party_size_glb * id_size_glb
 * @param index pertinent message index
 * @param partySize party size
 * @param prepare if the encrypted version of ID is sent, prepare = true, and we minus q/4 beforehand for b since message = 1
 * @return true/false
 */
bool verify(unsigned char* expected_key, const PVWParam& params, const vector<int>& id, int index, int partySize = party_size_glb, bool prepare = false) {
    vector<uint64_t> polyFlat = loadDataSingle(index, "cluePoly", (params.n + params.ell) * (partySize + secure_extra_length_glb));
    unsigned char* key = (unsigned char *) malloc(sizeof(unsigned char) * AES_KEY_SIZE);
    loadSingleAESKey(key, index);

    vector<vector<int>> ids(1);
    ids[0] = id;
    vector<vector<int>> compressed_id = compressVectorByAES(params, key, ids, party_size_glb + secure_extra_length_glb);

    vector<vector<long>> cluePolynomial(params.n + params.ell, vector<long>(compressed_id[0].size()));
    vector<long> res(params.n + params.ell, 0);

    for (int i = 0; i < params.n + params.ell; i++) {
      for(int j = 0; j < (int)compressed_id[0].size(); j++) {
            res[i] = (res[i] + polyFlat[i * compressed_id[0].size() + j] * compressed_id[0][j]) % params.q;
            res[i] = res[i] < 0 ? res[i] + params.q : res[i];
        }
    }

    vector<uint64_t> expected = loadDataSingle(index, "clues", params.n + params.ell);

    for (int i = 0; i < params.n + params.ell; i++) {
        long temp = expected[i] - 16384;
        temp = temp < 0 ? temp + params.q : temp % params.q;
        if ((prepare && i >= params.n && temp != res[i]) || (prepare && i < params.n && expected[i] - res[i]) != 0 || (!prepare && expected[i] - res[i] != 0)) {
            return false;
        }
    }

    free(key);
    return true;
}

// similar to preparingTransactionsFormal but for gOMR with Oblivious Multiplexer.
vector<vector<uint64_t>> preparingGroupCluePolynomial(vector<int>& pertinentMsgIndices, PVWpk& pk, int numOfTransactions,int pertinentMsgNum,
                                  const PVWParam& params, const vector<int>& targetId, bool prepare = false, int clueLength = 454,
                                  int partySize = party_size_glb) {
    vector<vector<int>> ids;
    vector<PVWCiphertext> clues(party_size_glb);
    bool check = false;
    vector<int> zeros(params.ell, 0);
    vector<vector<uint64_t>> ret;

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    choosePertinentMsg(numOfTransactions, pertinentMsgNum, pertinentMsgIndices, seed);

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_time = 0;

    unsigned char* key = (unsigned char *) malloc(sizeof(unsigned char) * AES_KEY_SIZE);

    for(int i = 0; i < numOfTransactions; i++) {
        vector<PVWsk> impert_sk(party_size_glb);
        vector<PVWpk> impert_pk(party_size_glb-1);
        for (int p = 0; p < party_size_glb; p++) {
            impert_sk[p] = PVWGenerateSecretKey(params);
        }

        while (true) {
            random_bytes(key, AES_KEY_SIZE);
            
            if (find(pertinentMsgIndices.begin(), pertinentMsgIndices.end(), i) != pertinentMsgIndices.end()) {
                ret.push_back(loadDataSingle(i));
                check = true;
                ids = initializeRecipientId(params, party_size_glb - 1, id_size_glb);
                ids.push_back(targetId);
                for (int p = 0; p < party_size_glb-1; p++) {
                    impert_pk[p] = PVWGeneratePublicKey(params, impert_sk[p]);
                }
            } else {
                ids = initializeRecipientId(params, party_size_glb, id_size_glb);
            }

            time_start = chrono::high_resolution_clock::now();
            if (check) {
                for (int p = 0; p < party_size_glb - 1; p++) {
                    PVWEncPK(clues[p], zeros, impert_pk[p], params);
                }
                PVWEncPK(clues[party_size_glb-1], zeros, pk, params);
                saveClues(clues[party_size_glb-1], i);
            } else {
                for (int p = 0; p < party_size_glb; p++) {
                    PVWEncSK(clues[p], zeros, impert_sk[p], params);
                }
            }

            vector<vector<int>> compressed_ids = compressVectorByAES(params, key, ids, party_size_glb + secure_extra_length_glb);

            vector<vector<long>> cluePolynomial = agomr::generateClue(params, clues, compressed_ids, prepare);
            saveGroupClues(cluePolynomial, key, i);
            time_end = chrono::high_resolution_clock::now();

            if (check) {
                check = false;
                if (verify(key, params, ids[partySize-1], i, partySize, prepare)) {
                    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
                    break;
                } else {
//                     cout << "Mismatch detected, regenerating clue poly for msg: " << i << endl;
                }
            } else {
                break;
            }
        }
    }
    free(key);
    cout << "\nSender average running time: " << total_time / pertinentMsgIndices.size() << "us." << "\n";

    return ret;
}

void verify_fg(const PVWParam& params, const PVWsk& target_secretSK, const mre::MREsharedSK& target_sharedSK) {
    vector<vector<int>> ct = loadFixedGroupClues(0, 1, params);

    vector<int> targetCT = ct[0];
    
    for (int l = 0; l < params.ell; l++) {
        long result = 0;
        for (int i = 0; i < params.n; i++) {
            result = ( result + targetCT[i] * target_secretSK[l][i].ConvertToInt() ) % params.q;
            result = result < 0 ? result + params.q : result;
        }
        for (int i = 0; i < partial_size_glb; i++) {
            result = ( result + targetCT[params.n + l * partial_size_glb + i] * target_sharedSK[i].ConvertToInt()) % params.q;
            result = result < 0 ? result + params.q : result;
        }
    }
}

// similar to preparingTransactionsFormal but for fixed group GOMR which requires a MREGroupPK for each message.
// pertinentMsgIndices, groupPK, numOfTransactions, num_of_pertinent_msgs_glb, params, mreseed);
vector<vector<uint64_t>> preparingMREGroupClue(vector<int>& pertinentMsgIndices, int numOfTransactions, int pertinentMsgNum, const PVWParam& params,
                                               const PVWsk& target_secretSK, const mre::MREsharedSK& target_sharedSK, prng_seed_type& seed,
                                               const int partialSize = partial_size_glb, const int partySize = party_size_glb) {

    vector<vector<uint64_t>> ret;
    vector<int> zeros(params.ell, 0);
    PVWsk sk;
    vector<fgomr::FixedGroupSecretKey> groupSK;
    fgomr::FixedGroupSharedKey gPK;

    choosePertinentMsg(numOfTransactions, pertinentMsgNum, pertinentMsgIndices, seed);

    prng_seed_type mreseed, compress_seed;
    for (auto &i : mreseed) { // the seed to randomly sample A1, and b in secret key
        i = random_uint64();
    }

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_time = 0;

    for(int i = 0; i < numOfTransactions; i++){
        PVWCiphertext tempclue;

        for (auto &i : compress_seed) { // the seed to perform exponential extension for the sharedSK
            i = random_uint64();
        }
        if (find(pertinentMsgIndices.begin(), pertinentMsgIndices.end(), i) != pertinentMsgIndices.end()) {
            groupSK = fgomr::secretKeyGen(params, target_secretSK, target_sharedSK);
            gPK = fgomr::groupKeyGenAux(params, groupSK, mreseed);

            time_start = chrono::high_resolution_clock::now();
            tempclue = fgomr::genClue(params, zeros, gPK, compress_seed);
            time_end = chrono::high_resolution_clock::now();
            total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

            ret.push_back(loadDataSingle(i));
            saveCluesWithRandomness(tempclue, i, compress_seed);
            verify_fg(params, target_secretSK, target_sharedSK);
        } else {
            auto non_pert_params = PVWParam(params.n + (partySize + secure_extra_length_glb) * params.ell, params.q, params.std_dev, params.m, params.ell);
            sk = PVWGenerateSecretKey(non_pert_params);
            PVWEncSK(tempclue, zeros, sk, non_pert_params);
            saveCluesWithRandomness(tempclue, i, compress_seed);
        }  
    }

    cout << "\nSender average running time: " << total_time / pertinentMsgIndices.size() << "us." << "\n";

    return ret;
}
