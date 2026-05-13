#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>

/**
 * Пространство с тестами NIST STS
 */
namespace tests {
    /**
     * @brief Runs Test из NIST STS
     */
    void Runs(const std::vector<uint32_t>& sample) {
        int     S, k;
        double  pi, V, erfc_arg, p_value;
        
        int n = sample.size() * 32;
        const double ALPHA = 0.01;

        S = 0;
        for ( k=0; k<n; k++ )
            if ( ((sample[k/32] >> (k%32)) & 1) )
                S++;
        pi = (double)S / (double)n;

        if ( std::fabs(pi - 0.5) > (2.0 / std::sqrt(n)) ) {
            printf("\t\t\t\tRUNS TEST\n");
            printf("\t\t------------------------------------------\n");
            printf("\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
            p_value = 0.0;
            printf("FAILURE\t\tp_value = %f\n\n", p_value);
        }
        else {

            V = 1;
            for ( k=1; k<n; k++ ) {
                int current_bit = (sample[k/32] >> (k%32)) & 1;
                int prev_bit = (sample[(k-1)/32] >> ((k-1)%32)) & 1;
                if ( current_bit != prev_bit )
                    V++;
            }
        
            erfc_arg = std::fabs(V - 2.0 * n * pi * (1.0 - pi)) / (2.0 * pi * (1.0 - pi) * std::sqrt(2.0 * n));
            
            p_value = std::erfc(erfc_arg);
            
            printf("\t\t\t\tRUNS TEST\n");
            printf("\t\t------------------------------------------\n");
            printf("\t\tCOMPUTATIONAL INFORMATION:\n");
            printf("\t\t------------------------------------------\n");
            printf("\t\t(a) Pi                        = %f\n", pi);
            printf("\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
            printf("\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
            printf("\t\t    -----------------------   = %f\n", erfc_arg);
            printf("\t\t      2 sqrt(2n) pi (1-pi)\n");
            printf("\t\t------------------------------------------\n");
            
            if ( p_value < 0.0 || p_value > 1.0 )
                printf("WARNING:  P_VALUE IS OUT OF RANGE.\n");

            printf("%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
        }
    }

    /**
     * Longest Run Of Ones Test из NIST STS
     */
    void LongestRunOfOnes(const std::vector<uint32_t>& sample) {
        double          chi2, pi[7];
        int             run, v_n_obs, N, i, j, K, M, V[7];
        unsigned int    nu[7] = { 0, 0, 0, 0, 0, 0, 0 };
        
        int n = sample.size() * 32;

        if ( n < 128 ) {
            printf("\t\t\t  LONGEST RUNS OF ONES TEST\n");
            printf("\t\t---------------------------------------------\n");
            printf("\t\t   n=%d is too short\n", n);
            return;
        }
        if ( n < 6272 ) {
            K = 3;
            M = 8;
            V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
            pi[0] = 0.21484375;
            pi[1] = 0.3671875;
            pi[2] = 0.23046875;
            pi[3] = 0.1875;
        }
        else if ( n < 750000 ) {
            K = 5;
            M = 128;
            V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
            pi[0] = 0.1174035788;
            pi[1] = 0.242955959;
            pi[2] = 0.249363483;
            pi[3] = 0.17517706;
            pi[4] = 0.102701071;
            pi[5] = 0.112398847;
        }
        else {
            K = 6;
            M = 10000;
                V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
            pi[0] = 0.0882;
            pi[1] = 0.2092;
            pi[2] = 0.2483;
            pi[3] = 0.1933;
            pi[4] = 0.1208;
            pi[5] = 0.0675;
            pi[6] = 0.0727;
        }
        
        N = n/M;
        for ( i=0; i<N; i++ ) {
            v_n_obs = 0;
            run = 0;
            for ( j=0; j<M; j++ ) {
                int global_idx = i*M+j;
                if ( ((sample[global_idx/32] >> (global_idx%32)) & 1) == 1 ) {
                    run++;
                    if ( run > v_n_obs )
                        v_n_obs = run;
                }
                else
                    run = 0;
            }
            if ( v_n_obs < V[0] )
                nu[0]++;
            for ( j=0; j<=K; j++ ) {
                if ( v_n_obs == V[j] )
                    nu[j]++;
            }
            if ( v_n_obs > V[K] )
                nu[K]++;
        }

        chi2 = 0.0;
        for ( i=0; i<=K; i++ )
            chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

        double chi2_crit = (K == 3) ? 11.345 : ((K == 5) ? 15.086 : 16.812);
        int is_success = (chi2 <= chi2_crit);

        printf("\t\t\t  LONGEST RUNS OF ONES TEST\n");
        printf("\t\t---------------------------------------------\n");
        printf("\t\tCOMPUTATIONAL INFORMATION:\n");
        printf("\t\t---------------------------------------------\n");
        printf("\t\t(a) N (# of substrings)  = %d\n", N);
        printf("\t\t(b) M (Substring Length) = %d\n", M);
        printf("\t\t(c) Chi^2                = %f\n", chi2);
        printf("\t\t---------------------------------------------\n");
        printf("\t\t      F R E Q U E N C Y\n");
        printf("\t\t---------------------------------------------\n");

        if ( K == 3 ) {
            printf("\t\t  <=1     2     3    >=4   Chi^2      Assignment");
            printf("\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
        }
        else if ( K == 5 ) {
            printf("\t\t<=4  5  6  7  8  >=9 Chi^2      Assignment");
            printf("\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5]);
        }
        else {
            printf("\t\t<=10  11  12  13  14  15 >=16 Chi^2      Assignment");
            printf("\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], nu[6]);
        }

        printf("\n%s\t\tChi^2 = %f (Crit: %f)\n\n", is_success ? "SUCCESS" : "FAILURE", chi2, chi2_crit);
    }

    double psi2(int m, int n, const std::vector<uint32_t>& sample) {
        if (m == 0 || m == -1)
            return 0.0;
            
        double numOfBlocks = n;
        int powLen = (int)std::pow(2, m + 1) - 1;
        
        std::vector<unsigned int> P(powLen, 0);

        for (int i = 0; i < numOfBlocks; i++) {
            int k = 1;
            for (int j = 0; j < m; j++) {
                int global_idx = (i + j) % n;
                int bit = (sample[global_idx / 32] >> (global_idx % 32)) & 1;
                
                if (bit == 0)
                    k *= 2;
                else if (bit == 1)
                    k = 2 * k + 1;
            }
            P[k - 1]++;
        }
        
        double sum = 0.0;
        for (int i = (int)std::pow(2, m) - 1; i < (int)std::pow(2, m + 1) - 1; i++) {
            sum += std::pow(P[i], 2);
        }
        sum = (sum * std::pow(2, m) / (double)n) - (double)n;
        
        return sum;
    }

    /**
     * @brief Serial Test из NIST STS
     */
    void Serial(const std::vector<uint32_t>& sample, int m = 5) {
        int n = sample.size() * 32;
        double psim0, psim1, psim2, del1, del2;
        
        psim0 = psi2(m, n, sample);
        psim1 = psi2(m - 1, n, sample);
        psim2 = psi2(m - 2, n, sample);
        
        del1 = psim0 - psim1;
        del2 = psim0 - 2.0 * psim1 + psim2;
        
        double crit1 = 32.000; // df = 16, alpha = 0.01
        double crit2 = 20.090; // df = 8, alpha = 0.01

        const double ALPHA = 0.01;

        std::cout << "  [NIST] SERIAL TEST\n";
        std::cout << "    ---------------------------------------------\n";
        std::cout << "    (a) Block length    (m) = " << m << "\n";
        std::cout << "    (b) Sequence length (n) = " << n << "\n";
        std::cout << "    (c) Psi_m               = " << psim0 << "\n";
        std::cout << "    (d) Psi_m-1             = " << psim1 << "\n";
        std::cout << "    (e) Psi_m-2             = " << psim2 << "\n";
        std::cout << "    (f) Del_1               = " << del1 << "\n";
        std::cout << "    (g) Del_2               = " << del2 << "\n";
        std::cout << "    ---------------------------------------------\n";

        if (del1 > crit1) {
            std::cout << "    FAILURE\t\tdel1_stat = " << del1 << " (Crit: " << crit1 << ")\n";
        } else {
            std::cout << "    SUCCESS\t\tdel1_stat = " << del1 << "\n";
        }

        if (del2 > crit2) {
            std::cout << "    FAILURE\t\tdel2_stat = " << del2 << " (Crit: " << crit2 << ")\n\n";
        } else {
            std::cout << "    SUCCESS\t\tdel2_stat = " << del2 << "\n\n";
        }
    }

    /**
     * @brief Approximate Entropy Test из NIST STS
     */
    void ApproximateEntropy(const std::vector<uint32_t>& sample, int m = 5) {
        int i, j, k, r, blockSize, seqLength, powLen, index;
        double sum, numOfBlocks, ApEn[2], apen, chi_squared;
        
        seqLength = sample.size() * 32;

        std::cout << "  [NIST] APPROXIMATE ENTROPY TEST\n";
        std::cout << "    --------------------------------------------\n";
        std::cout << "    COMPUTATIONAL INFORMATION:\n";
        std::cout << "    --------------------------------------------\n";
        std::cout << "    (a) m (block length)    = " << m << "\n";

        r = 0;
        
        for (blockSize = m; blockSize <= m + 1; blockSize++) {
            if (blockSize == 0) {
                ApEn[0] = 0.00;
                r++;
            } else {
                numOfBlocks = (double)seqLength;
                powLen = (int)std::pow(2, blockSize + 1) - 1;
                
                std::vector<unsigned int> P(powLen, 0);
                
                for (i = 0; i < numOfBlocks; i++) {
                    k = 1;
                    for (j = 0; j < blockSize; j++) {
                        k <<= 1;
                        
                        int global_idx = (i + j) % seqLength;
                        int bit = (sample[global_idx / 32] >> (global_idx % 32)) & 1;
                        
                        if (bit == 1)
                            k++;
                    }
                    P[k - 1]++;
                }
                
                sum = 0.0;
                index = (int)std::pow(2, blockSize) - 1;
                for (i = 0; i < (int)std::pow(2, blockSize); i++) {
                    if (P[index] > 0)
                        sum += P[index] * std::log(P[index] / numOfBlocks);
                    index++;
                }
                sum /= numOfBlocks;
                ApEn[r] = sum;
                r++;
            }
        }
        
        apen = ApEn[0] - ApEn[1];
        chi_squared = 2.0 * seqLength * (std::log(2.0) - apen);
        
        double chi2_crit = 53.486; 

        std::cout << "    (b) n (sequence length) = " << seqLength << "\n";
        std::cout << "    (c) Chi^2               = " << chi_squared << "\n";
        std::cout << "    (d) Phi(m)              = " << ApEn[0] << "\n";
        std::cout << "    (e) Phi(m+1)            = " << ApEn[1] << "\n";
        std::cout << "    (f) ApEn                = " << apen << "\n";
        std::cout << "    (g) Log(2)              = " << std::log(2.0) << "\n";
        std::cout << "    --------------------------------------------\n";

        if (m > (int)(std::log(seqLength) / std::log(2.0) - 5)) {
            int recommended = std::max(1, (int)(std::log(seqLength) / std::log(2.0) - 5));
            std::cout << "    Note: The blockSize = " << m << " exceeds recommended value of " << recommended << "\n";
            std::cout << "    Results are inaccurate!\n";
            std::cout << "    --------------------------------------------\n";
        }
        
        if (chi_squared > chi2_crit) {
            std::cout << "    FAILURE\t\tchi2_stat = " << chi_squared << " (Crit: " << chi2_crit << ")\n\n";
        } else {
            std::cout << "    SUCCESS\t\tchi2_stat = " << chi_squared << "\n\n";
        }
    }

    /**
     * @brief Frequency Test из NIST STS
     */
    void Frequency(const std::vector<uint32_t>& sample) {
        int     i;
        double  f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
        
        int n = sample.size() * 32;
        const double ALPHA = 0.01;
        
        sum = 0.0;
        for ( i=0; i<n; i++ ) {
            int bit = (sample[i/32] >> (i%32)) & 1;
            sum += 2 * bit - 1;
        }
        
        s_obs = std::fabs(sum)/std::sqrt(n);
        f = s_obs/sqrt2;
        p_value = std::erfc(f);

        printf("\t\t\t      FREQUENCY TEST\n");
        printf("\t\t---------------------------------------------\n");
        printf("\t\tCOMPUTATIONAL INFORMATION:\n");
        printf("\t\t---------------------------------------------\n");
        printf("\t\t(a) The nth partial sum = %d\n", (int)sum);
        printf("\t\t(b) S_n/n               = %f\n", sum/n);
        printf("\t\t---------------------------------------------\n");

        printf("%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
    }
    /* ------------------- */
    void runAllNistTests(const std::vector<uint32_t>& sample) {
        Frequency(sample);
        Runs(sample);
        LongestRunOfOnes(sample);
        Serial(sample);
        ApproximateEntropy(sample);
    }
}