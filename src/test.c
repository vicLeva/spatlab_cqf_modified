/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"

#include <stdint.h>

void print_bits(uint64_t value) {
    for (int i = 63; i >= 0; i--) {
        printf("%d", (int)((value >> i) & 1));
    }
    printf("\n");
}


void print_mer(const char *kmer){
    printf("32-mer : ");
    
    for (int i = 0; kmer[i] != '\0'; i++) {
        printf("%c", kmer[i]);
    }
    
    printf("\n");
}

const uint8_t rev_table[256] = {
  0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f, 0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0xf, 
  0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b, 0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0xb, 
  0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27, 0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x7, 
  0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23, 0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x3, 
  0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e, 0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0xe, 
  0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a, 0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0xa, 
  0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26, 0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x6, 
  0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22, 0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x2, 
  0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d, 0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0xd, 
  0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29, 0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x9, 
  0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25, 0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x5, 
  0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21, 0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x1, 
  0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c, 0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0xc, 
  0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28, 0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x8, 
  0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24, 0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x4, 
  0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20, 0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x0};

uint64_t revcomp64 (const uint64_t v, size_t bitsize){
  return (((uint64_t)rev_table[v & 0xff] << 56) | 
    ((uint64_t)rev_table[(v >> 8) & 0xff] << 48) | 
    ((uint64_t)rev_table[(v >> 16) & 0xff] << 40) | 
    ((uint64_t)rev_table[(v >> 24) & 0xff] << 32) | 
    ((uint64_t)rev_table[(v >> 32) & 0xff] << 24) | 
    ((uint64_t)rev_table[(v >> 40) & 0xff] << 16) |
    ((uint64_t)rev_table[(v >> 48) & 0xff] << 8) |
    ((uint64_t)rev_table[(v >> 56) & 0xff])) >> (64-bitsize);
}

uint64_t canonical(uint64_t smer, size_t size){
    uint64_t revcomp = revcomp64(smer, size);
    if (revcomp < smer) { return revcomp; }
    else { return smer; }
}

uint64_t encode_for_canonical(const char* kmer) {
    uint64_t encoded = 0;
    
    for (int i = 0; kmer[i] != '\0'; i++) {
        char c = kmer[i];
        
        if (c == 'A') {
            encoded <<= 2;
        } else if (c == 'C') {
            encoded <<= 2;
            encoded |= 1;
        } else if (c == 'G') {
            encoded <<= 2;
            encoded |= 2;
        } else {
            encoded <<= 2;
            encoded |= 3;
        }
    }
    
    return encoded;
}

uint64_t encode(const char* kmer) {
    uint64_t encoded = 0;
    
    for (int i = 0; kmer[i] != '\0'; i++) {
        char c = kmer[i];
        
        if (c == 'A') {
            encoded <<= 2;
            encoded |= 3;
        } else if (c == 'C') {
            encoded <<= 2;
            encoded |= 2;
        } else if (c == 'G') {
            encoded <<= 2;
            encoded |= 1;
        } else {
            encoded <<= 2;
        }
    }
    
    return encoded;
}


uint64_t bfc_hash_64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}



uint64_t kmer_to_hash(const char* kmer, uint64_t k) {
    uint64_t encoded = encode(kmer);
    uint64_t mask = (1ULL << (k * 2)) - 1;
    return bfc_hash_64(encoded, mask);
}

uint64_t canonical_code_to_hash(uint64_t canonical, uint64_t k) {
    uint64_t mask = (1ULL << (k * 2)) - 1;
    return bfc_hash_64(~canonical, mask);
}



char generate_random_base() {
    char bases[] = {'A', 'C', 'G', 'T'};
    int index = rand() % 4;
    return bases[index];
}

void generate_random_kmer(char* kmer) {
    for (int i = 0; i < 32; i++) {
        kmer[i] = generate_random_base();
    }
    kmer[32] = '\0';  // Null-terminate the string
}


int main(int argc, char **argv)
{
	if (argc < 3) {
		fprintf(stderr, "Please specify the log of the number of slots and the number of remainder bits in the CQF.\n");
		exit(1);
	}
	QF qf;
	uint64_t qbits = atoi(argv[1]);
	uint64_t rbits = atoi(argv[2]);
	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);

	char posQueries[100000][33];
	char negQueries[100000][33];

	clock_t start;
	clock_t end;
	double elapsed_time;

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    // BUILD

	start = clock();  // Start the timer
	qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_NONE, 0);


	

    uint64_t i = 0;

	const char* kmer;
	uint64_t count;

    fp = fopen("/scratch/vlevallois/data/AHX_ACXIOSF_6_1_32_2andmore.txt", "r");

    while ((read = getline(&line, &len, fp)) != -1) {
		kmer = strtok(line, "\t");
		count = strtoull(strtok(NULL, "\t"), NULL, 0);
		qf_insert(&qf, kmer_to_hash(kmer, 32), 0, count, QF_NO_LOCK);
        i++;
        if (i % 500000 == 0){
            printf("%lu / 514 -> %lu\n", i, qf_get_num_occupied_slots(&qf));
        }
        
        //printf("%s->%d\n", kmer, count);
    }
	end = clock();  // Stop the timer


	elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time build + inserts: %.6f seconds\n", elapsed_time);

    fclose(fp);

	printf("%" PRIu64 "(at the end)\n", qf_get_num_occupied_slots(&qf));

    // QUERIES

    start = clock();  // Start the timer

    fp = fopen("/home/genouest/genscale/vlevallois/data/queries.fasta", "r"); 
    //this file is basically reads from the original fastq, used to build the index 

    uint64_t query;

    while ((read = getline(&line, &len, fp)) != -1) {
        line[strlen(line) - 1] = '\0';
        //printf("%s\n", line);
        //printf("%d\n", strlen(line));

        for (int i = 0; i <= strlen(line) - 32; i++) {
            char kmer[33];
            strncpy(kmer, line + i, 32);
            kmer[32] = '\0';
            query = qf_count_key_value(&qf, canonical_code_to_hash(canonical(encode_for_canonical(kmer), 64), 32), 0, 0);
            //printf("32-mer : %s, %d\n", kmer, query);
        }
    }
	end = clock();  // Stop the timer

	elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time queries: %.6f seconds\n", elapsed_time);

    fclose(fp);

    start = clock();  // Start the timer

    fp = fopen("/home/genouest/genscale/vlevallois/data/neg_queries.fasta", "r");
    //this file is a file of randomly generated reads of length 70-110 nt


    while ((read = getline(&line, &len, fp)) != -1) {
        line[strlen(line) - 1] = '\0';
        //printf("%s\n", line);
        //printf("%d\n", strlen(line));

        for (int i = 0; i <= strlen(line) - 32; i++) {
            char kmer[33];
            strncpy(kmer, line + i, 32);
            kmer[32] = '\0';
            query = qf_count_key_value(&qf, canonical_code_to_hash(canonical(encode_for_canonical(kmer), 64), 32), 0, 0);
            //printf("32-mer : %s, %d\n", kmer, query);
        }
    }
	end = clock();  // Stop the timer

	elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time queries: %.6f seconds\n", elapsed_time);

    fclose(fp);



	/*printf("start veriff\n");

	FILE * fp2;

	int surestim = 0;
    int totsurestim = 0;
    int bug = 0;
    
    i = 0;
    uint64_t query;
    int count_limit = 0;

	

	 fp2 = fopen("/scratch/vlevallois/data/AHX_ACXIOSF_6_1_32_all.txt", "r");


	while ((read = getline(&line, &len, fp2)) != -1) {
		kmer = strtok(line, "\t");
		count = strtoull(strtok(NULL, "\t"), NULL, 0);

		if (kmer != NULL && i < 100000) {
            strncpy(posQueries[i], kmer, 32);
			posQueries[i][32] = '\0';
		}

		i++;
		if (i%10000000 == 0){
			printf("%f / 158\n", i/10000000.0);
		}
		

		query = qf_count_key_value(&qf, kmer_to_hash(kmer, 32), 0, 0);

		if (query > count){
            surestim ++;
            totsurestim = totsurestim + query - count;
        } 
        if (query < count){

            printf("%s kmer--count %d  vs query : %s\n", kmer, count, query);
            bug ++;
        } 
    }

    fclose(fp2);
    if (line)
        free(line);

    printf("nb elems inserted : %d \n", i);
    printf("nb surestim : %d \n", surestim);
    printf("avg surestim : %f \n", totsurestim / (double)surestim);
    printf("nb bug : %d \n", bug);
    printf("end verif\n"); */


	/* srand(time(NULL));  // Seed the random number generator
	for (int i = 0; i < 100000; i++) {
        generate_random_kmer(negQueries[i]);
    }

	
	start = clock();
    for (int k = 0; k < 100; k++) {
        for (int i = 0; i < 100000; i++) {
                count_limit = qf_count_key_value(&qf, kmer_to_hash(negQueries[i], 32), 0, 0);
            }
    }
	end = clock();  // Stop the timer

	elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time 10.000.000 neg queries: %.6f seconds\n", elapsed_time);


	start = clock();
	for (int i = 0; i < 100000; i++) {
        count_limit = qf_count_key_value(&qf, kmer_to_hash(posQueries[i], 32), 0, 0);
    }
	end = clock();  // Stop the timer

	elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time 100k pos queries: %.6f seconds\n", elapsed_time); */





	qf_free(&qf);

	fprintf(stdout, "Validated the CQF.\n");

	#include "/home/genouest/genscale/vlevallois/my_env/include/openssl/opensslconf.h"
}
