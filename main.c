#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "temp.h"

/* PREPROCESSOR BEGIN                   */

#define N_SYMBOLS   4
#define DATA_SIZE   28
/* PREPROCESSOR END                     */


/* LOCAL OBJECTS BEGIN                  */

const static uint8_t poly_generator[N_SYMBOLS + 1] = { 0x1, 0xf, 0x36, 0x78, 0x40 }; /* Precomputed generator polynomial for 4 error correction symbols */
const static uint16_t poly_prime = 0x11d;
/* LOCAL OBJECTS END                    */


/* LOCAL FUNCTIONS DECLARATIONS BEGIN   */

inline static uint8_t gf_pow2(uint8_t pow);

static uint8_t  gf_mult (uint8_t x, uint8_t y);  
static uint8_t  gf_poly_evaluate (const uint8_t *polynomial, uint8_t size, uint8_t x);
static void     gf_poly_scale (uint8_t *new_poly, const uint8_t *old_poly, uint8_t size, uint8_t scale);

static void     compute_poly_syndromes (uint8_t *poly_syndromes, const uint8_t *encoded_message);
static uint8_t  compute_poly_err_locator (uint8_t *poly_err_locator, const uint8_t *poly_syndromes);
static uint8_t  compute_poly_err_evaluator (uint8_t *poly_err_evaluator, const uint8_t *poly_err_locator);

static void     print_hex_n (const uint8_t *str, uint8_t length);
/* LOCAL FUNCTIONS DECLARATIONS END     */


int main(int argc, char *argv[])
{
    uint8_t message_in[DATA_SIZE] = { 0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
                                             0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec,
                                             0x37, 0x17, 0x17, 0x73, 0x12, 0x91, 0x37, 0xab,
                                             0x1b, 0x3d, 0xd7, 0xe2 };
    uint8_t encoded_message[DATA_SIZE + N_SYMBOLS] = { 0 };
    uint8_t coefficient;


    /**********************\
     *      ENCODING      *
    \**********************/

    memcpy(encoded_message, message_in, DATA_SIZE);
   
    for (uint8_t i = 0; i < DATA_SIZE; ++i)
    {
        coefficient = encoded_message[i];

        if (coefficient != 0)
        {
            for (uint8_t j = 1; j < N_SYMBOLS + 1; ++j)
            {
                encoded_message[i+j] ^= gf_mult(poly_generator[j], coefficient);
            }
        }
    }

    memcpy(encoded_message, message_in, DATA_SIZE);

    printf("Encoded 28-byte data to 32-byte data with last 4 bytes being error correction symbols: \n\n");
    print_hex_n(encoded_message, DATA_SIZE + N_SYMBOLS);
    printf("\n");
 

    /**********************\
     *      DECODING      *
    \**********************/
    ++encoded_message[0];
    ++encoded_message[1];
    //++encoded_message[2];
    uint8_t poly_syndromes[N_SYMBOLS];
    compute_poly_syndromes(poly_syndromes, encoded_message);

    printf("Computed polynomial syndromes for encoded message: \n\n");
    print_hex_n(poly_syndromes, N_SYMBOLS);
    printf("\n");

    uint8_t poly_err_locator[N_SYMBOLS];
    uint8_t poly_err_locator_size = compute_poly_err_locator(poly_err_locator, poly_syndromes);

    printf("Computed polynomial error locator with size %d for encoded message: \n\n", poly_err_locator_size);
    print_hex_n(poly_err_locator, N_SYMBOLS);

    /* Invert order of coefficients in poly_err_locator */
    for (uint8_t i = 0; i < N_SYMBOLS / 2; ++i)
    {
        poly_err_locator[i] ^= poly_err_locator[N_SYMBOLS - i - 1];
        poly_err_locator[N_SYMBOLS - i - 1] ^= poly_err_locator[i];
        poly_err_locator[i] ^= poly_err_locator[N_SYMBOLS - i - 1];
    }
    printf("\n");
    printf("Computed polynomial error locator with size %d for encoded message: \n\n", poly_err_locator_size);
    print_hex_n(poly_err_locator, N_SYMBOLS);

    uint8_t poly_err_evaluator[N_SYMBOLS];
    uint8_t poly_err_evaluator_size = compute_poly_err_evaluator(poly_err_evaluator, poly_err_locator);

    printf("Computed polynomial error evaluator with size %d for encoded message: \n\n", poly_err_evaluator_size);
    print_hex_n(poly_err_evaluator, N_SYMBOLS);
    printf("\n");

    return 0;
}





/* LOCAL FUNCTIONS DEFINITIONS BEGIN    */

/* Function for multiplying two numbers in GF(2^8) modulo poly_prime */
static uint8_t gf_mult(uint8_t x, uint8_t y)
{
    uint8_t result = 0;
    uint16_t temp_x = x;    /* uint16_t because we need to check if the value passed 256 (2^8) to divide it modulo poly_prime*/

    while (y > 0)
    {
        if (y & 0x01)
        {
            result ^= temp_x;
        }

        y >>= 1;
        temp_x <<= 1;

        if (temp_x & 0x0100)
        {
            temp_x ^= poly_prime;
        }
    }

    return result;
}

inline static uint8_t gf_pow2(uint8_t pow)
{
    uint16_t result = 1;
    
    while (pow-- > 0)
    {
        result <<= 1;
        if (result & 0x0100)
        {
            result ^= poly_prime;
        }
    }

    return (uint8_t)result;
}

/* Function for evaluating value of polynomial at given x based on Horner's scheme */
static uint8_t gf_poly_evaluate(const uint8_t *polynomial, uint8_t size, uint8_t x)
{
    uint8_t y = polynomial[0];
    for (uint8_t i = 1; i < size; ++i)
    {
        y = gf_mult(y, x) ^ polynomial[i];
    }
    return y;
}

/* Function for scaling polynomial coefficients */
static void gf_poly_scale(uint8_t *new_poly, const uint8_t *old_poly, uint8_t size, uint8_t scale)
{
    for (uint8_t i = 0; i < size; ++i)
    {
        new_poly[i] = gf_mult(old_poly[i], scale);
    }
}

/* Function for computing syndromes polynomial (treating message as polynomial and evaluating its value at points which should be its roots if the data is not damaged) */
static void compute_poly_syndromes(uint8_t *poly_syndromes, const uint8_t *encoded_message)
{
    for (uint8_t i = 0; i < N_SYMBOLS; ++i)
    {
        poly_syndromes[i] = gf_poly_evaluate(encoded_message, DATA_SIZE + N_SYMBOLS, 1 << i);       /* because our N_SYMBOLS is 4 we can treat 2^i in GF(2^8) as simple 1 << i */
    }
}

/* Function for computing error locator polynomial (Berlekamp-Massey algorithm) */
static uint8_t compute_poly_err_locator(uint8_t *poly_err_locator, const uint8_t *poly_syndromes)
{
    uint8_t curr_poly_err_locator[N_SYMBOLS] = { 1, 0, 0, 0 };
    uint8_t old_poly_err_locator[N_SYMBOLS]  = { 1, 0, 0, 0 };
    uint8_t curr_size = 1;
    uint8_t old_size = 1;
    uint8_t discrepancy = 0;
    //L = 0 # update flag variable, not needed here because we use an alternative equivalent way of checking if update is needed (but using the flag could potentially be faster depending on if using length(list) is taking linear time in your language, here in Python it's constant so it's as fast.

    for (uint8_t i = 0; i < N_SYMBOLS; ++i)
    {
        printf("\ni=%d\n", i);

        // K = i
        discrepancy = poly_syndromes[i];

        printf("discrepancy: %d\n", discrepancy);
        for (uint8_t j = 1; j < curr_size; ++j)
        {
            discrepancy ^= gf_mult(curr_poly_err_locator[curr_size-j-1], poly_syndromes[i-j]);
        }
        printf("discrepancy: %d\n", discrepancy);

        ++old_size;

        if (discrepancy != 0)
        {
            uint8_t temp_poly_err_locator[N_SYMBOLS];

            printf("old_size=%d\ncurr_size=%d\n", old_size, curr_size);


            if (old_size > curr_size)
            {
                gf_poly_scale(temp_poly_err_locator, old_poly_err_locator, N_SYMBOLS, discrepancy);
                gf_poly_scale(old_poly_err_locator, curr_poly_err_locator, N_SYMBOLS, gf_inv(discrepancy));

                /* swap old_size and curr_size */
                old_size ^= curr_size;
                curr_size ^= old_size;
                old_size ^= curr_size;

                memcpy(curr_poly_err_locator, temp_poly_err_locator, N_SYMBOLS);

                printf("\ntemp_loc:\n");
                print_hex_n(temp_poly_err_locator, N_SYMBOLS);
                printf("\n");
                printf("\nold_loc:\n");
                print_hex_n(old_poly_err_locator, N_SYMBOLS);
                printf("\n");
                printf("\ncurr_loc:\n");
                print_hex_n(curr_poly_err_locator, N_SYMBOLS);
                printf("\n");

            }

            gf_poly_scale(temp_poly_err_locator, old_poly_err_locator, N_SYMBOLS, discrepancy);

            /* old_size <= curr_size */
            printf("old_size=%d\ncurr_size=%d\n", old_size, curr_size);

            uint8_t size_shift = curr_size - old_size;

            for (uint8_t iter = size_shift; iter < N_SYMBOLS; ++iter)
            {
                curr_poly_err_locator[iter] ^= temp_poly_err_locator[iter - size_shift];
            }

            printf("\ncurr_loc:\n");
            print_hex_n(curr_poly_err_locator, N_SYMBOLS);
            printf("\n");

        }
    }

    memcpy(poly_err_locator, curr_poly_err_locator, N_SYMBOLS);
    return curr_size;
}


static uint8_t compute_poly_err_evaluator (uint8_t *poly_err_evaluator, const uint8_t *poly_err_locator)
{
    uint8_t error_counter = 0;
    /* Brute force checking */
    for (uint8_t i = 0; i < DATA_SIZE + N_SYMBOLS; ++i)
    {
        /* Check every byte for condition to be root of polynomial error locator (the place where the erroc has happened) */
        uint8_t poly_value = gf_poly_evaluate(poly_err_locator, N_SYMBOLS, gf_pow2(i));
        printf("Poly value for gf_pow2(i)<i> [%d<%d>] is %d\n", gf_pow2(i), i, poly_value);
        if (poly_value == 0)
        {
            printf("Found error\n");
            if (2*error_counter >= N_SYMBOLS)
            {
                /* Found too many errors */
                break;
            }
            poly_err_evaluator[error_counter++] = DATA_SIZE + N_SYMBOLS - 1 - i;
        }
    }

    return error_counter;
}

/* Function for displaying uint8_t arrays in hex */
static void print_hex_n (const uint8_t *str, uint8_t length)
{
    while (length-- > 0)
    {
        printf("%02x ", *str++);
    }
}


/* LOCAL FUNCTIONS DEFINITIONS END      */
   
