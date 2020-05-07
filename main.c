#include <stdio.h>
#include <stdint.h>
#include <string.h>

/* PREPROCESSOR BEGIN                   */

#define N_SYMBOLS   4
#define DATA_SIZE   28

#define ERR_CODE    0xff
/* PREPROCESSOR END                     */


/* LOCAL OBJECTS BEGIN                  */

const static uint8_t poly_generator[N_SYMBOLS + 1] = { 0x1, 0xf, 0x36, 0x78, 0x40 }; /* Precomputed generator polynomial for 4 error correction symbols */
const static uint16_t poly_prime = 0x11d;


static uint8_t poly_syndromes[N_SYMBOLS];
static uint8_t poly_err_locator[N_SYMBOLS];
static uint8_t poly_err_evaluator[N_SYMBOLS];
static uint8_t poly_corruption[N_SYMBOLS/2];

static uint8_t poly_err_locator_size;
static uint8_t poly_err_evaluator_size;

/* LOCAL OBJECTS END                    */


/* LOCAL FUNCTIONS DECLARATIONS BEGIN   */

static uint8_t  gf_pow2(uint8_t pow);
static uint8_t  gf_mult (uint8_t x, uint8_t y);  
static uint8_t  gf_inv (uint8_t x);
static uint8_t  gf_poly_evaluate (const uint8_t *polynomial, uint8_t size, uint8_t x);
static void     gf_poly_scale (uint8_t *new_poly, const uint8_t *old_poly, uint8_t size, uint8_t scale);

static void     compute_poly_syndromes (const uint8_t *encoded_message);
static void     compute_poly_err_locator ();
static void     compute_poly_err_evaluator ();
static uint8_t  compute_poly_corruption ();

static void     print_hex_n (const uint8_t *str, uint8_t length);
/* LOCAL FUNCTIONS DECLARATIONS END     */


int main(int argc, char *argv[])
{
    uint8_t message_in[DATA_SIZE] = { 0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
                                             0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec,
                                             0x37, 0x17, 0x17, 0x73, 0x12, 0x91, 0x37, 0xab,
                                             0x1b, 0x3d, 0xd7, 0xe2 };
    uint8_t encoded_message[DATA_SIZE + N_SYMBOLS] = { 0 };


    /**********************\
     *      ENCODING      *
    \**********************/

    memcpy(encoded_message, message_in, DATA_SIZE);
   
    for (uint8_t i = 0; i < DATA_SIZE; ++i)
    {
        uint8_t coefficient = encoded_message[i];

        if (coefficient != 0)
        {
            for (uint8_t j = 1; j < N_SYMBOLS + 1; ++j)
            {
                encoded_message[i+j] ^= gf_mult(poly_generator[j], coefficient);
            }
        }
    }

    memcpy(encoded_message, message_in, DATA_SIZE);

    printf("Encoded 28-byte data to 32-byte data with last 4 bytes being error correction symbols: \n");
    print_hex_n(encoded_message, DATA_SIZE + N_SYMBOLS);
    printf("\n");
 

    /**********************\
     *      DECODING      *
    \**********************/
    /* Simulating error in message */
    encoded_message[5] += 20;
    encoded_message[10] -= 52;
    encoded_message[20] -= 52;


    compute_poly_syndromes(encoded_message);

    /* Check if poly_syndromes is empty */
    if (*((uint32_t *)(poly_syndromes)) == 0)
    {
        printf("Message is not corrupted\n");
    }
    else
    {

        printf("\nCorrupted message:\n");
        print_hex_n(encoded_message, DATA_SIZE);
        printf("\n");

        /* Compute error locator polynomial */
        compute_poly_err_locator();
        
        /* Invert order of coefficients in poly_err_locator for easier calculations later */
        for (uint8_t i = 0; i < N_SYMBOLS / 2; ++i)
        {
            poly_err_locator[i] ^= poly_err_locator[N_SYMBOLS - i - 1];
            poly_err_locator[N_SYMBOLS - i - 1] ^= poly_err_locator[i];
            poly_err_locator[i] ^= poly_err_locator[N_SYMBOLS - i - 1];
        }

        /* Find indices of errors */
        compute_poly_err_evaluator();
        if (poly_err_evaluator_size == ERR_CODE)
        {
            printf("\nFound too many errors, message unrecoverable\n");
        }
        else
        {
            /* Find the error/deviation from original message */
            if (compute_poly_corruption() == ERR_CODE)
            {
                printf("\nError locators error, message unrecoverable\n");
            }
            else
            {
                /* Subtract previously found deviation */
                encoded_message[poly_err_evaluator[0]] ^= poly_corruption[0];
                encoded_message[poly_err_evaluator[1]] ^= poly_corruption[1];

                printf("\nRecovered message:\n");
                print_hex_n(encoded_message, DATA_SIZE);
                printf("\n");
            }
        }
        
    }

    return 0;
}





/* LOCAL FUNCTIONS DEFINITIONS BEGIN    */

/* Function for computing 2^pow */
static uint8_t gf_pow2(uint8_t pow)
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

/* Function for calculating inverse through brute-force*/
static uint8_t gf_inv(uint8_t x)
{
    /* Ignore edge case */
    if (x == 0)     return 0;

    uint8_t y = 1;
    while (gf_mult(x, y) != 1) ++y;

    return y;
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
static void compute_poly_syndromes(const uint8_t *encoded_message)
{
    for (uint8_t i = 0; i < N_SYMBOLS; ++i)
    {
        poly_syndromes[i] = gf_poly_evaluate(encoded_message, DATA_SIZE + N_SYMBOLS, 1 << i);       /* because our N_SYMBOLS is 4 we can treat 2^i in GF(2^8) as simple 1 << i */
    }
}

/* Function for computing error locator polynomial (Berlekamp-Massey algorithm) */
static void compute_poly_err_locator()
{
    uint8_t curr_poly_err_locator[N_SYMBOLS] = { 1, 0, 0, 0 };
    uint8_t old_poly_err_locator[N_SYMBOLS]  = { 1, 0, 0, 0 };
    uint8_t curr_size = 1;
    uint8_t old_size = 1;
    uint8_t discrepancy = 0;

    for (uint8_t i = 0; i < N_SYMBOLS; ++i)
    {
        discrepancy = poly_syndromes[i];

        for (uint8_t j = 1; j < curr_size; ++j)
        {
            discrepancy ^= gf_mult(curr_poly_err_locator[curr_size-j-1], poly_syndromes[i-j]);
        }
        ++old_size;

        if (discrepancy != 0)
        {
            uint8_t temp_poly_err_locator[N_SYMBOLS];

            if (old_size > curr_size)
            {
                gf_poly_scale(temp_poly_err_locator, old_poly_err_locator, N_SYMBOLS, discrepancy);
                gf_poly_scale(old_poly_err_locator, curr_poly_err_locator, N_SYMBOLS, gf_inv(discrepancy));

                /* swap old_size and curr_size */
                old_size ^= curr_size;
                curr_size ^= old_size;
                old_size ^= curr_size;

                memcpy(curr_poly_err_locator, temp_poly_err_locator, N_SYMBOLS);
            }

            gf_poly_scale(temp_poly_err_locator, old_poly_err_locator, N_SYMBOLS, discrepancy);

            uint8_t size_shift = curr_size - old_size;

            for (uint8_t iter = size_shift; iter < N_SYMBOLS; ++iter)
            {
                curr_poly_err_locator[iter] ^= temp_poly_err_locator[iter - size_shift];
            }
        }
    }

    memcpy(poly_err_locator, curr_poly_err_locator, N_SYMBOLS);
    poly_err_locator_size = curr_size;
}


static void compute_poly_err_evaluator ()
{
    uint8_t error_counter = 0;
    /* Brute force checking */
    for (uint8_t i = 0; i < DATA_SIZE + N_SYMBOLS; ++i)
    {
        /* Check every byte for condition to be root of polynomial error locator (the place where the erroc has happened) */
        uint8_t poly_value = gf_poly_evaluate(poly_err_locator, N_SYMBOLS, gf_pow2(i));
        if (poly_value == 0)
        {
            printf("Found error\n");
            if (2*error_counter >= N_SYMBOLS)
            {
                error_counter = ERR_CODE;
                /* Found too many errors */
                break;
            }
            poly_err_evaluator[error_counter++] = DATA_SIZE + N_SYMBOLS - 1 - i;
        }
    }

    poly_err_evaluator_size = error_counter;
}


static uint8_t compute_poly_corruption ()
{
    uint8_t ret_val = 0;
    if (poly_err_evaluator_size == 1)
    {
        uint8_t err_eval[2] = {
                                gf_mult(poly_syndromes[0], poly_err_locator[2]),
                                0
                              };

        uint8_t X_inv       = gf_pow2(256 - DATA_SIZE - N_SYMBOLS + poly_err_evaluator[0]);
        uint8_t X           = gf_inv(X_inv);
        uint8_t y           = gf_poly_evaluate(err_eval, 2, X_inv);

        poly_corruption[0]  = gf_mult(X, y);

    }
    else if (poly_err_evaluator_size == 2)
    {
        uint8_t err_eval[3] = { 
                                gf_mult(poly_syndromes[1], poly_err_locator[1]) ^ gf_mult(poly_syndromes[0], poly_err_locator[2]),
                                gf_mult(poly_syndromes[0], poly_err_locator[1]),
                                0
                              };

        uint8_t X_inv[2] = { 
                             gf_pow2(256 - DATA_SIZE - N_SYMBOLS + poly_err_evaluator[0]),
                             gf_pow2(256 - DATA_SIZE - N_SYMBOLS + poly_err_evaluator[1])
                           };
        uint8_t X[2] = {
                         gf_inv(X_inv[0]),
                         gf_inv(X_inv[1])
                       };
        uint8_t scale_adjustment;
        uint8_t y;

        /* index 0 */
        scale_adjustment = 0x01 ^ gf_mult(X_inv[0], X[1]);
        if (scale_adjustment == 0)
        {
            ret_val = ERR_CODE;
        }

        y = gf_poly_evaluate(err_eval, 3, X_inv[0]);
        y = gf_mult(X[0], y);
        poly_corruption[0] = gf_mult(y, gf_inv(scale_adjustment));

        /* index 1 */
        scale_adjustment = 0x01 ^ gf_mult(X_inv[1], X[0]);
        if (scale_adjustment == 0)
        {
            ret_val = ERR_CODE;
        }

        y = gf_poly_evaluate(err_eval, 3, X_inv[1]);
        y = gf_mult(X[1], y);
        poly_corruption[1] = gf_mult(y, gf_inv(scale_adjustment));
    }
    else
    {
        ret_val = ERR_CODE;
    }
    return ret_val;
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

