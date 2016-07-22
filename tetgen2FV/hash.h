#ifndef HASH_H
#define HASH_H

#include "tetgen2FV.h"

typedef unsigned int ub4;   /* unsigned 4-byte quantities */
typedef unsigned char ub1;        /* unsigned 1-byte quantities */

struct hash_table {
        struct hash_table *next;  /**< next item in hash table*/
        char *key;               /**< item key*/
        void *contents;             /**< ptr to item contents*/
};

ub4 jerkins_hash(ub1 *key);
unsigned long hash(char *key);
struct hash_table *retrieve_key(char *key,
                                ub4 mask,
                                struct hash_table **hashtab);
struct hash_table *store_key(char *key,
                             ub4 mask,
                             struct hash_table **hashtab);
struct hash_table **init_hashtab(ub4 size);


#endif
