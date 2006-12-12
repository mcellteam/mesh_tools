#ifndef TILE_HASH_H
#define TILE_HASH_H

void free_tile_hash_table();

void init_tile_hash_table(int size);

static int calculate_tile_hash_value(int addr1, int addr2);

int is_in_tile_hash_table(int addr1, int addr2);

void add_an_entrance_to_tile_hash(int addr1, int addr2);

#endif
