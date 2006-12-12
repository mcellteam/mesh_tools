#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "common.h"
#include "myutil.h"
#include "tile_hash.h"

static int Hash_table_size;
static short *Hash_addr1;
static short *Hash_addr2;
static int Item_num;
extern int STATIC_CHQ;

void my_clear_tile_hash()
{
	Hash_table_size =0;
	Hash_addr1 = NULL;
	Hash_addr2 = NULL;
	Item_num =0;
}

void free_tile_hash_table()
{
    if(Hash_addr1)free(Hash_addr1);
    if(Hash_addr2)free(Hash_addr2);
    Hash_addr1=NULL;
    Hash_addr2=NULL;
}

void init_tile_hash_table(int size)
{
    int i;
    Hash_table_size=size*2+1;
    if(Hash_addr1)free(Hash_addr1);
    if(Hash_addr2)free(Hash_addr2);

    Hash_addr1=(short*)mymalloc(Hash_table_size*sizeof(short));
    Hash_addr2=(short*)mymalloc(Hash_table_size*sizeof(short));
    
	for(i=Hash_table_size-1; i>=0; i--)
	{
		Hash_addr1[i]=-1;
		Hash_addr2[i]=-1;
    }
    Item_num=0;
}

static int calculate_tile_hash_value(int addr1, int addr2)
{
    int val;
    
	if(addr1<0 || addr2<0)
	{
		printf("*** Error calculate_tile_hash_value(),addr1=%d addr2=%d\n",
			addr1,addr2);
		return -1;
    }
    
	val=(addr1*10001+addr2*7)%Hash_table_size;
    return val;
}


int is_in_tile_hash_table(int addr1, int addr2)
{
    int hash_val;
    hash_val=calculate_tile_hash_value(addr1, addr2);
    
	while(Hash_addr1[hash_val]>=0 && 
		!((Hash_addr1[hash_val]==addr1 && Hash_addr2[hash_val]==addr2))) 
		if(++hash_val>=Hash_table_size)hash_val=0;
		if(Hash_addr1[hash_val]<0)return 0;
		else return 1;
}


void add_an_entrance_to_tile_hash(int addr1, int addr2)
{ 
    int hash_val;
    hash_val=calculate_tile_hash_value(addr1, addr2);

    while(Hash_addr1[hash_val]>=0 && 
		!((Hash_addr1[hash_val]==addr1 && Hash_addr2[hash_val]==addr2))) 
		if(++hash_val>=Hash_table_size)hash_val=0;
		if(Hash_addr1[hash_val]>=0)
		{
		/*
		fprintf(stderr,"Err. add_an_entrance_to_tile_hash(), %d %d is in\n",
		addr1,addr2);
			*/
			return ;
		}
		
		if(Item_num>Hash_table_size-3)
		{
			printf("*** Error, add_..tile_hash(), over%d >= %d\n",
				Item_num,Hash_table_size-3);
			return;
		}
		
		Hash_addr1[hash_val]=addr1;
		Hash_addr2[hash_val]=addr2;
		Item_num++;
}
